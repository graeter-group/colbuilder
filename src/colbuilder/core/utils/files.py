"""
File management utilities for the ColBuilder system.

This module provides utilities for handling file operations, resource management, and progress
tracking throughout the ColBuilder pipeline. It includes context managers for resource tracking,
output suppression, PDB file manipulation, and a centralized `FileManager` class for consistent
file handling.
"""

import contextlib
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
from typing import (
    Generator,
    Optional,
    Protocol,
    TypeVar,
    Any,
    Set,
    Dict,
    List,
    Callable,
)
import functools
import io
import time
import shutil
import os
from dataclasses import dataclass

from colbuilder.core.utils.exceptions import SequenceGenerationError, SystemError
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)

OperationType = TypeVar("OperationType")


@dataclass
class OperationContext:
    """
    Base context for all operations.

    Attributes:
        config: Configuration for the operation
        working_dir: Working directory path
    """

    config: ColbuilderConfig
    working_dir: Path


class Operation(Protocol[OperationType]):
    """
    Base protocol for all operations.

    This protocol defines the common interface that all operations must implement,
    allowing for consistent execution patterns throughout the system.
    """

    async def execute(self, input_data: Optional[Any] = None) -> OperationType:
        """
        Execute the operation.

        Args:
            input_data: Optional input data for the operation

        Returns:
            The operation result of type OperationType
        """
        pass


@contextmanager
def managed_resources(resource_name=None):
    """
    Decorator for methods that need temporary resource management.

    This decorator ensures that a method properly initializes temporary
    resources and cleans them up after execution, even if an exception occurs.

    Args:
        resource_name: Optional name for the resource category

    Returns:
        Decorated method
    """

    def decorator(func):
        @functools.wraps(func)
        async def wrapper(self, *args, **kwargs):
            resource_context = resource_name or func.__name__
            LOG.debug(f"Managing resources for {resource_context}")

            original_dir = None
            temp_dir = None

            try:
                # Store original directory
                original_dir = Path.cwd()

                # Create a temporary directory if needed
                if hasattr(self, "file_manager") and self.file_manager:
                    temp_dir = self.file_manager.get_temp_path(
                        resource_context, create_dir=True
                    )
                    if hasattr(self, "temp_dir"):
                        self.temp_dir = temp_dir
                    LOG.debug(
                        f"Created temporary directory for {resource_context}: {temp_dir}"
                    )

                # Call the decorated method
                return await func(self, *args, **kwargs)

            except Exception as e:
                LOG.error(f"Error in {resource_context}: {str(e)}")
                raise

            finally:
                # Return to original directory
                if original_dir:
                    os.chdir(original_dir)
                    LOG.debug(f"Returned to original directory: {original_dir}")

                # Don't clean up if debug mode is enabled
                if (
                    hasattr(self, "config")
                    and hasattr(self.config, "debug")
                    and self.config.debug
                ):
                    LOG.info(
                        f"Debug mode enabled, preserving temporary directory: {temp_dir}"
                    )
                    # Create a marker file in temp dir with info
                    if temp_dir and temp_dir.exists():
                        try:
                            with open(temp_dir / "_DEBUG_INFO.txt", "w") as f:
                                f.write(f"Debug information for {resource_context}\n")
                                f.write(f"Created at: {time.ctime()}\n")
                                f.write(f"Function: {func.__name__}\n")
                                if hasattr(self, "config"):
                                    f.write(f"Configuration:\n")
                                    for key, value in vars(self.config).items():
                                        if not key.startswith("_"):
                                            f.write(f"  {key}: {value}\n")
                        except Exception as e:
                            LOG.warning(f"Failed to create debug info file: {e}")
                else:
                    # Clean up unless debugging
                    if hasattr(self, "file_manager") and self.file_manager:
                        LOG.debug(f"Cleaning up resources for {resource_context}")
                        if resource_context in [
                            "geometry_operation",
                            "replacement",
                            "replace_generation",
                        ]:
                            LOG.debug(
                                f"Skipping cleanup of important directories: {resource_context}"
                            )
                        else:
                            self.file_manager.cleanup()

        return wrapper

    return decorator


@contextmanager
def suppress_output() -> Generator[None, None, None]:
    """
    Context manager to suppress stdout and stderr output.

    This utility captures and discards all standard output and error streams
    during its execution scope, which is useful when calling noisy external
    libraries or tools where their console output is not relevant.

    Yields:
        None

    Raises:
        SystemError: If there's an error managing output streams
    """
    try:
        with io.StringIO() as stdout_buf, io.StringIO() as stderr_buf:
            with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                yield
    except Exception as e:
        raise SystemError(
            message="Failed to manage output streams",
            original_error=e,
            error_code="SYS_ERR_001",
            context={"action": "suppress_output"},
        )


class ProgressTracker:
    """
    Tracks progress of multi-step operations.

    This class provides a simple way to track and report progress through
    a sequence of steps, maintaining state about how far along the process is.
    """

    def __init__(self, total_steps: int) -> None:
        """
        Initialize the progress tracker.

        Args:
            total_steps: Total number of steps in the operation
        """
        self.total_steps: int = total_steps
        self.current_step: int = 0

    def update(self, message: str) -> None:
        """
        Update progress with a step completion message.

        Increments the step counter and logs the progress message with
        the current step number and total steps.

        Args:
            message: Progress message to log
        """
        self.current_step += 1
        LOG.info(f"Step {self.current_step}/{self.total_steps}: {message}")


def update_pdb_header(pdb_file: Path, first_line: str = "") -> None:
    """
    Update PDB file header, removing duplicates and unnecessary lines.

    This function reads a PDB file, extracts CRYST1 and ATOM records,
    and writes them back with an optional custom first line. It helps
    standardize PDB files and ensure they have the correct header structure.

    Args:
        pdb_file: Path to PDB file
        first_line: Optional first line to add to the file

    Raises:
        SequenceGenerationError: If update fails
    """
    try:
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        cryst_lines = []
        atom_lines = []
        in_atom_section = False

        for line in lines:
            if line.startswith("CRYST1"):
                cryst_lines.append(line)
            elif line.startswith("ATOM"):
                in_atom_section = True
                atom_lines.append(line)
            elif in_atom_section:
                atom_lines.append(line)

        output_lines = []

        if first_line and not first_line.isspace() and len(first_line) > 0:
            output_lines.append(first_line + "\n")

        if cryst_lines:
            output_lines.append(cryst_lines[0])

        output_lines.extend(atom_lines)

        with open(pdb_file, "w") as f:
            f.writelines(output_lines)

    except Exception as e:
        LOG.error(f"Error updating PDB header: {str(e)}")
        raise SequenceGenerationError(
            "Failed to update PDB header",
            error_code="SEQ_ERR_004",
            context={"pdb_file": str(pdb_file), "error": str(e)},
        )


class FileManager:
    """
    Centralized utility for file management in Colbuilder.

    This class provides consistent methods for managing temporary and permanent files,
    ensuring proper cleanup and respecting debug settings. It handles file paths,
    temporary directories, copying files, and directory creation with consistent
    logging and error handling.

    Attributes:
        config: Colbuilder configuration
        temp_files: Set of temporary files to track
        temp_dirs: Set of temporary directories to track
    """

    TEMP_DIR = ".tmp"
    TEMP_SEQUENCE_DIR = "sequence_gen"
    TEMP_GEOMETRY_DIR = "geometry_gen"
    TEMP_REPLACEMENT_DIR = "replace_crosslinks"
    TEMP_MIXING_DIR = "mixing_crosslinks"
    TEMP_TOPOLOGY_DIR = "topology_gen"

    def __init__(self, config: ColbuilderConfig) -> None:
        """
        Initialize with configuration settings.

        Args:
            config: Colbuilder configuration
        """
        self.config: ColbuilderConfig = config
        self.project_root = config.PROJECT_ROOT  # Use the project root from config
        self.working_dir = Path(config.working_directory).resolve()
        self.temp_base = self.working_dir / self.TEMP_DIR
        self.geometry_dir = self.temp_base / self.TEMP_GEOMETRY_DIR
        self.replacement_dir = self.temp_base / self.TEMP_REPLACEMENT_DIR
        self.mixing_dir = self.temp_base / self.TEMP_MIXING_DIR
        self.topology_dir = self.temp_base / self.TEMP_TOPOLOGY_DIR
        self.sequence_dir = self.temp_base / self.TEMP_SEQUENCE_DIR
        self.temp_files: Set[Path] = set()
        self.temp_dirs: Set[Path] = set()

        # Create and track standard directories
        self._ensure_standard_dirs()

    def _ensure_standard_dirs(self) -> None:
        """Create and track standard directories."""
        self.temp_base.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.temp_base)

        self.geometry_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.geometry_dir)

        self.replacement_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.replacement_dir)

        self.mixing_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.mixing_dir)

        self.topology_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.topology_dir)

        self.sequence_dir.mkdir(exist_ok=True, parents=True)
        self.temp_dirs.add(self.sequence_dir)

    def get_temp_path(
        self, basename: str, suffix: Optional[str] = None, create_dir: bool = False
    ) -> Path:
        """
        Get a path for a temporary file or directory.

        Creates a path for a temporary file or directory, optionally creating
        the directory structure. The path is tracked for later cleanup.

        Args:
            basename: Base name for the file or directory
            suffix: Optional suffix to add (e.g., file extension)
            create_dir: Whether to create a directory at this path

        Returns:
            Path object for the temporary file or directory
        """
        LOG.debug(f"Getting temp path for basename: {basename}")

        # Sanitize the basename to prevent nested .tmp directories
        clean_basename = str(basename)
        if clean_basename.startswith(".tmp/"):
            LOG.debug(f"Removing leading .tmp/ from path: {clean_basename}")
            clean_basename = clean_basename[5:]

        if ".tmp/.tmp" in clean_basename:
            LOG.debug(f"Removing nested .tmp from path: {clean_basename}")
            clean_basename = clean_basename.replace(".tmp/.tmp", ".tmp")

        # Handle standard directory prefixes
        if clean_basename.startswith(self.TEMP_GEOMETRY_DIR):
            LOG.debug(f"Using standard geometry directory for: {clean_basename}")
            path_parts = Path(clean_basename).parts[1:]  # Skip the geometry_gen part
            if path_parts:
                result_path = self.geometry_dir.joinpath(*path_parts)
            else:
                result_path = self.geometry_dir

        elif clean_basename.startswith(self.TEMP_REPLACEMENT_DIR):
            LOG.debug(f"Using standard replacement directory for: {clean_basename}")
            path_parts = Path(clean_basename).parts[
                1:
            ]  # Skip the replace_crosslinks part
            if path_parts:
                result_path = self.replacement_dir.joinpath(*path_parts)
            else:
                result_path = self.replacement_dir

        elif clean_basename.startswith(self.TEMP_MIXING_DIR):
            LOG.debug(f"Using standard mixing directory for: {clean_basename}")
            path_parts = Path(clean_basename).parts[
                1:
            ]  # Skip the mixing_crosslinks part
            if path_parts:
                result_path = self.mixing_dir.joinpath(*path_parts)
            else:
                result_path = self.mixing_dir

        elif clean_basename.startswith(self.TEMP_TOPOLOGY_DIR):
            LOG.debug(f"Using standard topology directory for: {clean_basename}")
            path_parts = Path(clean_basename).parts[1:]  # Skip the topology_gen part
            if path_parts:
                result_path = self.topology_dir.joinpath(*path_parts)
            else:
                result_path = self.topology_dir

        elif clean_basename.startswith(self.TEMP_SEQUENCE_DIR):
            path_parts = Path(clean_basename).parts[1:]
            if path_parts:
                result_path = self.sequence_dir.joinpath(*path_parts)
            else:
                result_path = self.sequence_dir

        else:
            # Handle other paths - determine if it's a nested path or simple filename
            path_obj = Path(clean_basename)
            if len(path_obj.parts) > 1:
                # It's a nested path, create in temp base
                result_path = self.temp_base.joinpath(*path_obj.parts)
                LOG.info(f"Using nested path under temp base: {result_path}")
            else:
                # It's a simple filename
                result_path = self.temp_base / clean_basename
                LOG.info(f"Using simple path under temp base: {result_path}")

        # Add suffix if provided
        if suffix:
            result_path = Path(str(result_path) + suffix)

        # Create directory if requested
        if create_dir:
            result_path.mkdir(exist_ok=True, parents=True)
            self.temp_dirs.add(result_path)
        else:
            # Make sure parent exists
            result_path.parent.mkdir(exist_ok=True, parents=True)
            self.temp_files.add(result_path)
            LOG.info(f"Prepared file path: {result_path}")

        # Final sanity check
        if ".tmp/.tmp" in str(result_path):
            LOG.warning(f"Found nested .tmp in final path: {result_path}")
            fixed_path = Path(str(result_path).replace(".tmp/.tmp", ".tmp"))
            LOG.info(f"Fixed path to: {fixed_path}")

            # Update tracking if needed
            if create_dir:
                self.temp_dirs.remove(result_path)
                self.temp_dirs.add(fixed_path)
                fixed_path.mkdir(exist_ok=True, parents=True)
            else:
                self.temp_files.remove(result_path)
                self.temp_files.add(fixed_path)
                fixed_path.parent.mkdir(exist_ok=True, parents=True)

            return fixed_path

        return result_path

    def get_output_path(
        self, basename: str, suffix: Optional[str] = None, mkdir: bool = False
    ) -> Path:
        """
        Get a path for a permanent output file or directory.

        Creates a path for an output file or directory in the working directory,
        optionally creating the directory structure.

        Args:
            basename: Base name for the file or directory
            suffix: Optional suffix to add (e.g., file extension)
            mkdir: Whether to create a directory at this path

        Returns:
            Path object for the output file or directory
        """
        full_name = f"{basename}{suffix if suffix else ''}"
        full_path = self.config.working_directory / full_name

        if mkdir:
            full_path.mkdir(exist_ok=True, parents=True)

        return full_path

    def ensure_geometry_dir(self) -> Path:
        """
        Ensure geometry directory exists and return its path.

        Returns:
            Path: The path to the geometry directory
        """
        if not self.geometry_dir.exists():
            self.geometry_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dirs.add(self.geometry_dir)
            LOG.info(f"Created geometry directory: {self.geometry_dir}")

        return self.geometry_dir

    def ensure_replacement_dir(self) -> Path:
        """
        Ensure replacement directory exists and return its path.

        Returns:
            Path: The path to the replacement directory
        """
        if not self.replacement_dir.exists():
            self.replacement_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dirs.add(self.replacement_dir)
            LOG.info(f"Created replacement directory: {self.replacement_dir}")

        return self.replacement_dir

    def ensure_mixing_dir(self) -> Path:
        """
        Ensure mixing directory exists and return its path.

        Returns:
            Path: The path to the mixing directory
        """
        if not self.mixing_dir.exists():
            self.mixing_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dirs.add(self.mixing_dir)
            LOG.info(f"Created mixing directory: {self.mixing_dir}")

        return self.mixing_dir

    def cleanup(self, category: Optional[str] = None, force: bool = False) -> None:
        """
        Clean up temporary files and directories.

        Removes all tracked temporary files and directories, unless
        in debug mode and not forced. Always preserves the .tmp directory.

        Args:
            category: Optional category of resources to clean up
            force: Whether to clean up even in debug mode
        """
        # Skip cleanup in debug mode
        if self.config.debug and not force:
            LOG.info("Debug mode enabled, skipping cleanup")

            # Create a marker file in temp directories with info
            for dir_path in self.temp_dirs:
                try:
                    if dir_path.exists():
                        marker_file = dir_path / "_DEBUG_INFO.txt"
                        with open(marker_file, "w") as f:
                            f.write(f"Debug information for temp directory\n")
                            f.write(f"Created at: {time.ctime()}\n")
                            f.write(
                                f"Contains temporary files and directories for debugging\n"
                            )
                            f.write(f"Configuration:\n")
                            for key, value in vars(self.config).items():
                                if not key.startswith("_"):
                                    f.write(f"  {key}: {value}\n")
                except Exception as e:
                    LOG.warning(f"Failed to create debug info file in {dir_path}: {e}")
            return

        # Normal cleanup (when debug is False)
        for file_path in self.temp_files:
            try:
                if file_path.exists():
                    file_path.unlink()
                    LOG.debug(f"Removed temporary file: {file_path}")
            except Exception as e:
                LOG.warning(f"Failed to remove temporary file {file_path}: {str(e)}")

        for dir_path in sorted(self.temp_dirs, key=lambda p: -len(str(p))):
            try:
                if dir_path.exists():
                    # ALWAYS preserve the main directories
                    if dir_path in [
                        self.temp_base,
                        self.geometry_dir,
                        self.replacement_dir,
                        self.mixing_dir,
                        self.topology_dir,
                        self.sequence_dir,
                    ]:
                        continue
                    if dir_path.name == "replace_manual":
                        continue

                    # Skip certain directories in debug mode
                    if self.config.debug and any(
                        x in str(dir_path) for x in ["geometry_gen", "replacement", "D"]
                    ):
                        LOG.info(f"Preserving directory in debug mode: {dir_path}")
                        continue

                    shutil.rmtree(dir_path)
                    LOG.debug(f"Removed temporary directory: {dir_path}")
            except Exception as e:
                LOG.warning(
                    f"Failed to remove temporary directory {dir_path}: {str(e)}"
                )

        # Don't clear temp_files and temp_dirs lists to keep tracking them
        preserved_files = {f for f in self.temp_files if f.exists()}
        preserved_dirs = {d for d in self.temp_dirs if d.exists()}

        # Only remove tracked entries for items that were actually deleted
        self.temp_files = preserved_files
        self.temp_dirs = preserved_dirs

    def copy_to_output(self, source: Path, dest_name: Optional[str] = None) -> Path:
        """
        Copy a file to the output directory.

        Copies a source file to the working directory with the given
        destination name, or the original name if not specified.

        Args:
            source: Source file path
            dest_name: Optional destination name (uses source name if not provided)

        Returns:
            Path to the destination file

        Raises:
            Various exceptions from shutil.copy2 if copy fails
        """
        dest_name = dest_name or source.name
        dest_path = self.config.working_directory / dest_name

        try:
            shutil.copy2(source, dest_path)
            LOG.debug(f"Copied {source} to {dest_path}")
            return dest_path
        except Exception as e:
            LOG.error(f"Failed to copy {source} to {dest_path}: {str(e)}")
            raise

    def ensure_dir(self, dirname: str) -> Path:
        """
        Ensure a directory exists in the working directory.

        Creates the directory if it doesn't already exist.

        Args:
            dirname: Name of the directory

        Returns:
            Path to the directory
        """
        dir_path = self.config.working_directory / dirname
        dir_path.mkdir(exist_ok=True, parents=True)
        return dir_path

    @contextmanager
    def temp_file_context(
        self, basename: str, suffix: Optional[str] = None
    ) -> Generator[Path, None, None]:
        """
        Context manager for temporary file operations.

        Creates a temporary file path, yields it for use, and ensures
        it gets tracked for cleanup later.

        Args:
            basename: Base name for the file
            suffix: Optional suffix/extension

        Yields:
            Path to the temporary file
        """
        temp_path = self.get_temp_path(basename, suffix)
        try:
            yield temp_path
        finally:
            self.temp_files.add(temp_path)

    @contextmanager
    def temp_dir_context(self, dirname: str) -> Generator[Path, None, None]:
        """
        Context manager for temporary directory operations.

        Creates a temporary directory, yields it for use, and ensures
        it gets tracked for cleanup later.

        Args:
            dirname: Name for the directory

        Yields:
            Path to the temporary directory
        """
        temp_dir = self.get_temp_path(dirname, create_dir=True)
        try:
            yield temp_dir
        finally:
            self.temp_dirs.add(temp_dir)

    def find_file(
        self, filename_or_path: str, search_paths: Optional[List[Path]] = None
    ) -> Optional[Path]:
        """
        Find a file by searching in multiple locations.

        Searches for a file in the provided search paths, standard Colbuilder locations,
        and relative to the project root.

        Args:
            filename_or_path: Name of the file or relative path to find
            search_paths: Optional list of paths to search (added to standard locations)

        Returns:
            Path to the file if found, None otherwise
        """
        filepath = Path(filename_or_path)

        if filepath.is_absolute() and filepath.exists():
            return filepath

        if search_paths is None:
            search_paths = []

        # Add standard search locations
        search_paths.extend(
            [
                self.working_dir,  # User's current working directory
                Path.cwd(),  # Python process working directory
                self.project_root,  # Project root from config
                self.project_root / "data",  # Project data directory
                self.project_root / "colbuilder" / "data",  # Package data directory
                self.config.DATA_DIR,  # Data directory from config
                self.config.HOMOLOGY_LIB_DIR,  # Homology lib directory from config
                self.config.FORCE_FIELD_DIR,  # Force field directory from config
            ]
        )

        # Filter out None values
        search_paths = [p for p in search_paths if p is not None]

        for base_path in search_paths:
            full_path = base_path / filepath
            if full_path.exists():
                return full_path

            if len(filepath.parts) > 1:
                filename_only = filepath.name
                name_only_path = base_path / filename_only
                if name_only_path.exists():
                    return name_only_path

        return None

    def copy_to_directory(
        self,
        source: Path,
        dest_dir: Optional[Path] = None,
        dest_name: Optional[str] = None,
    ) -> Path:
        """
        Copy a file to a specified directory.

        Args:
            source (Path): Source file path.
            dest_dir (Optional[Path]): Destination directory (defaults to the working directory).
            dest_name (Optional[str]): Optional destination name (uses source name if not provided).

        Returns:
            Path: Path to the copied file in the destination directory.

        Raises:
            FileNotFoundError: If the source file does not exist.
            Exception: If the copy operation fails.
        """
        if not source.exists():
            raise FileNotFoundError(f"Source file not found: {source}")

        dest_dir = dest_dir or self.config.working_directory
        dest_dir.mkdir(parents=True, exist_ok=True)

        dest_name = dest_name or source.name
        dest_path = dest_dir / dest_name

        try:
            shutil.copy2(source, dest_path)
            LOG.debug(f"Copied {source} to {dest_path}")
            return dest_path
        except Exception as e:
            LOG.error(f"Failed to copy {source} to {dest_path}: {str(e)}")
            raise

    def get_type_dir(self, model_type: str, parent_dir: Optional[Path] = None) -> Path:
        """
        Get or create the directory for a specific model type.

        Args:
            model_type (str): The type of the model (e.g., "D", "NC").
            parent_dir (Optional[Path]): Optional parent directory (defaults to geometry_dir).

        Returns:
            Path: Path to the directory for the specified model type.

        Raises:
            FileNotFoundError: If the directory cannot be created or accessed.
        """
        try:
            base_dir = parent_dir or self.geometry_dir

            type_dir = base_dir / model_type

            type_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dirs.add(type_dir)
            LOG.info(f"Using type directory: {type_dir}")

            return type_dir
        except Exception as e:
            raise FileNotFoundError(
                f"Failed to locate or create type directory for model type: {model_type}"
            ) from e

    def get_temp_dir(self, dirname: str) -> Path:
        """
        Get a temporary directory path and create the directory.

        This is a convenience wrapper around get_temp_path with create_dir=True.

        Args:
            dirname: Name for the directory

        Returns:
            Path to the created temporary directory
        """
        temp_dir = self.get_temp_path(dirname, create_dir=True)
        return temp_dir
