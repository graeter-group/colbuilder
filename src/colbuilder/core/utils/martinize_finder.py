"""
This module provides utilities for managing and integrating the Martinize2 tool and custom force fields
into the ColBuilder pipeline. It includes functionality for locating the Martinize2 executable,
installing custom force fields and mappings, and constructing Conda-based commands for execution.

Key Features:
--------------
1. **Conda Environment Detection**:
   - Identify the active Conda environment and retrieve its name and path.
   - Support for environments with Martinize2 installed.

2. **Martinize2 Executable Finder**:
   - Locate the Martinize2 executable in the system PATH or active Conda environment.
   - Provide fallback to Conda-based execution if the executable is not directly found.

3. **Custom Force Field Installation**:
   - Install custom force fields (e.g., amber99, martini300C) and mappings into Vermouth's directories.
   - Handle non-standard files like `modifications.mapping` and `selectors.py`.

4. **Conda Command Construction**:
   - Generate Conda-based commands for running Martinize2 with the active environment.

Usage:
------
This module is designed to be used as part of the ColBuilder pipeline to integrate Martinize2 and
custom force fields. It can also be used independently to manage Martinize2-related configurations.

Example:
--------
```python
from colbuilder.core.utils.martinize_finder import (
    find_martinize2_executable,
    find_and_install_custom_force_field,
    get_conda_command_with_path
)

# Find Martinize2 executable
executable, use_conda, conda_env = find_martinize2_executable()
print(f"Executable: {executable}, Use Conda: {use_conda}, Conda Env: {conda_env}")

# Install custom force fields
source_dir = "/path/to/force_fields"
success = find_and_install_custom_force_field(source_dir)
print(f"Force field installation successful: {success}")

# Construct Conda command
command = get_conda_command_with_path("martinize2", "-f input.pdb -o output.itp")
print(f"Conda command: {command}")
```
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
import subprocess
import shutil
from pathlib import Path
import json
from typing import Tuple, Optional, Union

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


def get_active_conda_env() -> Tuple[Optional[str], Optional[str]]:
    """
    Get information about the currently active conda environment.

    Returns
    -------
    Tuple[Optional[str], Optional[str]]
        Tuple containing (env_name, env_path) - name and filesystem path of the active environment
    """
    try:
        env_name = os.environ.get("CONDA_DEFAULT_ENV")
        env_path = os.environ.get("CONDA_PREFIX")

        if env_name and env_path:
            return env_name, env_path

        result = subprocess.run(
            ["conda", "info", "--envs", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        env_data = json.loads(result.stdout)

        for env in env_data.get("envs", []):
            if os.path.exists(env) and env == os.environ.get("CONDA_PREFIX"):
                env_name = os.path.basename(env)
                return env_name, env

        LOG.warning("Could not determine active conda environment")
        return None, None
    except Exception as e:
        LOG.error(f"Error determining active conda environment: {e}")
        return None, None


def find_and_install_custom_force_field(source_dir: Union[str, Path]) -> bool:
    """
    Find vermouth's directories and install custom force fields and mappings.

    Parameters
    ----------
    source_dir : Union[str, Path]
        Path to directory containing force field files and directories to copy

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        try:
            import vermouth

            vermouth_dir = Path(vermouth.__file__).parent
            LOG.debug(f"Found vermouth at: {vermouth_dir}")
        except ImportError:
            LOG.error(
                "Could not import vermouth module - is martinize2 installed correctly?"
            )
            return False

        # Set source path
        source_dir = Path(source_dir)
        if not source_dir.exists():
            LOG.error(f"Source directory does not exist: {source_dir}")
            return False

        # Define target directories
        force_fields_dir = vermouth_dir / "data" / "force_fields"
        mappings_dir = vermouth_dir / "data" / "mappings"

        # Ensure target directories exist
        force_fields_dir.mkdir(exist_ok=True, parents=True)
        mappings_dir.mkdir(exist_ok=True, parents=True)

        # Copy amber99 force field directory
        amber_source = source_dir / "amber99"
        amber_dest = force_fields_dir / "amber99"
        if amber_source.exists() and amber_source.is_dir():
            if amber_dest.exists():
                shutil.rmtree(amber_dest)
            shutil.copytree(amber_source, amber_dest)
            LOG.debug(f"Copied amber99 force field from {amber_source} to {amber_dest}")
        else:
            LOG.warning(f"Amber99 force field directory not found at {amber_source}")

        # Copy martini300C force field directory
        martini_source = source_dir / "martini300C-ff"
        martini_dest = force_fields_dir / "martini300C"
        if martini_source.exists() and martini_source.is_dir():
            if martini_dest.exists():
                shutil.rmtree(martini_dest)
            shutil.copytree(martini_source, martini_dest)
            LOG.debug(
                f"Copied martini300C force field from {martini_source} to {martini_dest}"
            )
        else:
            LOG.warning(
                f"Martini300C force field directory not found at {martini_source}"
            )

        # Copy martini300C-mapping directory to correct location
        mapping_source = source_dir / "martini300C-mapping"
        mapping_dest = mappings_dir / "martini300C"
        if mapping_source.exists() and mapping_source.is_dir():
            if mapping_dest.exists():
                shutil.rmtree(mapping_dest)
            shutil.copytree(mapping_source, mapping_dest)
            LOG.debug(
                f"Copied martini300C mappings from {mapping_source} to {mapping_dest}"
            )
        else:
            LOG.warning(f"Martini300C mapping directory not found at {mapping_source}")

        # Copy modifications.mapping file
        mods_source = source_dir / "modifications.mapping"
        mods_dest = mappings_dir / "modifications.mapping"
        if mods_source.exists():
            shutil.copy2(mods_source, mods_dest)
            LOG.debug(f"Copied modifications.mapping from {mods_source} to {mods_dest}")
        else:
            LOG.warning(f"Modifications mapping file not found at {mods_source}")

        # Copy selectors.py file
        selectors_source = source_dir / "selectors.py"
        selectors_dest = vermouth_dir / "selectors.py"
        if selectors_source.exists():
            shutil.copy2(selectors_source, selectors_dest)
            LOG.debug(
                f"Copied selectors.py from {selectors_source} to {selectors_dest}"
            )
        else:
            LOG.warning(f"Selectors.py file not found at {selectors_source}")

        LOG.debug("Successfully installed custom force fields and mappings")
        return True

    except Exception as e:
        LOG.error(f"Failed to install custom force fields: {e}")
        LOG.debug(f"Error details: {str(e)}", exc_info=True)
        return False


def find_martinize2_executable() -> Tuple[str, bool, Optional[str]]:
    """
    Find the martinize2 executable in a way that works across different environments.

    Returns
    -------
    Tuple[str, bool, Optional[str]]
        Tuple containing:
        - executable_path: Path to martinize2 or 'martinize2' if using conda run
        - use_conda: Whether to use conda run
        - conda_env: Name of conda environment (if use_conda is True)
    """
    LOG.debug("Looking for martinize2 executable")

    direct_path = shutil.which("martinize2")
    if direct_path:
        LOG.debug(f"Found martinize2 in PATH: {direct_path}")
        return direct_path, False, None

    env_name, env_path = get_active_conda_env()
    if env_path:
        for subdir in ["bin", "Scripts"]:
            exe_name = "martinize2.exe" if subdir == "Scripts" else "martinize2"
            possible_path = Path(env_path) / subdir / exe_name
            if possible_path.exists():
                LOG.info(f"Found martinize2 in current conda env: {possible_path}")
                return str(possible_path), False, None

        LOG.info(f"Using conda run with environment: {env_name}")
        return "martinize2", True, env_name

    LOG.warning(
        "Could not find martinize2. Defaulting to conda run with base environment"
    )
    return "martinize2", True, "base"


def get_conda_command_with_path(command: str, args: str) -> str:
    """
    Get a conda run command using the active environment path to avoid env resolution issues.

    Parameters
    ----------
    command : str
        The command to run
    args : str
        Command line arguments

    Returns
    -------
    str
        Full conda run command with path-based environment
    """
    env_name, env_path = get_active_conda_env()

    if env_path:
        LOG.debug(f"Using conda run with environment path: {env_path}")
        return f"conda run -p {env_path} {command} {args}"
    else:
        env_name = env_name or "base"
        LOG.warning(
            f"Could not find conda environment path, using environment name: {env_name}"
        )
        return f"conda run -n {env_name} {command} {args}"
