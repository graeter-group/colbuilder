"""
Colbuilder Crystal Builder Module

This module provides services for building and managing crystal systems,
including initialization, contact distance calculations, and system optimization.

Key Features:
    - Crystal structure initialization
    - Contact distance calculations
    - System optimization
    - Cap addition
    - PDB file generation
"""

from pathlib import Path
from typing import Tuple, Optional, List, Dict, Any, Union
import shutil
import os
from colorama import Fore, Style

from colbuilder.core.utils.exceptions import GeometryGenerationError
from colbuilder.core.utils.logger import setup_logger
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.files import FileManager

from .crystal import Crystal
from .system import System
from .model import Model
from .crystalcontacts import CrystalContacts
from .connect import Connect
from .caps import Caps
from .optimize import Optimizer
from .chimera import Chimera
from . import crosslink

LOG = setup_logger(__name__)


class CrystalBuilder:
    """
    Service for crystal system operations.

    This class encapsulates all operations related to building and managing
    crystal systems, including initialization, optimization, and file handling.
    """

    def __init__(self, file_manager: Optional[FileManager] = None):
        """Initialize the crystal builder."""
        self.steps = 7
        self.file_manager = file_manager

    def set_file_manager(self, file_manager: FileManager) -> None:
        """Set the file manager for consistent file handling."""
        self.file_manager = file_manager

    @timeit
    def matrixset_system(
        self, system: System, crystalcontacts_file: Union[str, Path]
    ) -> System:
        """
        Process and update system after matrix generation.

        Args:
            system: System to process
            crystalcontacts_file: Path to crystal contacts file

        Returns:
            System: Processed system

        Raises:
            FileNotFoundError: If required files are missing
            GeometryGenerationError: If processing fails
        """
        LOG.debug(
            f"Setting system after cutting fibril. crystalcontacts_file type: {type(crystalcontacts_file)}"
        )
        crystalcontacts_file = Path(crystalcontacts_file)

        if crystalcontacts_file.suffix == ".txt":
            crystalcontacts_file = crystalcontacts_file.with_suffix("")

        id_file = crystalcontacts_file.with_name(f"{crystalcontacts_file.name}_id.txt")
        LOG.debug(f"Looking for crystal contacts ID file at: {id_file}")

        if not id_file.exists():
            raise FileNotFoundError(f"Crystal contacts ID file not found: {id_file}")

        try:
            with open(id_file, "r") as f:
                contacts = [float(i.split(" ")[1]) for i in f.readlines()]
        except Exception as e:
            LOG.error(f"Error reading crystal contacts ID file: {str(e)}")
            raise

        models_before = len(list(system.get_models()))
        for model in list(system.get_models()):
            if model not in contacts:
                system.delete_model(model_id=model)
            elif system.get_model(model_id=model).connect is not None:
                connect_before = len(list(system.get_model(model_id=model).connect))
                for connect in list(system.get_model(model_id=model).connect):
                    if connect not in contacts:
                        system.get_model(model_id=model).delete_connect(
                            connect_id=connect
                        )
                connect_after = len(list(system.get_model(model_id=model).connect))
                if connect_before != connect_after:
                    pass

        models_after = len(list(system.get_models()))
        return system

    async def build(self, config: ColbuilderConfig) -> System:
        """Build a crystal system from configuration."""
        if not self.file_manager:
            LOG.info(f"Creating new FileManager for CrystalBuilder")
            self.file_manager = FileManager(config)

        try:
            LOG.info(f"Step 1/{self.steps} Initializing crystal")
            crystal = await self._initialize_crystal(config)

            LOG.info(f"Step 2/{self.steps} Building initial system")
            system, crystalcontacts, connect = await self._build_initial_system(
                crystal, config
            )

            LOG.info(f"Step 3/{self.steps} Writing crystal contacts")
            crystalcontacts.write_crystalcontacts(
                system=system, crystalcontacts_file=crystalcontacts.crystalcontacts_file
            )
            LOG.debug(
                f"Crystal contacts written to: {crystalcontacts.crystalcontacts_file}"
            )

            LOG.info(f"Step 4/{self.steps} Generating system matrix")
            LOG.info(
                f"{Fore.BLUE}Please wait, this may take some time ...{Style.RESET_ALL}"
            )

            pdb_path = Path(config.pdb_file).resolve() if config.pdb_file else None
            LOG.debug(f"Using PDB file at: {pdb_path}")

            try:
                chimera_scripts_dir = self.file_manager.find_file("chimera_scripts")
            except FileNotFoundError:
                chimera_scripts_dir = self.file_manager.ensure_dir("chimera_scripts")
                LOG.warning(
                    f"Chimera scripts directory not found, created at: {chimera_scripts_dir}"
                )

            chimera = Chimera(config, str(pdb_path))
            chimera.matrixset(
                pdb=str(pdb_path),
                crystalcontacts=crystalcontacts.crystalcontacts_file,
                system_size=system.get_size(),
                fibril_length=config.fibril_length,
            )

            system = self.matrixset_system(
                system=system, crystalcontacts_file=crystalcontacts.crystalcontacts_file
            )

            # Refresh crosslinks from the actually generated per-model PDBs in the geometry directory.
            geometry_dir = self.file_manager.ensure_geometry_dir()
            for model_id in system.get_models():
                model_obj = system.get_model(model_id=model_id)
                pdb_path = geometry_dir / f"{int(model_id)}.pdb"

                model_obj.crosslink = []
                model_obj.pdb_file = str(pdb_path)

                if not pdb_path.exists():
                    LOG.warning(
                        f"Cannot refresh crosslinks for model {model_id}: {pdb_path} missing"
                    )
                    continue

                try:
                    cls = crosslink.read_crosslink(pdb_file=pdb_path)
                    for c in cls:
                        c.model_id = model_obj.id
                    model_obj.crosslink = cls
                except Exception as e:
                    LOG.warning(
                        f"Failed to reload crosslinks for model {model_id} from {pdb_path}: {e}"
                    )

            # Recompute connectivity using refreshed crosslinks
            connect = Connect(system=system, connect_file=connect.connect_file)
            system_connect = connect.run_connect(system=system)
            for mid in system.get_models():
                model_obj = system.get_model(model_id=mid)
                if system_connect and mid in system_connect:
                    model_obj.add_connect(connect_id=mid, connect=system_connect[mid])
                else:
                    model_obj.connect = None            

            LOG.info(f"Step 5/{self.steps} Writing models' connectivity")
            connect.write_connect(system=system, connect_file=connect.connect_file)
            LOG.debug(f"Connect file written to: {connect.connect_file}")

            LOG.info(f"Step 6/{self.steps} Adding caps")
            model_type = system.get_model(model_id=0.0).type
            LOG.info(f"{Fore.BLUE}Model type: {model_type}{Style.RESET_ALL}")

            geometry_dir = self.file_manager.ensure_geometry_dir()
            caps_dir = geometry_dir / model_type
            caps_dir.mkdir(parents=True, exist_ok=True)
            LOG.debug(f"Using caps directory: {caps_dir}")

            has_crosslinks = self._needs_optimization(system)
            if has_crosslinks:
                LOG.debug(f"System has crosslinks, using crosslink capping")
                caps = Caps(system=system)
                LOG.debug(f"Writing capped PDB files to {geometry_dir}")
                for idx in system.get_models():
                    pdb_id = int(idx)
                    pdb_file = f"{pdb_id}.pdb"
                    pdb_path = Path(pdb_file)

                    if not pdb_path.exists():
                        LOG.warning(
                            f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}"
                        )
                        continue
                    caps.read_residues(pdb_id=pdb_id)

                    caps.add_caps(
                        pdb_id=pdb_id, crosslink_type=model_type, temp_dir=geometry_dir
                    )
            else:
                caps = Caps(system=system)

                LOG.debug(f"Writing capped PDB files to {geometry_dir}")
                for idx in system.get_models():
                    pdb_id = int(idx)
                    pdb_file = f"{pdb_id}.pdb"
                    pdb_path = Path(pdb_file)

                    if not pdb_path.exists():
                        LOG.warning(
                            f"PDB file {pdb_file} not found. Model IDs: {list(system.get_models())}"
                        )
                        continue

                    caps.read_residues(pdb_id=pdb_id)
                    caps.add_caps(
                        pdb_id=pdb_id, crosslink_type="NC", temp_dir=geometry_dir
                    )

            LOG.info(f"Step 7/{self.steps} Writing final structure")
            output_path = self.file_manager.get_output_path(config.output, ".pdb")

            system.write_pdb(
                pdb_out=output_path,
                fibril_length=config.fibril_length,
                cleanup=False,
                temp_dir=geometry_dir,  # Pass the dedicated caps directory to write_pdb
            )
            LOG.info(
                f"{Fore.BLUE}Final structure written successfully{Style.RESET_ALL}"
            )

            return system

        except Exception as e:
            LOG.error(f"Failed to build crystal system: {str(e)}", exc_info=True)
            raise GeometryGenerationError(
                message="Failed to build crystal system",
                original_error=e,
                error_code="GEO_ERR_001",
                context={"config": config.model_dump()},
            )

    async def _initialize_crystal(self, config: ColbuilderConfig) -> Crystal:
        """
        Initialize crystal structure from PDB file.

        Args:
            config: Configuration settings

        Returns:
            Crystal: Initialized crystal

        Raises:
            GeometryGenerationError: If initialization fails
        """
        try:
            if not config.pdb_file:
                raise GeometryGenerationError(
                    message="PDB file not specified",
                    error_code="GEO_ERR_005",
                    context={"config": config.model_dump()},
                )

            pdb_path = Path(config.pdb_file).resolve()
            if not pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"PDB file not found: {pdb_path}", error_code="GEO_ERR_005"
                )

            pdb_str = str(pdb_path)
            if pdb_str.endswith(".pdb"):
                pdb_str = pdb_str[:-4]

            geometry_dir = self.file_manager.ensure_geometry_dir()
            original_pdb = geometry_dir / f"{os.path.basename(pdb_str)}_original.pdb"

            try:
                shutil.copy(f"{pdb_str}.pdb", original_pdb)
            except (shutil.Error, IOError) as e:
                raise GeometryGenerationError(
                    message="Failed to create backup of original PDB file",
                    original_error=e,
                    error_code="GEO_ERR_005",
                    context={
                        "source": f"{pdb_str}.pdb",
                        "destination": str(original_pdb),
                    },
                )

            crystal = Crystal(pdb_str)
            crystal.translate_crystal(pdb=pdb_str, translate=[0, 0, 4000])
            return crystal

        except Exception as e:
            LOG.error(f"Failed to initialize crystal structure: {str(e)}")
            raise GeometryGenerationError(
                message="Failed to initialize crystal structure",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"pdb_file": str(pdb_path) if "pdb_path" in locals() else None},
            )

    async def _build_initial_system(
        self, crystal: Crystal, config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        """
        Build initial system from crystal structure.

        Args:
            crystal: Crystal structure
            config: Configuration settings

        Returns:
            Tuple containing:
                - System: Built system
                - CrystalContacts: Crystal contacts information
                - Connect: System connections

        Raises:
            GeometryGenerationError: If system building fails
        """
        if config.contact_distance and not config.crystalcontacts_file:
            LOG.debug(
                f"Building system from contact distance: {config.contact_distance}"
            )
            return await self._build_from_contact_distance(crystal, config)
        elif not config.contact_distance and config.crystalcontacts_file:
            LOG.debug(
                f"Building system from crystal contacts file: {config.crystalcontacts_file}"
            )
            return await self._build_from_crystal_contacts(crystal, config)
        else:
            raise GeometryGenerationError(
                message="Invalid configuration for system building",
                error_code="GEO_ERR_001",
                context={
                    "contact_distance": config.contact_distance,
                    "crystalcontacts_file": (
                        str(config.crystalcontacts_file)
                        if config.crystalcontacts_file
                        else None
                    ),
                },
            )

    async def _build_from_contact_distance(
        self, crystal: Crystal, config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        try:
            if not self.file_manager:
                raise GeometryGenerationError(
                    message="FileManager not initialized in CrystalBuilder",
                    error_code="GEO_ERR_002",
                )

            pdb_path = Path(config.pdb_file).resolve() if config.pdb_file else None

            if not pdb_path or not pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"PDB file not found: {pdb_path}", error_code="GEO_ERR_005"
                )

            geometry_dir = self.file_manager.ensure_geometry_dir()
            LOG.debug(f"Using geometry directory: {geometry_dir}")

            # Create files directly within the geometry directory
            crystalcontacts_file = geometry_dir / "crystalcontacts_from_colbuilder"
            connect_file = geometry_dir / "connect_from_colbuilder"

            try:
                chimera_scripts_dir = self.file_manager.find_file("chimera_scripts")
            except FileNotFoundError:
                chimera_scripts_dir = self.file_manager.ensure_dir("chimera_scripts")
                LOG.warning(
                    f"Chimera scripts directory not found, created one at: {chimera_scripts_dir}"
                )

            chimera = Chimera(config, str(pdb_path))
            chimera.matrixget(
                pdb=str(pdb_path),
                contact_distance=config.contact_distance,
                crystalcontacts=str(crystalcontacts_file),
            )

            # Check if a .txt extension was added
            if crystalcontacts_file.with_suffix(".txt").exists():
                crystalcontacts_file = crystalcontacts_file.with_suffix(".txt")

            if not crystalcontacts_file.exists():
                raise GeometryGenerationError(
                    message=f"Failed to generate crystalcontacts file at: {crystalcontacts_file}",
                    error_code="GEO_ERR_002",
                )

            crystalcontacts = CrystalContacts(str(crystalcontacts_file))
            system = await self._build_system_structure(crystal, crystalcontacts)

            system, connect = await self._connect_system(system, str(connect_file))

            if self._needs_optimization(system):
                LOG.info(f"{Fore.BLUE}Optimizing system{Style.RESET_ALL}")
                system = await self._optimize_system(
                    system, connect, config.solution_space
                )

                optimized_crystalcontacts_file = crystalcontacts_file.with_name(
                    f"{crystalcontacts_file.stem}_opt"
                )
                crystalcontacts.crystalcontacts_file = str(
                    optimized_crystalcontacts_file
                )
                LOG.debug(
                    f"Updated CrystalContacts file path for optimization: {crystalcontacts.crystalcontacts_file}"
                )

            return system, crystalcontacts, connect

        except Exception as e:
            LOG.error(f"Error in _build_from_contact_distance: {str(e)}")
            raise GeometryGenerationError(
                message="Failed to build system from contact distance",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "contact_distance": config.contact_distance,
                    "pdb_file": str(pdb_path) if "pdb_path" in locals() else None,
                    "working_dir": str(config.working_directory),
                },
            )

    async def _build_from_crystal_contacts(
        self, crystal: Crystal, config: ColbuilderConfig
    ) -> Tuple[System, CrystalContacts, Connect]:
        """
        Build system from existing crystal contacts file.

        Args:
            crystal: Crystal structure
            config: Configuration settings

        Returns:
            Tuple containing system components

        Raises:
            GeometryGenerationError: If building fails
        """
        try:
            crystalcontacts_file = Path(config.crystalcontacts_file)
            if (
                not crystalcontacts_file.exists()
                and crystalcontacts_file.with_suffix(".txt").exists()
            ):
                crystalcontacts_file = crystalcontacts_file.with_suffix(".txt")

            LOG.info(f"Using CrystalContacts file path: {crystalcontacts_file}")

            # Setup connect file
            geometry_dir = self.file_manager.ensure_geometry_dir()
            if config.connect_file:
                connect_file = Path(config.connect_file)
                LOG.info(f"Using provided connect file: {connect_file}")
            else:
                connect_file = geometry_dir / "connect_from_colbuilder"
                LOG.info(f"Using default connect file path: {connect_file}")

            try:
                crystalcontacts = CrystalContacts(str(crystalcontacts_file))
                LOG.info(f"Created CrystalContacts object successfully")
            except Exception as e:
                raise GeometryGenerationError(
                    message="Failed to load crystal contacts file",
                    original_error=e,
                    error_code="GEO_ERR_002",
                    context={"crystalcontacts_file": str(crystalcontacts_file)},
                )

            LOG.info(f"Building system structure")
            system = await self._build_system_structure(crystal, crystalcontacts)
            LOG.info(f"Connecting system")
            system, connect = await self._connect_system(system, str(connect_file))

            if config.crystalcontacts_optimize:
                LOG.info(f"{Fore.BLUE}Optimizing system{Style.RESET_ALL}")
                system = await self._optimize_system(
                    system, connect, config.solution_space
                )
                optimized_path = str(crystalcontacts_file) + "_opt"
                crystalcontacts.crystalcontacts_file = optimized_path
                LOG.info(f"Using optimized crystalcontacts file: {optimized_path}")

            if config.connect_file:
                try:
                    LOG.info(f"Using external connect file: {config.connect_file}")
                    system = connect.get_external_connect_file(
                        system=system, connect_file=str(connect_file)
                    )
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to process external connect file",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"connect_file": str(connect_file)},
                    )

            return system, crystalcontacts, connect

        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Failed to build system from crystal contacts: {str(e)}")
            raise GeometryGenerationError(
                message="Failed to build system from crystal contacts",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "crystalcontacts_file": str(config.crystalcontacts_file),
                    "connect_file": (
                        str(config.connect_file) if config.connect_file else None
                    ),
                },
            )

    async def _build_system_structure(
        self, crystal: Crystal, crystalcontacts: CrystalContacts
    ) -> System:
        """
        Build core system structure from crystal and contacts.

        Args:
            crystal: Crystal structure
            crystalcontacts: Crystal contacts information

        Returns:
            System: Built system

        Raises:
            GeometryGenerationError: If structure building fails
        """
        try:
            system = System(crystal=crystal, crystalcontacts=crystalcontacts)
            transformation = system.crystalcontacts.read_t_matrix()

            unit_cell = {
                model_id: system.crystal.get_s_matrix(t_matrix=transformation[model_id])
                for model_id in transformation
            }

            for model_id in transformation:
                try:
                    model = Model(
                        id=model_id,
                        transformation=transformation[model_id],
                        unit_cell=unit_cell[model_id],
                        pdb_file=crystal.pdb_file,
                    )
                    system.add_model(model=model)
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to create and add model",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"model_id": model_id},
                    )

            return system

        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Failed to build system structure: {str(e)}")
            raise GeometryGenerationError(
                message="Failed to build system structure",
                original_error=e,
                error_code="GEO_ERR_002",
                context={
                    "crystal_pdb": crystal.pdb_file,
                    "crystalcontacts_file": crystalcontacts.crystalcontacts_file,
                },
            )

    async def _connect_system(
        self, system: System, connect_file: str
    ) -> Tuple[System, Connect]:
        """
        Establish connections in the system.

        Args:
            system: System to connect
            connect_file: Connection file path

        Returns:
            Tuple containing connected system and connection information

        Raises:
            GeometryGenerationError: If connection fails
        """
        try:
            connect = Connect(system=system, connect_file=connect_file)
            system_connect = connect.run_connect(system=system)

            connection_count = 0
            for model_id in system_connect:
                try:
                    model_connections = system_connect[model_id]
                    system.get_model(model_id=model_id).add_connect(
                        connect_id=model_id, connect=model_connections
                    )
                    connection_count += (
                        len(model_connections)
                        if isinstance(model_connections, dict)
                        else 1
                    )
                except Exception as e:
                    raise GeometryGenerationError(
                        message="Failed to add connection to model",
                        original_error=e,
                        error_code="GEO_ERR_002",
                        context={"model_id": model_id},
                    )

            LOG.debug(f"Added {connection_count} connections to system")
            return system, connect

        except GeometryGenerationError:
            raise
        except Exception as e:
            LOG.error(f"Failed to connect system: {str(e)}")
            raise GeometryGenerationError(
                message="Failed to connect system",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"connect_file": connect_file},
            )

    def _needs_optimization(self, system: System) -> bool:
        """
        Check if system needs optimization.

        Args:
            system: System to check

        Returns:
            bool: True if system has crosslinks
        """
        needs_opt = any(
            hasattr(system.get_model(model_id), "crosslink")
            and system.get_model(model_id).crosslink
            for model_id in system.get_models()
        )

        if not needs_opt:
            LOG.debug("System has no crosslinks, optimization not needed")

        return needs_opt

    async def _optimize_system(
        self, system: System, connect: Connect, solution_space: List[float]
    ) -> System:
        """
        Optimize system structure.

        Args:
            system: System to optimize
            connect: System connections
            solution_space: Optimization parameters

        Returns:
            System: Optimized system

        Raises:
            GeometryGenerationError: If optimization fails
        """
        try:
            optimizer = Optimizer(system=system, solution_space=solution_space)
            optimized_system = optimizer.run_optimize(system=system, connect=connect)
            optimized_system, _ = await self._connect_system(
                optimized_system, connect.connect_file
            )
            LOG.debug(f"System optimization complete")
            return optimized_system

        except Exception as e:
            LOG.error(f"System optimization failed: {str(e)}")
            raise GeometryGenerationError(
                message="System optimization failed",
                original_error=e,
                error_code="GEO_ERR_002",
                context={"solution_space": solution_space},
            )
