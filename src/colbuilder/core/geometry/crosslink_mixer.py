from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Set
import shutil
from colorama import Fore, Style
import os
import traceback

from colbuilder.core.geometry.crystal import Crystal
from colbuilder.core.geometry.system import System
from colbuilder.core.geometry.chimera import Chimera
from colbuilder.core.geometry.connect import Connect
from colbuilder.core.geometry.caps import Caps
from colbuilder.core.geometry.crystalcontacts import CrystalContacts
from colbuilder.core.geometry.mix import Mix
from colbuilder.core.geometry.optimize import Optimizer
from colbuilder.core.geometry.unpaired_crosslinks import UnpairedCrosslinkFinder
from colbuilder.core.utils.config import ColbuilderConfig
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class CrosslinkMixer:
    """
    A class for mixing and managing different crosslink types in collagen microfibrils.

    The CrosslinkMixer handles the generation, modification, and mixing of collagen systems
    with different crosslink types (e.g., divalent and trivalent). It manages:
    - Crystal structure transformations
    - System building from PDB files
    - Crosslink mixing according to specified ratios
    - Terminal capping of collagen models
    - Length control of the resulting fibril

    Attributes
    ----------
    path_wd : Optional[Path]
        Working directory for temporary and output files
    fibril_length : Optional[float]
        Length of the collagen microfibril in nanometers
    contact_distance : Optional[float]
        Distance threshold for identifying crystal contacts in Angstroms
    crystalcontacts_file : str
        Base filename for crystal contacts (default: 'crystalcontacts_from_colbuilder')

    Methods
    -------
    mix(system, config, temp_dir)
        Main method to mix different crosslink types according to ratios
    _build_system(crystal, crystalcontacts)
        Build a system from crystal structure and contacts
    _cap_system(system, crosslink_type, temp_dir)
        Add terminal caps to collagen models
    _matrixset_system(system)
        Update system coordinates and cut to specified length

    Examples
    --------
    >>> mixer = CrosslinkMixer()
    >>> config = ColbuilderConfig(
    ...     ratio_mix={"D": 80, "T": 20},
    ...     files_mix=["human-D.pdb", "human-T.pdb"],
    ...     fibril_length=40.0
    ... )
    >>> system = await mixer.mix(None, config)

    Notes
    -----
    - Requires proper PDB files with collagen triple helix structures
    - Crystal contacts are handled for space group 1 only
    - Fibril length is specified in nanometers
    - Contact distance is specified in Angstroms
    """

    def __init__(self):
        """Initialize CrosslinkMixer with default values."""
        self.path_wd: Optional[Path] = None
        self.fibril_length: Optional[float] = None
        self.contact_distance: Optional[float] = None
        self.crystalcontacts_file: str = "crystalcontacts_from_colbuilder"
        self.config: Optional[ColbuilderConfig] = None

    @staticmethod
    def _ensure_pdb_extension(filename: str) -> str:
        """Ensure filename has .pdb extension."""
        if not filename.endswith(".pdb"):
            return filename + ".pdb"
        return filename

    @staticmethod
    def _build_system(
        crystal: Crystal, crystalcontacts: Optional[CrystalContacts] = None
    ) -> System:
        """Build a system from crystal and crystal contacts."""
        system = System(crystal=crystal, crystalcontacts=crystalcontacts)

        if crystalcontacts is None:
            LOG.warning("No crystal contacts provided. Adding a default model.")
            from colbuilder.core.geometry.model import Model

            default_transformation = crystal.get_default_transformation()
            default_unit_cell = crystal.get_s_matrix(t_matrix=default_transformation)
            default_model = Model(
                id=0,
                transformation=default_transformation,
                unit_cell=default_unit_cell,
                pdb_file=crystal.pdb_file,
            )
            system.add_model(model=default_model)
        else:
            transformation = system.crystalcontacts.read_t_matrix()
            unit_cell: Dict[float, Any] = {
                k: system.crystal.get_s_matrix(t_matrix=transformation[k])
                for k in transformation
            }

            from colbuilder.core.geometry.model import Model

            for key_m in transformation:
                model = Model(
                    id=key_m,
                    transformation=transformation[key_m],
                    unit_cell=unit_cell[key_m],
                    pdb_file=crystal.pdb_file,
                )
                system.add_model(model=model)

        LOG.debug(f"Built system with {len(system.get_models())} models")
        return system

    def _build_from_contactdistance(
        self,
        path_wd: Path,
        pdb_file: str,
        contact_distance: float,
        solution_space: List[float],
        crystalcontacts_file: Optional[str],
        chimera: Chimera,
        crystal: Crystal,
    ) -> Tuple[System, CrystalContacts, Connect]:
        """Build system from contact distance."""
        path_pdb_file = path_wd / pdb_file
        connect_file = "connect_from_colbuilder"

        LOG.debug(
            f"     Getting CrystalContacts for contact distance {contact_distance} Ang"
        )
        chimera.matrixget(
            pdb=str(path_pdb_file),
            contact_distance=contact_distance,
            crystalcontacts=self.crystalcontacts_file,
        )

        LOG.debug(f"     Writing {self.crystalcontacts_file}")
        crystalcontacts = CrystalContacts(self.crystalcontacts_file)

        LOG.info("     Building system")
        system = self._build_system(crystal=crystal, crystalcontacts=crystalcontacts)

        LOG.info("     Connecting system")
        system, connect = self._connect_system(system=system, connect_file=connect_file)

        has_crosslinks = any(
            hasattr(system.get_model(model_id=model_id), "crosslink")
            and system.get_model(model_id=model_id).crosslink
            for model_id in system.get_models()
        )

        if has_crosslinks:
            LOG.info("     Optimizing system")
            optimizer = Optimizer(system=system, solution_space=solution_space)
            system = optimizer.run_optimize(system=system, connect=connect)
            system, connect = self._connect_system(
                system=system, connect_file=connect_file
            )
            crystalcontacts.crystalcontacts_file = self.crystalcontacts_file + "_opt"
        else:
            LOG.info("     Skipping optimization for non-crosslinked system")

        return system, crystalcontacts, connect

    @staticmethod
    def _connect_system(
        system: System, connect_file: Optional[str] = None
    ) -> Tuple[System, Connect]:
        """Connect models in the system."""
        LOG.debug("Identifying crosslink connections")
        connect = Connect(system=system, connect_file=connect_file)
        system_connect = connect.run_connect(system=system)

        for key_m in system_connect:
            system.get_model(model_id=key_m).add_connect(
                connect_id=key_m, connect=system_connect[key_m]
            )

        return system, connect

    def _auto_fix_unpaired(
        self,
        geom_dir: Path,
        base_dir: Path,
        allowed_resnames: Optional[Set[str]] = None,
        config: Optional[ColbuilderConfig] = None,
    ) -> None:
        """
        Detect and mutate unpaired markers in ``geom_dir`` caps files.
        """
        try:
            fixer = UnpairedCrosslinkFinder(
                base_dir=base_dir,
                geom_dir=geom_dir,
                allowed_resnames=allowed_resnames,
            )
            entries, _ = fixer.run()
            if entries:
                LOG.info("Auto-fixing %d unpaired markers under %s", len(entries), geom_dir)
                if not self._apply_chimera_swaps(entries, geom_dir, config=config):
                    raise RuntimeError(f"Chimera swapaa failed for unpaired markers in {geom_dir}")
                LOG.debug("Chimera swapaa applied for unpaired markers in %s", geom_dir)
            else:
                LOG.info("No unpaired markers found in %s", geom_dir)
        except Exception as e:
            LOG.error("Auto-fix for unpaired markers in %s failed: %s", geom_dir, e)
            raise

    def _apply_chimera_swaps(
        self, entries: List[str], system_dir: Path, config: Optional[ColbuilderConfig] = None
    ) -> bool:
        """
        Use Chimera swapaa to rebuild sidechains for unpaired markers.

        Returns True on success, False otherwise.
        """
        try:
            cfg = config or self.config
            if cfg is None:
                return False

            replace_file = system_dir / "replace_unpaired.txt"
            lines: List[str] = []
            for line in entries:
                parts = line.split()
                if len(parts) < 4:
                    continue
                pdb_name, _, resid, chain = parts[:4]
                orig_res = parts[1].upper()
                target_res = "ARG" if orig_res in {"AGS", "APD"} else "LYS"
                lines.append(f"{pdb_name} {target_res} {resid} {chain}")

            if not lines:
                return False

            replace_file.write_text("\n".join(lines) + "\n")

            chim = Chimera(cfg, pdb=str(system_dir))
            result = chim.swapaa(replace=str(replace_file), system_type=str(system_dir))
            if result.returncode != 0:
                LOG.error("Chimera swapaa failed for %s: %s", system_dir, result.stderr.decode() if hasattr(result, "stderr") else result.stderr)
                return False
            return True
        except Exception as e:
            LOG.warning("Chimera swapaa application failed: %s", e)
            return False    

    def _cap_system(
        self, system: System, crosslink_type: str, temp_dir: Optional[Path] = None
    ) -> Caps:
        """Add caps to the system."""

        if temp_dir is None:
            LOG.warning(
                "No temporary directory provided for capping. Using current directory."
            )
            temp_dir = Path.cwd()

        caps_dir = temp_dir / str(crosslink_type)
        if not caps_dir.exists():
            caps_dir.mkdir(parents=True, exist_ok=True)
            LOG.debug(f"Created caps directory: {caps_dir}")

        original_dir = os.getcwd()
        try:
            os.chdir(caps_dir)

            caps = Caps(system=system)

            model_ids = list(system.get_models())

            LOG.debug(f"    Writing capped PDB files to {caps_dir}")
            for idx in model_ids:
                try:
                    pdb_id = int(float(idx))
                    pdb_file = f"{pdb_id}.pdb"

                    if not os.path.exists(pdb_file):
                        LOG.warning(f"PDB file {pdb_file} not found in {caps_dir}.")
                        continue
                    caps.read_residues(pdb_id=pdb_id)
                    caps.add_caps(pdb_id=pdb_id, crosslink_type=crosslink_type)
                except Exception as e:
                    LOG.error(f"Error processing caps for model {idx}: {str(e)}")

            return caps

        except Exception as e:
            LOG.error(f"Error in capping operation: {str(e)}")
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise
        finally:
            os.chdir(original_dir)

    def _matrixset_system(self, system: System) -> System:
        """Update system after cutting to specified length."""

        id_file = Path(f"{self.crystalcontacts_file}_opt_id.txt")

        if not id_file.exists():
            raise FileNotFoundError(f"Crystal contacts ID file not found: {id_file}")

        try:
            with open(id_file, "r") as f:
                contacts = [float(i.split(" ")[1]) for i in f.readlines()]
        except Exception as e:
            LOG.error(f"Error reading crystal contacts ID file: {str(e)}")
            raise

        for model in list(system.get_models()):
            if model not in contacts:
                system.delete_model(model_id=model)
            elif system.get_model(model_id=model).connect is not None:
                for connect in list(system.get_model(model_id=model).connect):
                    if connect not in contacts:
                        system.get_model(model_id=model).delete_connect(
                            connect_id=connect
                        )

        return system

    async def mix(
        self,
        system: Optional[System],
        config: ColbuilderConfig,
        temp_dir: Optional[Path] = None,
    ) -> Tuple[Optional[System], Optional[Path]]:
        """Main method to mix different crosslink types."""
        original_dir = os.getcwd()

        try:
            if not system:
                system = System()

            self.path_wd = temp_dir if temp_dir else Path(config.working_directory)
            if not self.path_wd.exists():
                self.path_wd.mkdir(parents=True, exist_ok=True)

            os.chdir(self.path_wd)

            self.fibril_length = config.fibril_length
            self.contact_distance = config.contact_distance

            if isinstance(config.ratio_mix, str):
                ratio_dict = {}
                for part in config.ratio_mix.split():
                    if ":" in part:
                        key, value = part.split(":")
                        try:
                            ratio_dict[key] = int(value)
                        except ValueError:
                            LOG.error(f"Invalid ratio value in {part}")
                            ratio_dict[key] = 0
                config.ratio_mix = ratio_dict

            if not config.ratio_mix or not isinstance(config.ratio_mix, dict):
                LOG.error(f"Invalid ratio_mix format in config: {config.ratio_mix}")
                raise ValueError(
                    f"ratio_mix must be a dictionary, got {type(config.ratio_mix)}"
                )

            LOG.debug(
                f"Mix setup: {Fore.MAGENTA}{config.ratio_mix} ('TypeA': %, 'TypeB': %){Style.RESET_ALL}"
            )

            if not config.files_mix or len(config.files_mix) < len(config.ratio_mix):
                LOG.error(f"Not enough files in files_mix: {config.files_mix}")
                raise ValueError(
                    f"files_mix must contain at least {len(config.ratio_mix)} files"
                )

            mix_pdb = {}
            for i, key in enumerate(config.ratio_mix.keys()):
                if i < len(config.files_mix):
                    mix_pdb[key] = config.files_mix[i]

            LOG.info("Step 1/2 Generating mix setup")

            if system.get_size() == 0:
                first_pdb = str(config.files_mix[0])
                if first_pdb.endswith(".pdb"):
                    first_pdb = first_pdb[:-4]

                LOG.debug(
                    f"No system provided, building initial system from {first_pdb}"
                )
                crystal = Crystal(first_pdb)
                crystal.translate_crystal(pdb=first_pdb, translate=[0, 0, 4000])

                chimera = Chimera(config, str(self.path_wd / first_pdb))

                system, crystalcontacts, connect = self._build_from_contactdistance(
                    self.path_wd,
                    first_pdb,
                    self.contact_distance,
                    config.solution_space,
                    None,
                    chimera,
                    crystal,
                )

                LOG.debug(f"Writing {crystalcontacts.crystalcontacts_file}")
                crystalcontacts.write_crystalcontacts(
                    system=system,
                    crystalcontacts_file=crystalcontacts.crystalcontacts_file,
                )

            else:
                LOG.debug("Using provided system")

            system_size = system.get_size()

            AGE_RESN = {"AGS", "LGX", "APD", "LPS"}

            if hasattr(system, "crystalcontacts") and system.crystalcontacts:
                LOG.debug(
                    f"Crystalcontacts file: {system.crystalcontacts.crystalcontacts_file}"
                )

            for key in list(config.ratio_mix.keys()):
                if key not in mix_pdb:
                    LOG.warning(f"No PDB file mapped for type {key}, skipping")
                    continue

                pdb_file = mix_pdb[key]

                pdb_path = Path(pdb_file)
                if not pdb_path.exists():
                    pdb_path = self.path_wd / pdb_file
                    if not pdb_path.exists():
                        LOG.error(f"PDB file not found: {pdb_file}")
                        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

                crystal = Crystal(pdb=str(pdb_path))
                crystal.translate_crystal(pdb=str(pdb_path), translate=[0, 0, 4000])

                chimera = Chimera(config, str(pdb_path))

                LOG.info(f" - System {key}:")
                LOG.info(f"    Generating system from {pdb_file}")

                crystalcontacts_file = (
                    system.crystalcontacts.crystalcontacts_file
                    if hasattr(system, "crystalcontacts") and system.crystalcontacts
                    else self.crystalcontacts_file
                )

                if (
                    not Path(crystalcontacts_file).exists()
                    and not Path(f"{crystalcontacts_file}.txt").exists()
                ):
                    LOG.error(f"Crystalcontacts file not found: {crystalcontacts_file}")
                    raise FileNotFoundError(
                        f"Crystalcontacts file not found: {crystalcontacts_file}"
                    )

                fibril_length_nm = float(self.fibril_length)

                chimera.matrixset(
                    pdb=str(pdb_path),
                    crystalcontacts=str(crystalcontacts_file),
                    system_size=system_size,
                    fibril_length=fibril_length_nm,
                )

                LOG.info(f"     Cutting system to {fibril_length_nm} nm")
                system = self._matrixset_system(system)

                type_dir = temp_dir / str(key)
                type_dir.mkdir(parents=True, exist_ok=True)

                # Copy opt_id file into this type dir for downstream unpaired detection
                opt_id_src = self.path_wd / "crystalcontacts_from_colbuilder_opt_id.txt"
                if opt_id_src.exists():
                    try:
                        shutil.copy2(opt_id_src, type_dir / opt_id_src.name)
                    except Exception:
                        pass

                for model_id in system.get_models():
                    model_id_int = int(float(model_id))
                    source_pdb = Path(f"{model_id_int}.pdb")
                    target_pdb = type_dir / f"{model_id_int}.pdb"

                    if source_pdb.exists():
                        try:
                            shutil.copy2(source_pdb, target_pdb)
                        except Exception as e:
                            LOG.error(
                                f"Failed to copy model PDB file for model {model_id}: {str(e)}"
                            )
                    else:
                        LOG.warning(
                            f"Transformed PDB file not found for model {model_id}: {source_pdb}"
                        )

                LOG.info("     Adding caps")
                self._cap_system(system, key, temp_dir)

                # Write connect file for this type (per-type auto-fix handled after mixing)
                Connect(system=system).write_connect(
                    system=system, connect_file=type_dir / "connect_from_colbuilder"
                )

                # Auto-fix unpaired crosslinks per-type before mixing
                self._auto_fix_unpaired(
                    geom_dir=type_dir,
                    base_dir=type_dir,
                    allowed_resnames=None,
                    config=config,
                )

            # After fixing per-type caps, compute average AGE markers per model for each type
            per_type_age_avg: Dict[str, float] = {}
            for key in mix_pdb.keys():
                type_dir = temp_dir / str(key)
                caps_files = list(type_dir.glob("*.caps.pdb"))
                if not caps_files:
                    per_type_age_avg[key] = 1.0
                    continue
                marker_count = 0
                for cf in caps_files:
                    try:
                        with cf.open("r") as fh:
                            for line in fh:
                                if not line.startswith(("ATOM", "HETATM")):
                                    continue
                                if line[17:20].strip() in AGE_RESN:
                                    marker_count += 1
                    except Exception:
                        continue
                per_type_age_avg[key] = marker_count / float(len(caps_files) or 1) or 1.0

            LOG.info("Step 2/2 Mixing systems")

            mix_ = Mix(ratio_mix=config.ratio_mix, system=system)
            system = mix_.add_mix(system=system)

            type_counts = {}
            for model_id in system.get_models():
                model = system.get_model(model_id=model_id)
                model_type = getattr(model, "type", "UNKNOWN")
                type_counts[model_type] = type_counts.get(model_type, 0) + 1

                if not hasattr(model, "transformation") or model.transformation is None:
                    LOG.warning(
                        f"Model {model_id} has no transformation matrix after mixing"
                    )

            connect_file_path = self.path_wd / "connect_from_colbuilder.txt"
            LOG.debug(f"Writing connect file {connect_file_path}")
            Connect(system=system).write_connect(
                system=system, connect_file=connect_file_path
            )

            # Auto-fix unpaired AGE markers in the mixed caps before final output
            try:
                fixer = UnpairedCrosslinkFinder(
                    base_dir=self.path_wd,
                    geom_dir=temp_dir,
                    allowed_resnames=None,
                )
                entries, _ = fixer.run()
                if entries:
                    LOG.info("Removing %d unpaired markers in mixed caps", len(entries))
                    if not self._apply_chimera_swaps(entries, temp_dir, config=config):
                        raise RuntimeError("Chimera swapaa failed for unpaired markers in mixed caps")
            except Exception as e:
                LOG.error("Auto-fix for unpaired markers failed: %s", e)
                raise
            
            temp_pdb = temp_dir / f"{config.output or 'collagen_fibril'}.pdb"

            LOG.debug(f"Writing temporary output PDB to {temp_pdb}")
            system.write_pdb(
                pdb_out=temp_pdb,
                fibril_length=self.fibril_length,
                temp_dir=temp_dir,
            )

            return system, temp_pdb

        except Exception as e:
            LOG.error(f"Exception in mix: {str(e)}")
            LOG.error(f"Exception type: {type(e)}")
            LOG.debug(f"Traceback: {traceback.format_exc()}")
            raise

        finally:
            os.chdir(original_dir)
