# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Union
import os

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class IncompleteCrosslinkError(ValueError):
    """Raised when marker residues cannot form a complete crosslink topology."""


class Crosslink:
    """
    Setup crosslink topology for collagen models.

    This class handles the identification and parameterization of crosslinks
    in collagen molecular models. It processes merged PDB files to identify
    crosslink sites and generates the necessary bonded parameters (bonds,
    angles, dihedrals) for crosslinks in Martini coarse-grained models.

    The class supports four types of crosslinks:
    - Divalent HLKNL-crosslinks between L4Y and L5Y residues
    - Trivalent PYD-crosslinks between LYX, LY2, and LY3 residues
    - Glucosepane crosslinks between the LGX and AGS marker residues
    - MOLD crosslinks between the LZS and LZD marker residues
    """

    # Glucosepane terms that cross the LGX/AGS marker boundary, per Martini 3 parameter set
    # (glucosepane_martini_timestep).  All other glucosepane interactions are defined inside the
    # marker templates (data/topology/martini300C-glucosepane/glucosepane_<N>fs.ff).
    # s1_r1_r4 = None: the set has no S1--R1--R4 angle.
    GLUCOSEPANE_CROSS_TERMS = {
        # 20 fs: model A_b2 with two virtual sites (R2 = AGS:SC6, R5 = AGS:SC3), fitted on the
        # all-atom collagen fibril (GCP_A_2VS_fib_b2.itp, 17/10/2026); no S1--R1--R4 angle.
        20: dict(s1_r1=('0.2609', '7335.8'), b1_s1_r1=('129.113', '201.92'),
                 s1_r1_r2=('122.006', '44.43'), s1_r1_r4=None),
        # G21 20 fs model, commented out 17/10/2026 (replaced by A_b2)
        # 20: dict(s1_r1=('0.295', '63395.824'), b1_s1_r1=('141.879', '127.334'),
        #          s1_r1_r2=('89.976', '48.526'), s1_r1_r4=('132.752', '89.996')),
        # 10 fs: model without virtual sites (all ten beads carry mass)
        10: dict(s1_r1=('0.297', '57650.518'), b1_s1_r1=('138.898', '112.282'),
                 s1_r1_r2=('106.601', '52.839'), s1_r1_r4=('152.855', '150')),
    }

    def __init__(self, cnt_model: Optional[int] = None, glucosepane_timestep: int = 20) -> None:
        """
        Initialize the Crosslink object with model parameters.

        Parameters
        ----------
        cnt_model : Optional[int]
            Model counter used for file naming. If provided, the merged PDB
            file will be named '{cnt_model}.merge.pdb'
        glucosepane_timestep : int
            Glucosepane Martini 3 parameter set: 20 (two virtual sites, default) or
            10 (no virtual sites); must match the installed glucosepane.ff.
        """
        self.file: str = f"{int(cnt_model)}.merge.pdb" if cnt_model is not None else ""
        self.crosslink_coords: List[List[float]] = []
        self.crosslink_pdb: List[List[Any]] = []  # Mix of str and float
        self.residue_beads: Dict[int, Dict[str, str]] = {}
        self.crosslink_neighbors: List[Any] = []
        self.crosslink_connect: List[List[List[Any]]] = []
        self.crosslink_pairs: List[Tuple[List[Any], List[Any]]] = []  # Store valid pairs
        self.crosslink_bonded: Dict[str, List[List[Any]]] = {
            'bonds': [],
            'angles': [],
            'dihedrals': []
        }

        # Distance thresholds for different crosslink types (in Å).
        # Equilibrium bead-bead distances are ~2.3-4.2 A (see dlyxly3/dlyxly2/dly45
        # below). A loose 15 A cutoff lets one bead pair with several partners in
        # crosslink-dense regions, creating spurious bonds. Keep the cutoff modest
        # so only genuine partners (close in the assembled fibril) are bonded.
        self.crosslink_thresholds = {
            'LYX_LY2': 6.0,  # Maximum distance for LYX SC4/SC5 - LY2 SC1 crosslinks
            'LYX_LY3': 6.0,  # Maximum distance for LYX SC5 - LY3 SC1 crosslinks
            'L4Y_L5Y': 6.0,  # Maximum distance for L4Y SC1 - L5Y SC2 crosslinks
            'LGX_AGS': 6.0,  # Maximum distance for glucosepane S1 - R1
            'LZS_LZD': 6.0   # Maximum distance for MOLD R2 - S2
        }

        # HLKNL-crosslink parameters
        self.dly45: str = '0.415'    # L4Y-L5Y bond equilibrium distance (nm)
        self.kly45: str = '7000'     # L4Y-L5Y bond force constant (kJ/mol/nm^2)
        self.al45y_1: str = '140'    # L4Y-L5Y SC1-SC2-SC1(L4Y) angle (degrees)
        self.al45y_2: str = '140'    # L4Y-L5Y SC2(L5Y)-SC1-BB angle (degrees)
        self.k_angle: str = '153'    # Universal angle force constant (kJ/mol/rad^2)
        
        # PYD-crosslink parameters
        self.klyxly2: str = '11000'   # LYX-LY2 bond force constant (kJ/mol/nm^2)
        self.klyxly3: str = '12000'  # LYX-LY3 bond force constant (kJ/mol/nm^2)
        self.klyx5ly2: str = '12000'  # LYX-LY2 bond force constant (kJ/mol/nm^2) (conserve ring)
        self.dlyxly2: str = '0.290'  # LYX-LY2 bond equilibrium distance (nm)
        self.dlyxly3: str = '0.230'  # LYX-LY3 bond equilibrium distance (nm)
        self.dlyx5ly2: str = '0.370'  # LYX-LY2 bond equilibrium distance (nm) (added)
        
        # PYD-crosslink angle parameters (degrees)
        self.al2yx_1: str = '60'    # LY2-LYX TP1q-TC6q-TC4 angle R2 R4 R3
        self.al2yx_2: str = '100'     # LY2-LYX TQ2p-TP1q-TC4 angle R1 R2 R3
        self.al2yx_3: str = '100'    # LY2-LYX TP1q-TC4-SP2 angle R2 B2 R3
        self.al3yx_1: str = '100'    # LY3-LYX TP1q-TC6q-TC4 angle R2 R4 R3/S3
        self.al3yx_2: str = '140'    # LY3-LYX TQ2p-TC6q-TC4 angle R1 R4 S3
        self.al3yx_3: str = '130'    # LY3-LYX TC6q-TC4-SP2 angle R4 S3 B3

        # Glucosepane terms that cross the LGX/AGS marker boundary (selected parameter set).
        if int(glucosepane_timestep) not in self.GLUCOSEPANE_CROSS_TERMS:
            raise ValueError(f"glucosepane_timestep must be 10 or 20, got {glucosepane_timestep}")
        self.glucosepane_timestep: int = int(glucosepane_timestep)
        gcp = self.GLUCOSEPANE_CROSS_TERMS[self.glucosepane_timestep]
        self.dgcp_s1_r1: str = gcp['s1_r1'][0]
        self.kgcp_s1_r1: str = gcp['s1_r1'][1]
        self.agcp_b1_s1_r1: Tuple[str, str] = gcp['b1_s1_r1']
        self.agcp_s1_r1_r2: Tuple[str, str] = gcp['s1_r1_r2']
        self.agcp_s1_r1_r4: Optional[Tuple[str, str]] = gcp['s1_r1_r4']

        # MOLD terms crossing the LZS/LZD marker boundary.  The renamed
        # standalone aliases are R2=LZS:SC3 and S2=LZD:SC1.
        self.dmold_r2_s2: str = '0.310'
        self.kmold_r2_s2: str = '18770.020'
        self.amold_r1_r2_s2: Tuple[str, str] = ('117.168', '61.307')
        self.amold_r3_r2_s2: Tuple[str, str] = ('90.680', '48.532')
        self.amold_r2_s2_b2: Tuple[str, str] = ('137.413', '105.966')

        LOG.debug(f"Initialized Crosslink for model counter {cnt_model}")

    def get_crosslink_coords(self, cnt_model: Optional[int] = None) -> List[List[float]]:
        """
        Extract coordinates of crosslink sites from a PDB file.

        Following the working version pattern: stores coordinates as floats.
        """
        if cnt_model is None:
            file = self.file
        else:
            file = f"{int(cnt_model)}.merge.pdb"

        LOG.debug(f"Reading crosslink coordinates from {file}")

        if not os.path.exists(file):
            LOG.error(f"PDB file not found: {file}")
            return self.crosslink_coords

        self.crosslink_coords = []
        self.crosslink_pdb = []
        self.residue_beads = {}
        it_pdb = 0
        residue_occurrence = -1
        previous_residue = None

        try:
            with open(file, 'r') as f:
                for line in f:
                    if not line.startswith('ATOM'):
                        continue

                    it_pdb += 1
                    resname = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    residue_key = (resname, line[21:26])
                    if residue_key != previous_residue:
                        residue_occurrence += 1
                        previous_residue = residue_key
                    self.residue_beads.setdefault(residue_occurrence, {})[
                        atom_name
                    ] = str(it_pdb)

                    if ((resname == 'LYX' and atom_name in ['SC4', 'SC5']) or
                        (resname in ['LY2', 'LY3'] and atom_name == 'SC1') or
                        (resname == 'L4Y' and atom_name == 'SC1') or
                        (resname == 'L5Y' and atom_name == 'SC2') or
                        (resname == 'LGX' and atom_name == 'SC1') or
                        (resname == 'AGS' and atom_name == 'SC5') or
                        (resname == 'LZS' and atom_name == 'SC3') or
                        (resname == 'LZD' and atom_name == 'SC1')):

                        coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

                        self.crosslink_pdb.append([
                            str(it_pdb),           # atom index as string
                            resname,               # residue name
                            atom_name,              # atom name
                            line[21:26],           # chain info
                            residue_occurrence,    # unambiguous residue occurrence
                            coords[0],             # x coordinate as float
                            coords[1],             # y coordinate as float
                            coords[2]              # z coordinate as float
                        ])
                        self.crosslink_coords.append(coords)

                        LOG.debug(f"Found crosslink site: {line[17:20].strip()} {line[12:15].strip()} (atom {it_pdb})")

            LOG.debug(f"Found {len(self.crosslink_coords)} crosslink sites")

        except Exception as e:
            LOG.error(f"Error reading PDB file {file}: {str(e)}")

        return self.crosslink_coords

    def _same_residue_bead(self, atom: List[Any], bead_name: str) -> Optional[str]:
        """Return a bead index from the same marker occurrence as ``atom``."""
        residue_occurrence = atom[4]
        return self.residue_beads.get(residue_occurrence, {}).get(bead_name)

    def get_crosslink_connect(self, cnt_model: Optional[int] = None) -> List[List[List[Any]]]:
        """
        Get nearest crosslinks to determine connections.

        This method is kept for backwards compatibility but now uses the improved
        pair-finding algorithm internally.
        """
        LOG.debug("Using improved crosslink pair detection algorithm")

        pairs = self.find_crosslink_pairs(cnt_model=cnt_model)

        self.crosslink_connect = []
        if pairs:
            all_atoms = []
            for pair in pairs:
                if pair[0] not in all_atoms:
                    all_atoms.append(pair[0])
                if pair[1] not in all_atoms:
                    all_atoms.append(pair[1])

            if all_atoms:
                self.crosslink_connect.append(all_atoms)
                LOG.debug(f"Converted {len(pairs)} pairs to connection group with {len(all_atoms)} atoms")

        return self.crosslink_connect

    def _match_nearest_pairs(
        self,
        a_atoms: List[List[Any]],
        b_atoms: List[List[Any]],
        threshold: float,
    ) -> List[Tuple[List[Any], List[Any], float]]:
        """
        Greedy one-to-one nearest-neighbour matching between two atom groups.

        Returns (a_atom, b_atom, distance) for the closest compatible pairs within
        ``threshold``, using each atom at most once (shortest pairs assigned first).
        Matching is per bead-category, so a bead that legitimately participates in
        two different categories (e.g. LYX SC5 bonds both LY2 and LY3) still forms
        both bonds — but within a single category a bead bonds only its one nearest
        partner, preventing spurious extra bonds in crosslink-dense regions.
        """
        candidates: List[Tuple[float, int, int]] = []
        for ai, a in enumerate(a_atoms):
            a_xyz = np.array(a[-3:], dtype=float)
            for bi, b in enumerate(b_atoms):
                dist = float(np.linalg.norm(a_xyz - np.array(b[-3:], dtype=float)))
                if dist <= threshold:
                    candidates.append((dist, ai, bi))

        candidates.sort(key=lambda t: t[0])
        used_a: set = set()
        used_b: set = set()
        matched: List[Tuple[List[Any], List[Any], float]] = []
        for dist, ai, bi in candidates:
            if ai in used_a or bi in used_b:
                LOG.debug(
                    f"    Skipping already-bonded pair {a_atoms[ai][0]}-{b_atoms[bi][0]} "
                    f"(d={dist:.3f} Å); each crosslink bead bonds a single partner per category."
                )
                continue
            used_a.add(ai)
            used_b.add(bi)
            matched.append((a_atoms[ai], b_atoms[bi], dist))
        return matched

    def find_crosslink_pairs(self, cnt_model: Optional[int] = None) -> List[Tuple[List[Any], List[Any]]]:
        """
        Find valid crosslink pairs based on distance and compatibility.

        Pairing is a one-to-one nearest-neighbour match per bead-category
        (LYX-SC4/LY2, LYX-SC5/LY2, LYX-SC5/LY3, L4Y/L5Y, LGX/AGS,
        LZS/LZD): the closest compatible
        pairs are bonded first and each atom is used at most once within a category.
        This avoids the spurious extra bonds that "bond every pair within the cutoff"
        produced, while still allowing the shared LYX SC5 bead to bond both LY2
        (conserve-ring) and LY3.
        """
        self.get_crosslink_coords(cnt_model=cnt_model)

        if not self.crosslink_coords:
            LOG.warning("No crosslink coordinates found")
            return []

        self.crosslink_pairs = []

        lyx_sc4_atoms = []
        lyx_sc5_atoms = []
        ly2_sc1_atoms = []
        ly3_sc1_atoms = []
        l4y_sc1_atoms = []
        l5y_sc2_atoms = []
        lgx_sc1_atoms = []
        ags_sc5_atoms = []
        lzs_sc3_atoms = []
        lzd_sc1_atoms = []

        for atom in self.crosslink_pdb:
            if atom[1] == 'LYX' and atom[2] == 'SC4':
                lyx_sc4_atoms.append(atom)
            elif atom[1] == 'LYX' and atom[2] == 'SC5':
                lyx_sc5_atoms.append(atom)
            elif atom[1] == 'LY2' and atom[2] == 'SC1':
                ly2_sc1_atoms.append(atom)
            elif atom[1] == 'LY3' and atom[2] == 'SC1':
                ly3_sc1_atoms.append(atom)
            elif atom[1] == 'L4Y' and atom[2] == 'SC1':
                l4y_sc1_atoms.append(atom)
            elif atom[1] == 'L5Y' and atom[2] == 'SC2':
                l5y_sc2_atoms.append(atom)
            elif atom[1] == 'LGX' and atom[2] == 'SC1':
                lgx_sc1_atoms.append(atom)
            elif atom[1] == 'AGS' and atom[2] == 'SC5':
                ags_sc5_atoms.append(atom)
            elif atom[1] == 'LZS' and atom[2] == 'SC3':
                lzs_sc3_atoms.append(atom)
            elif atom[1] == 'LZD' and atom[2] == 'SC1':
                lzd_sc1_atoms.append(atom)

        LOG.debug(f"Crosslink atoms:")
        LOG.debug(f"  LYX SC4: {len(lyx_sc4_atoms)}")
        LOG.debug(f"  LYX SC5: {len(lyx_sc5_atoms)}")
        LOG.debug(f"  LY2 SC1: {len(ly2_sc1_atoms)}")
        LOG.debug(f"  LY3 SC1: {len(ly3_sc1_atoms)}")
        LOG.debug(f"  L4Y SC1: {len(l4y_sc1_atoms)}")
        LOG.debug(f"  L5Y SC2: {len(l5y_sc2_atoms)}")
        LOG.debug(f"  LGX SC1: {len(lgx_sc1_atoms)}")
        LOG.debug(f"  AGS SC5: {len(ags_sc5_atoms)}")
        LOG.debug(f"  LZS SC3: {len(lzs_sc3_atoms)}")
        LOG.debug(f"  LZD SC1: {len(lzd_sc1_atoms)}")

        # LYX SC4 - LY2 SC1 (one-to-one nearest neighbour)
        for lyx_atom, ly2_atom, dist in self._match_nearest_pairs(
            lyx_sc4_atoms, ly2_sc1_atoms, self.crosslink_thresholds['LYX_LY2']
        ):
            self.crosslink_pairs.append((lyx_atom, ly2_atom))
            LOG.info(f" Added LYX-LY2 (SC4) pair: atoms {lyx_atom[0]} - {ly2_atom[0]} (distance: {dist:.3f} Å)")

        # LYX SC5 - LY2 SC1 (conserve-ring bond; one-to-one nearest neighbour)
        for lyx_atom, ly2_atom, dist in self._match_nearest_pairs(
            lyx_sc5_atoms, ly2_sc1_atoms, self.crosslink_thresholds['LYX_LY2']
        ):
            self.crosslink_pairs.append((lyx_atom, ly2_atom))
            LOG.info(f" Added LYX-LY2 (SC5) pair: atoms {lyx_atom[0]} - {ly2_atom[0]} (distance: {dist:.3f} Å)")

        # LYX SC5 - LY3 SC1 (one-to-one nearest neighbour)
        for lyx_atom, ly3_atom, dist in self._match_nearest_pairs(
            lyx_sc5_atoms, ly3_sc1_atoms, self.crosslink_thresholds['LYX_LY3']
        ):
            self.crosslink_pairs.append((lyx_atom, ly3_atom))
            LOG.info(f" Added LYX-LY3 (SC5) pair: atoms {lyx_atom[0]} - {ly3_atom[0]} (distance: {dist:.3f} Å)")

        # L4Y SC1 - L5Y SC2 (one-to-one nearest neighbour)
        for l4y_atom, l5y_atom, dist in self._match_nearest_pairs(
            l4y_sc1_atoms, l5y_sc2_atoms, self.crosslink_thresholds['L4Y_L5Y']
        ):
            self.crosslink_pairs.append((l4y_atom, l5y_atom))
            LOG.info(f" Added L4Y-L5Y pair: atoms {l4y_atom[0]} - {l5y_atom[0]} (distance: {dist:.3f} Å)")

        # Glucosepane LGX S1 - AGS R1 (the inter-marker R1 is hosted as AGS SC5).
        for lgx_atom, ags_atom, dist in self._match_nearest_pairs(
            lgx_sc1_atoms, ags_sc5_atoms, self.crosslink_thresholds['LGX_AGS']
        ):
            self.crosslink_pairs.append((lgx_atom, ags_atom))
            LOG.info(
                f" Added LGX-AGS glucosepane pair: atoms {lgx_atom[0]} - "
                f"{ags_atom[0]} (distance: {dist:.3f} Å)"
            )

        # Glucosepane is represented by two marker residues but is a single
        # chemical crosslink.  Never allow an orphan marker or a failed
        # distance match to degrade silently into two ordinary side chains.
        glucosepane_pairs = [
            pair for pair in self.crosslink_pairs
            if pair[0][1:3] == ['LGX', 'SC1']
            and pair[1][1:3] == ['AGS', 'SC5']
        ]
        if lgx_sc1_atoms or ags_sc5_atoms:
            if not (
                len(lgx_sc1_atoms)
                == len(ags_sc5_atoms)
                == len(glucosepane_pairs)
            ):
                raise IncompleteCrosslinkError(
                    "Incomplete glucosepane mapping: expected a one-to-one "
                    "LGX:SC1--AGS:SC5 match for every marker, but found "
                    f"LGX={len(lgx_sc1_atoms)}, AGS={len(ags_sc5_atoms)}, "
                    f"matched={len(glucosepane_pairs)} within "
                    f"{self.crosslink_thresholds['LGX_AGS']:.1f} Å."
                )

        # MOLD R2 is an inter-marker bead hosted on LZS:SC3; its bridge to
        # the second lysine arm (LZD:SC1/S2) completes the CG crosslink.
        for lzs_atom, lzd_atom, dist in self._match_nearest_pairs(
            lzs_sc3_atoms, lzd_sc1_atoms, self.crosslink_thresholds['LZS_LZD']
        ):
            self.crosslink_pairs.append((lzs_atom, lzd_atom))
            LOG.info(
                f" Added LZS-LZD MOLD pair: atoms {lzs_atom[0]} - "
                f"{lzd_atom[0]} (distance: {dist:.3f} Å)"
            )

        mold_pairs = [
            pair for pair in self.crosslink_pairs
            if pair[0][1:3] == ['LZS', 'SC3']
            and pair[1][1:3] == ['LZD', 'SC1']
        ]
        if lzs_sc3_atoms or lzd_sc1_atoms:
            if not (
                len(lzs_sc3_atoms)
                == len(lzd_sc1_atoms)
                == len(mold_pairs)
            ):
                raise IncompleteCrosslinkError(
                    "Incomplete MOLD mapping: expected a one-to-one "
                    "LZS:SC3--LZD:SC1 match for every marker, but found "
                    f"LZS={len(lzs_sc3_atoms)}, LZD={len(lzd_sc1_atoms)}, "
                    f"matched={len(mold_pairs)} within "
                    f"{self.crosslink_thresholds['LZS_LZD']:.1f} Å."
                )

        return self.crosslink_pairs

    def set_crosslink_bonded(self, cnt_model: Optional[int] = None,
                           crosslink_connect: Optional[List[List[List[Any]]]] = None) -> Dict[str, List[List[Any]]]:
        """
        Setup topology for crosslink bonded parameters: bonds, angles and dihedrals.

        Now uses the improved pair-finding algorithm for better crosslink detection.
        """
        LOG.debug(f"Setting up crosslink bonded parameters for model {cnt_model}")

        self.crosslink_bonded = {'bonds': [], 'angles': [], 'dihedrals': []}

        self.find_crosslink_pairs(cnt_model=cnt_model)

        if not self.crosslink_pairs:
            LOG.debug("No crosslink pairs found, returning empty parameters")
            return self.crosslink_bonded

        connections_found = 0

        try:
            for clx, cly in self.crosslink_pairs:
                dist = np.linalg.norm(np.array(clx[-3:]) - np.array(cly[-3:]))

                LOG.debug(f"Processing crosslink between {clx[1]}{clx[2]} and {cly[1]}{cly[2]}: {dist:.3f} Å")

                # LYX SC4 - LY2 SC1 crosslinks
                if (clx[1] == 'LYX' and clx[2] == 'SC4' and
                    cly[1] == 'LY2' and cly[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dlyxly2, f"{self.klyxly2}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])+1), clx[0], cly[0], '1', self.al2yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al2yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al2yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f"Added LYX-LY2 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")

                # LYX SC5 - LY2 SC1 crosslinks (conserve-ring bond)
                elif (clx[1] == 'LYX' and clx[2] == 'SC5' and
                      cly[1] == 'LY2' and cly[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dlyx5ly2, f"{self.klyx5ly2}\n"
                    ])
                    # TODO: no angle terms defined for the SC5-LY2 bond yet.
                    connections_found += 1
                    LOG.info(f"Added LYX-LY2 (SC5) crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")

                # LYX SC5 - LY3 SC1 crosslinks
                elif (clx[1] == 'LYX' and clx[2] == 'SC5' and
                      cly[1] == 'LY3' and cly[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dlyxly3, f"{self.klyxly3}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al3yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-2), clx[0], cly[0], '1', self.al3yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al3yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY3 crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")

                # L4Y SC1 - L5Y SC2 crosslinks
                elif (clx[1] == 'L4Y' and clx[2] == 'SC1' and
                      cly[1] == 'L5Y' and cly[2] == 'SC2'):

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dly45, f"{self.kly45}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        clx[0], cly[0], str(int(cly[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(clx[0])-1), clx[0], cly[0], '1', self.al45y_2, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added L4Y-L5Y crosslink between {clx[0]} and {cly[0]} (distance: {dist:.3f} Å)")

                # Glucosepane: S1 (LGX:SC1) -- R1 (AGS:SC5), plus the
                # angles that cross the marker boundary (B1--S1--R1,
                # S1--R1--R2 with R2 = AGS:SC6 and, in the 10 fs set only,
                # S1--R1--R4).
                elif (clx[1] == 'LGX' and clx[2] == 'SC1' and
                      cly[1] == 'AGS' and cly[2] == 'SC5'):
                    lgx_bb = self._same_residue_bead(clx, 'BB')
                    ags_r2 = self._same_residue_bead(cly, 'SC6')
                    ags_r4 = self._same_residue_bead(cly, 'SC4')
                    if not all((lgx_bb, ags_r2)) or (
                        self.agcp_s1_r1_r4 is not None and not ags_r4
                    ):
                        LOG.error(
                            "Cannot complete glucosepane cross-marker terms: expected "
                            "LGX:BB and AGS:SC4/SC6 in the paired marker residues."
                        )
                        continue

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dgcp_s1_r1,
                        f"{self.kgcp_s1_r1}\n"
                    ])
                    self.crosslink_bonded['angles'].extend([
                        [lgx_bb, clx[0], cly[0], '1',
                         self.agcp_b1_s1_r1[0], f"{self.agcp_b1_s1_r1[1]}\n"],
                        [clx[0], cly[0], ags_r2, '1',
                         self.agcp_s1_r1_r2[0], f"{self.agcp_s1_r1_r2[1]}\n"],
                    ])
                    if self.agcp_s1_r1_r4 is not None:
                        self.crosslink_bonded['angles'].append(
                            [clx[0], cly[0], ags_r4, '1',
                             self.agcp_s1_r1_r4[0], f"{self.agcp_s1_r1_r4[1]}\n"]
                        )
                    connections_found += 1
                    LOG.info(
                        f" Added LGX-AGS glucosepane topology between {clx[0]} "
                        f"and {cly[0]} (distance: {dist:.3f} Å)"
                    )

                # MOLD: R2/IM2 (LZS:SC3) -- S2/L1 (LZD:SC1), plus the
                # three source angles crossing the marker boundary.
                elif (clx[1] == 'LZS' and clx[2] == 'SC3' and
                      cly[1] == 'LZD' and cly[2] == 'SC1'):
                    lzs_r1 = self._same_residue_bead(clx, 'SC2')
                    lzs_r3 = self._same_residue_bead(clx, 'SC4')
                    lzd_bb = self._same_residue_bead(cly, 'BB')
                    if not all((lzs_r1, lzs_r3, lzd_bb)):
                        LOG.error(
                            "Cannot complete MOLD terms: expected LZS:SC2/SC4 "
                            "and LZD:BB in the paired marker residues."
                        )
                        continue

                    self.crosslink_bonded['bonds'].append([
                        clx[0], cly[0], '1', self.dmold_r2_s2,
                        f"{self.kmold_r2_s2}\n"
                    ])
                    self.crosslink_bonded['angles'].extend([
                        [lzs_r1, clx[0], cly[0], '1',
                         self.amold_r1_r2_s2[0], f"{self.amold_r1_r2_s2[1]}\n"],
                        [lzs_r3, clx[0], cly[0], '1',
                         self.amold_r3_r2_s2[0], f"{self.amold_r3_r2_s2[1]}\n"],
                        [clx[0], cly[0], lzd_bb, '1',
                         self.amold_r2_s2_b2[0], f"{self.amold_r2_s2_b2[1]}\n"],
                    ])
                    connections_found += 1
                    LOG.info(
                        f" Added LZS-LZD MOLD topology between {clx[0]} and "
                        f"{cly[0]} (distance: {dist:.3f} Å)"
                    )

                # NOTE: find_crosslink_pairs always emits pairs as (LYX/L4Y, arm),
                # so the reverse-order branches below are currently unreachable; they
                # are kept as a safeguard in case the pair ordering ever changes.
                elif (cly[1] == 'LYX' and cly[2] == 'SC4' and
                      clx[1] == 'LY2' and clx[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dlyxly2, f"{self.klyxly2}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])+1), cly[0], clx[0], '1', self.al2yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al2yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al2yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY2 crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")

                elif (cly[1] == 'LYX' and cly[2] == 'SC5' and
                      clx[1] == 'LY2' and clx[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dlyx5ly2, f"{self.klyx5ly2}\n"
                    ])
                    # TODO: no angle terms defined for the SC5-LY2 bond yet.
                    connections_found += 1
                    LOG.info(f" Added LYX-LY2 (SC5) crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")

                elif (cly[1] == 'LYX' and cly[2] == 'SC5' and
                      clx[1] == 'LY3' and clx[2] == 'SC1'):

                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dlyxly3, f"{self.klyxly3}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al3yx_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-2), cly[0], clx[0], '1', self.al3yx_2, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al3yx_3, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added LYX-LY3 crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")

                elif (cly[1] == 'L4Y' and cly[2] == 'SC1' and
                      clx[1] == 'L5Y' and clx[2] == 'SC2'):

                    self.crosslink_bonded['bonds'].append([
                        cly[0], clx[0], '1', self.dly45, f"{self.kly45}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        cly[0], clx[0], str(int(clx[0])-1), '1', self.al45y_1, f"{self.k_angle}\n"
                    ])
                    self.crosslink_bonded['angles'].append([
                        str(int(cly[0])-1), cly[0], clx[0], '1', self.al45y_2, f"{self.k_angle}\n"
                    ])
                    connections_found += 1
                    LOG.info(f" Added L4Y-L5Y crosslink between {cly[0]} and {clx[0]} (distance: {dist:.3f} Å)")
                else:
                    LOG.debug(f"    Distance {dist:.3f} Å between {clx[1]}{clx[2]} and {cly[1]}{cly[2]} - not a recognized crosslink type")

            LOG.info(f"Created {len(self.crosslink_bonded['bonds'])} bonds and "
                    f"{len(self.crosslink_bonded['angles'])} angles from "
                    f"{connections_found} crosslink connections")

        except Exception as e:
            LOG.error(f"Error setting up crosslink bonded parameters: {str(e)}")

        glucosepane_pair_count = sum(
            1 for clx, cly in self.crosslink_pairs
            if clx[1:3] == ['LGX', 'SC1']
            and cly[1:3] == ['AGS', 'SC5']
        )
        if glucosepane_pair_count:
            glucosepane_bond_count = sum(
                1 for bond in self.crosslink_bonded['bonds']
                if bond[3:5] == [self.dgcp_s1_r1, f"{self.kgcp_s1_r1}\n"]
            )
            expected_angles = {
                (*angle,) for angle in (
                    self.agcp_b1_s1_r1, self.agcp_s1_r1_r2, self.agcp_s1_r1_r4
                ) if angle is not None
            }
            glucosepane_angle_count = sum(
                1 for angle in self.crosslink_bonded['angles']
                if (angle[4], angle[5].rstrip("\n")) in expected_angles
            )
            if (
                glucosepane_bond_count != glucosepane_pair_count
                or glucosepane_angle_count
                != len(expected_angles) * glucosepane_pair_count
            ):
                raise IncompleteCrosslinkError(
                    "Incomplete glucosepane bonded topology: expected "
                    f"{glucosepane_pair_count} inter-marker bond(s) and "
                    f"{len(expected_angles) * glucosepane_pair_count} angle(s), generated "
                    f"{glucosepane_bond_count} bond(s) and "
                    f"{glucosepane_angle_count} angle(s)."
                )

        mold_pair_count = sum(
            1 for clx, cly in self.crosslink_pairs
            if clx[1:3] == ['LZS', 'SC3']
            and cly[1:3] == ['LZD', 'SC1']
        )
        if mold_pair_count:
            mold_bond_count = sum(
                1 for bond in self.crosslink_bonded['bonds']
                if bond[3:5] == [self.dmold_r2_s2, f"{self.kmold_r2_s2}\n"]
            )
            expected_angles = {
                (*self.amold_r1_r2_s2,),
                (*self.amold_r3_r2_s2,),
                (*self.amold_r2_s2_b2,),
            }
            mold_angle_count = sum(
                1 for angle in self.crosslink_bonded['angles']
                if (angle[4], angle[5].rstrip("\n")) in expected_angles
            )
            if (
                mold_bond_count != mold_pair_count
                or mold_angle_count != 3 * mold_pair_count
            ):
                raise IncompleteCrosslinkError(
                    "Incomplete MOLD bonded topology: expected "
                    f"{mold_pair_count} inter-marker bond(s) and "
                    f"{3 * mold_pair_count} angle(s), generated "
                    f"{mold_bond_count} bond(s) and "
                    f"{mold_angle_count} angle(s)."
                )

        return self.crosslink_bonded
