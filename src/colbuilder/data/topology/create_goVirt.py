#!/usr/bin/env python
# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import argparse
import subprocess
import numpy as np
from typing import List, Tuple, Dict, Optional, Any
import logging
import sys
import os
import traceback

# Setup logging with more detailed format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
LOG = logging.getLogger(__name__)

def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate Go-like model for protein structure",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-s', '--structure', 
        required=True, 
        help='File containing the coarse-grained structure of the protein in PDB format.'
    )
    parser.add_argument(
        '-f', '--contacts', 
        required=True, 
        help='File containing the contact analysis of the (atomistic) protein structure.'
    )
    parser.add_argument(
        '--moltype', 
        default='molecule_0', 
        help='Molecule name used as prefix in output files and virtual bead names.'
    )
    parser.add_argument(
        '--go_eps', 
        type=float, 
        default=9.414, 
        help='Dissociation energy [kJ/mol] of the Lennard-Jones potential.'
    )
    parser.add_argument(
        '--cutoff_short', 
        type=float, 
        default=0.3, 
        help='Lower cutoff distance [nm] for Go-like interactions.'
    )
    parser.add_argument(
        '--cutoff_long', 
        type=float, 
        default=1.1, 
        help='Upper cutoff distance [nm] for Go-like interactions.'
    )
    parser.add_argument(
        '--Natoms', 
        type=int, 
        help='Number of coarse-grained beads in the protein excluding virtual Go beads.'
    )
    parser.add_argument(
        '--debug', 
        action='store_true', 
        help='Enable debug logging'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.structure):
        parser.error(f"Structure file not found: {args.structure}")
    if not os.path.exists(args.contacts):
        parser.error(f"Contacts file not found: {args.contacts}")
    
    # Set debug logging if requested
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        LOG.debug("Debug logging enabled")
    
    return args

def get_settings() -> Tuple[str, str, str, int, int, List[int], int, int, int]:
    """
    Get default settings for the Go-model generation.
    
    Returns:
        Tuple containing file names, header lines, sequence distance,
        column indices, missing residues, missing atoms, and force field type.
    """
    # Settings for processing contact maps and structure files
    file_BB = 'BB.pdb'           # Temporary file for backbone beads
    file_OV = 'OV.map'           # Temporary file for overlap contacts
    file_rCSU = 'rCSU.map'       # Temporary file for residue contacts
    header_lines = 0             # Number of header lines to skip in contact files
    seqDist = 4                  # Minimum sequence distance for contacts
    cols = [5, 9, 10]            # Columns to extract from contact files
    missRes = 0                  # Missing residue offset
    missAt = 0                   # Missing atom offset
    c6c12 = 0                    # Flag for C6/C12 parameter format (0=sigma/epsilon, 1=C6/C12)
    
    LOG.debug(f"Using settings: header_lines={header_lines}, seqDist={seqDist}, cols={cols}, missRes={missRes}, missAt={missAt}, c6c12={c6c12}")
    
    return (file_BB, file_OV, file_rCSU, header_lines, seqDist, cols, missRes, missAt, c6c12)

def run_command(cmd: str, error_msg: str) -> None:
    """
    Run a shell command safely with proper error handling.
    
    Args:
        cmd: Command to run.
        error_msg: Error message to display if command fails.
        
    Raises:
        subprocess.CalledProcessError: If command fails.
    """
    try:
        LOG.debug(f"        Running command: {cmd}")
        result = subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True
        )
        
        if result.stdout:
            LOG.debug(f"    Command output: {result.stdout.strip()}")
        
        return result
    except subprocess.CalledProcessError as e:
        LOG.error(f"{error_msg}: {e}")
        LOG.error(f"Command stderr: {e.stderr}")
        raise

def read_data(struct_pdb: str, file_contacts: str, file_BB: str, file_OV: str, file_rCSU: str, 
              header_lines: int, cols: List[int], Natoms: Optional[int]) -> Tuple[np.ndarray, List[str], List[List[float]], int]:
    """
    Read and process structure and contact data.
    
    Args:
        struct_pdb: Path to structure PDB file.
        file_contacts: Path to contacts file.
        file_BB: Path to backbone temporary file.
        file_OV: Path to overlap contact map temporary file.
        file_rCSU: Path to residue contact map temporary file.
        header_lines: Number of header lines to skip.
        cols: Column indices to extract from contact files.
        Natoms: Number of atoms in the structure (if known).
        
    Returns:
        Tuple containing backbone indices and coordinates, amino acid names,
        contact map data, and total atom count.
        
    Raises:
        subprocess.CalledProcessError: If shell commands fail.
        IOError: If file operations fail.
        ValueError: If data format is incorrect.
    """
    LOG.debug(f"Reading structure from {struct_pdb} and contacts from {file_contacts}")
    
    try:
        # Extract overlap contacts
        run_command(
            f"grep '1 [01] [01] [01]' {file_contacts} > {file_OV}",
            "Failed to extract overlap contacts"
        )
        run_command(f"echo '' >> {file_OV}", "Failed to append to overlap file")
        
        # Extract residue contacts
        run_command(
            f"grep '0 [01] [01] 1' {file_contacts} > {file_rCSU}",
            "Failed to extract residue contacts"
        )
        run_command(f"echo '' >> {file_rCSU}", "Failed to append to residue file")
        
        # Extract backbone beads
        run_command(
            f"grep 'BB' {struct_pdb} > {file_BB}",
            "Failed to extract backbone beads"
        )
        run_command(f"echo '' >> {file_BB}", "Failed to append to backbone file")
        
        # Read backbone data
        indBB = []
        nameAA = []
        try:
            with open(file_BB, 'r') as fid:
                dat = fid.readlines()[header_lines:-1]
            
            LOG.debug(f'Found {len(dat)} coarse-grained BB beads in the protein')
            
            for line in dat:
                try:
                    atom_num = int(line[6:11].strip())
                    x_coord = float(line[26:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:55].strip())
                    residue = line[17:20].strip()
                    
                    indBB.append([atom_num, x_coord, y_coord, z_coord])
                    nameAA.append(residue)
                except (ValueError, IndexError) as e:
                    LOG.warning(f"Skipping malformed PDB line: {line.strip()} - {e}")
                    
            indBB = np.array(indBB)
            
            if len(indBB) == 0:
                raise ValueError("No valid BB beads found in structure")
                
        except IOError as e:
            LOG.error(f"Error reading backbone data: {e}")
            raise
            
        # Read contact map data
        map_OVrCSU = []
        
        for filename in [file_OV, file_rCSU]:
            try:
                with open(filename, 'r') as fid:
                    dat = fid.readlines()[header_lines:-1]
                
                LOG.debug(f'Found {len(dat)} contacts in {filename}')
                
                for line in dat:
                    try:
                        tmp = line.replace('\t', ' ').split()
                        # Extract columns specified by indices
                        contact_data = [float(tmp[col]) for col in cols]
                        map_OVrCSU.append(contact_data)
                    except (ValueError, IndexError) as e:
                        LOG.warning(f"Skipping malformed contact line in {filename}: {line.strip()} - {e}")
            except IOError as e:
                LOG.error(f"Error reading contact data from {filename}: {e}")
                raise
        
        # Use provided Natoms if available, otherwise use the number of atoms from BB data
        if Natoms is None:
            Natoms = len(indBB)
            LOG.debug(f"Using {Natoms} atoms from structure (no --Natoms provided)")
        else:
            LOG.debug(f"Using user-provided Natoms: {Natoms}")
        
        return indBB, nameAA, map_OVrCSU, Natoms
        
    except (subprocess.CalledProcessError, IOError) as e:
        LOG.error(f"Failed to read data: {e}")
        LOG.debug(traceback.format_exc())
        raise
    except Exception as e:
        LOG.error(f"Unexpected error while reading data: {e}")
        LOG.debug(traceback.format_exc())
        raise ValueError(f"Error processing input data: {e}")

def get_go(indBB: np.ndarray, nameAA: List[str], map_OVrCSU: List[List[float]], cutoff_short: float, 
           cutoff_long: float, go_eps: float, seqDist: int, missRes: int) -> List[List[float]]:
    """
    Calculate Go-like interaction parameters.
    
    Args:
        indBB: Array of backbone indices and coordinates.
        nameAA: List of amino acid names.
        map_OVrCSU: List of contact map data.
        cutoff_short: Lower cutoff distance for Go-like interactions.
        cutoff_long: Upper cutoff distance for Go-like interactions.
        go_eps: Dissociation energy of the Lennard-Jones potential.
        seqDist: Minimum sequence distance for contacts.
        missRes: Missing residue offset.
        
    Returns:
        List of symmetric pairs with Go-like interaction parameters.
        
    Raises:
        ValueError: If calculation fails.
    """
    LOG.debug(f"Calculating Go-like interactions using cutoffs {cutoff_short}-{cutoff_long} nm and epsilon {go_eps} kJ/mol")
    
    try:
        # Calculate distances for all contacts
        valid_contacts = 0
        invalid_indices = 0
        
        for k in range(len(map_OVrCSU)):
            try:
                # Get indices with offset
                i = int(map_OVrCSU[k][0]) - missRes - 1
                j = int(map_OVrCSU[k][1]) - missRes - 1
                
                # Check index bounds
                if i < 0 or i >= len(indBB) or j < 0 or j >= len(indBB):
                    LOG.warning(f"Contact indices out of range: {i+missRes+1}-{j+missRes+1}")
                    invalid_indices += 1
                    continue
                
                # Calculate distance vector
                dist_vec = indBB[j, 1:4] - indBB[i, 1:4]
                # Calculate and store distance in nm (convert from Å)
                dist = np.linalg.norm(dist_vec) / 10.0
                map_OVrCSU[k][2] = dist
                valid_contacts += 1
                
            except (ValueError, IndexError) as e:
                LOG.warning(f"Error calculating distance for contact {k}: {e}")
        
        LOG.debug(f"Calculated distances for {valid_contacts} contacts, skipped {invalid_indices} with invalid indices")
        
        # Apply distance and sequence separation filters
        pairs = []
        short_contacts = 0
        long_contacts = 0
        seq_filtered = 0
        
        for k in range(len(map_OVrCSU)):
            try:
                # Get residue indices
                res_i = map_OVrCSU[k][0]
                res_j = map_OVrCSU[k][1]
                dist = map_OVrCSU[k][2]
                
                # Apply filters
                if dist <= cutoff_short:
                    short_contacts += 1
                    continue
                    
                if dist >= cutoff_long:
                    long_contacts += 1
                    continue
                    
                if abs(res_j - res_i) < seqDist:
                    seq_filtered += 1
                    continue
                
                # Calculate LJ parameters
                sigma = dist / 1.12246204830  # 2^(1/6) = position of LJ minimum
                Vii = 4.0 * pow(sigma, 6) * go_eps      # C6 parameter
                Wii = 4.0 * pow(sigma, 12) * go_eps     # C12 parameter
                
                # Get atom indices
                atom_i = int(indBB[int(res_i) - missRes - 1, 0])
                atom_j = int(indBB[int(res_j) - missRes - 1, 0])
                
                # Store the pair with all parameters
                pairs.append([
                    atom_i, atom_j, Vii, Wii, res_i, res_j, dist, sigma
                ])
                
            except (ValueError, IndexError) as e:
                LOG.warning(f"Error processing contact for Go calculation: {e}")
        
        LOG.debug(f"Applied filters: removed {short_contacts} short contacts (<{cutoff_short} nm), " +
                f"{long_contacts} long contacts (>{cutoff_long} nm), and {seq_filtered} close in sequence (<{seqDist} apart)")
        LOG.debug(f"Generated {len(pairs)} Go-like interaction pairs")
        
        # Find symmetric pairs (both i-j and j-i exist)
        sym_pairs = []
        for k in range(len(pairs)):
            # Check if symmetric pair exists
            is_symmetric = False
            for l in range(k+1, len(pairs)):
                if pairs[k][0] == pairs[l][1] and pairs[k][1] == pairs[l][0]:
                    is_symmetric = True
                    break
                    
            # Only keep i-j where i < j and j-i also exists
            if pairs[k][0] < pairs[k][1] and is_symmetric:
                sym_pairs.append(pairs[k])
        
        LOG.debug(f"Found {len(sym_pairs)} symmetric Go-like interactions out of {len(pairs)} total")
        
        return sym_pairs
        
    except Exception as e:
        LOG.error(f"Error calculating Go-like interactions: {e}")
        LOG.debug(traceback.format_exc())
        raise ValueError(f"Go calculation failed: {e}")

def write_files(file_pref: str, sym_pairs: List[List[float]], missAt: int, indBB: np.ndarray, 
                missRes: int, Natoms: int, nameAA: List[str], go_eps: float, c6c12: int) -> None:
    """
    Write output files for the Go-like model.
    
    Args:
        file_pref: Prefix for output files.
        sym_pairs: List of symmetric pairs with Go-like parameters.
        missAt: Missing atom offset.
        indBB: Array of backbone indices and coordinates.
        missRes: Missing residue offset.
        Natoms: Number of atoms in the structure.
        nameAA: List of amino acid names.
        go_eps: Dissociation energy of the Lennard-Jones potential.
        c6c12: Flag for C6/C12 parameter format.
        
    Raises:
        IOError: If file operations fail.
        subprocess.CalledProcessError: If shell commands fail.
    """
    LOG.debug(f"Writing Go-like model files with prefix '{file_pref}'")
    
    try:
        # Write table of non-bonded parameters
        table_file = f'{file_pref}_go-table.itp'
        LOG.debug(f"Writing non-bonded parameters to {table_file}")
        
        with open(table_file, 'w') as f:
            f.write('[ nonbond_params ] \n')
            f.write('; i   j   func     V  W \n')
            f.write('; OV + symmetric rCSU contacts \n')
            
            for k in range(len(sym_pairs)):
                # Format depends on c6c12 flag
                if c6c12 == 1:
                    # Use C6/C12 format
                    line = (f" {file_pref}_{int(sym_pairs[k][4])}  {file_pref}_{int(sym_pairs[k][5])}    1  "
                          f"{sym_pairs[k][2]:.10f}  {sym_pairs[k][3]:.10f}  ;  "
                          f"{int(sym_pairs[k][0])+missAt}  {int(sym_pairs[k][1])+missAt}  {sym_pairs[k][6]:.3f} \n")
                else:
                    # Use sigma/epsilon format
                    line = (f" {file_pref}_{int(sym_pairs[k][4])}  {file_pref}_{int(sym_pairs[k][5])}    1  "
                          f"{sym_pairs[k][7]:.10f}  {go_eps:.10f}  ;  "
                          f"{int(sym_pairs[k][0])+missAt}  {int(sym_pairs[k][1])+missAt}  {sym_pairs[k][6]:.3f} \n")
                f.write(line)
        
        # Include file in main go-table.itp
        run_command(
            f"echo '#include \"{table_file}\"' >> go-table.itp ", 
            "Failed to update go-table.itp"
        )
        
        # Write atom types for virtual sites
        sites_file = f'{file_pref}_go-sites.itp'
        LOG.debug(f"Writing virtual site definitions to {sites_file}")
        
        with open(sites_file, 'w') as f:
            f.write('[ atomtypes ] \n')
            f.write('; protein BB virtual particles \n')
            for k in range(len(indBB)):
                f.write(f"{file_pref}_{k+1+missRes} 0.0 0.000 A 0.0 0.0 \n")
        
        # Include file in main go-sites.itp
        run_command(
            f"echo '#include \"{sites_file}\"' >> go-sites.itp ", 
            "Failed to update go-sites.itp"
        )
        
        # Write exclusions
        excl_file = f'{file_pref}_go-excl.itp'
        LOG.debug(f"Writing exclusions to {excl_file}")
        
        with open(excl_file, 'w') as f:
            f.write(';[ exclusions ] \n')
            f.write('; OV + symmetric rCSU contacts \n')
            for k in range(len(sym_pairs)):
                f.write(f" {int(sym_pairs[k][0])+missAt}  {int(sym_pairs[k][1])+missAt}  \t ;  "
                      f"{int(sym_pairs[k][4])}  {int(sym_pairs[k][5])} \n")
        
        # Write harmonic bonds
        harm_file = f'{file_pref}_go-harm.itp'
        LOG.debug(f"Writing harmonic bonds to {harm_file}")
        
        with open(harm_file, 'w') as f:
            f.write('; Go bonds as harmonic bonds between the virtual particles: \n')
            f.write('; OV + symmetric rCSU contacts \n')
            
            # Write bonds for pairs
            for k in range(len(sym_pairs)):
                f.write(f" {int(sym_pairs[k][4]+Natoms)}  {int(sym_pairs[k][5]+Natoms)}  1  "
                      f"{sym_pairs[k][6]:.3f}  1250  ; {file_pref}_{int(sym_pairs[k][4])}  "
                      f"{file_pref}_{int(sym_pairs[k][5])} \n")
            
            # Add bonds for lonely beads to help visualization
            sym_pairs_array = np.array(sym_pairs)
            lonely_beads = 0
            
            for k in range(len(indBB)):
                # Check if bead participates in any Go contacts
                has_contacts = False
                
                if len(sym_pairs) > 0:  # Only check if we have pairs
                    if (np.sum(sym_pairs_array[:, 4] == k+1) > 0 or 
                        np.sum(sym_pairs_array[:, 5] == k+1) > 0):
                        has_contacts = True
                
                # Add dummy bond for visualization if no contacts
                if not has_contacts:
                    f.write(f" {k+1+Natoms}  {k+Natoms}  1  1.  1     ; "
                          f"{file_pref}_{k+1}  {file_pref}_{k} --> added for vmd \n")
                    lonely_beads += 1
            
            if lonely_beads > 0:
                LOG.debug(f"Added {lonely_beads} dummy bonds for lonely beads (for visualization)")

        LOG.debug("Successfully wrote all Go-like model files!")
        
    except IOError as e:
        LOG.error(f"IOError while writing files: {e}")
        LOG.debug(traceback.format_exc())
        raise
    except subprocess.CalledProcessError as e:
        LOG.error(f"Subprocess error while writing files: {e}")
        LOG.debug(traceback.format_exc())
        raise
    except Exception as e:
        LOG.error(f"Unexpected error while writing files: {e}")
        LOG.debug(traceback.format_exc())
        raise ValueError(f"Failed to write Go-like model files: {e}")

def main() -> None:
    """
    Main function that orchestrates the Go-like model generation.
    """
    try:
        args = parse_arguments()
        
        LOG.debug(f"Starting Go-like model generation for {args.moltype}")
        LOG.debug(f"Using structure file: {args.structure}")
        LOG.debug(f"Using contacts file: {args.contacts}")
        LOG.debug(f"Go-model parameters: eps={args.go_eps} kJ/mol, cutoffs={args.cutoff_short}-{args.cutoff_long} nm")
        
        # Get default settings
        file_BB, file_OV, file_rCSU, header_lines, seqDist, cols, missRes, missAt, c6c12 = get_settings()

        # Read structure and contacts data
        indBB, nameAA, map_OVrCSU, Natoms = read_data(
            args.structure, args.contacts, file_BB, file_OV, file_rCSU, 
            header_lines, cols, args.Natoms
        )

        # Calculate Go-like interactions
        sym_pairs = get_go(
            indBB, nameAA, map_OVrCSU, args.cutoff_short, args.cutoff_long, 
            args.go_eps, seqDist, missRes
        )

        # Write output files
        write_files(
            args.moltype, sym_pairs, missAt, indBB, missRes, 
            Natoms, nameAA, args.go_eps, c6c12
        )
        
        # Clean up temporary files
        for temp_file in [file_BB, file_OV, file_rCSU]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                LOG.debug(f"Removed temporary file: {temp_file}")
        
        LOG.debug("Go-like model generation completed successfully!")
        
    except Exception as e:
        LOG.error(f"Go-like model generation failed: {e}")
        LOG.debug(traceback.format_exc())
        sys.exit(1)

if __name__ == '__main__':
    main()