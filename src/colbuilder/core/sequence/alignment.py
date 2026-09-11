"""
This module handles the alignment of input sequences with template sequences, specifically for
collagen triple helices. It includes functionality for processing hydroxyprolines, staggering
sequences, and generating outputs compatible with MUSCLE® and MODELLER®.

Key Features:
--------------
1. **Sequence Alignment**:
   - Align input sequences amongst themselves using MUSCLE.
   - Align input multiple sequence alignment (MSA) with template sequences.

2. **Hydroxyproline Processing**:
   - Identify and process hydroxyprolines (OLC: O) in sequences.
   - Restore hydroxyprolines in aligned sequences using position mapping.

3. **Staggering and Formatting**:
   - Add stagger to sequences based on collagen register configuration (A: 0, B: +1, C: +2).
   - Equalize sequence lengths by padding with gaps ('-').

4. **Output Generation**:
   - Generate MSA files in FASTA format.
   - Generate MODELLER-compatible alignment files for structure modeling.

5. **MUSCLE Integration**:
   - Automatically detect MUSCLE version and execute alignment commands.
   - Handle errors and ensure compatibility with different MUSCLE versions.

Usage:
------
This module is designed to be used as part of a pipeline for aligning collagen sequences.
The main entry point is the `align_sequences` function, which performs all alignment steps
and outputs the final MSA and MODELLER-formatted files.

Example:
--------
```python
from pathlib import Path
from colbuilder.core.sequence.alignment import align_sequences

# Define input files
input_fasta = Path("input_sequences.fasta")
template_fasta = Path("template_sequences.fasta")
output_prefix = "aligned_output"
template_pdb = Path("template_structure.pdb")

# Perform alignment
msa_output, modeller_output = align_sequences(input_fasta, template_fasta,
                                output_prefix, template_pdb)

print(f"MSA file saved to: {msa_output}")
print(f"MODELLER alignment file saved to: {modeller_output}")
```
"""

# Copyright (c) 2024, ColBuilder Development Team
# Distributed under the terms of the Apache License 2.0

import os
import shutil
import subprocess
import tempfile
from typing import List, Dict, Tuple
import typing as t
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

from colbuilder.core.utils.dec import timeit
from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


class Alignment:
    """
    Attributes:
        input_fasta (Path): Path to the input FASTA file.
        template_fasta (Path): Path to the template FASTA file.
        output_prefix (str): Prefix for output files.
        template_pdb (Path): Path to the template PDB file.
        input_sequences (List[SeqRecord]): List of input sequences.
        template_sequences (List[SeqRecord]): List of template sequences.
        hydroxyproline_positions (Dict): Dictionary to store hydroxyproline positions.
        temp_dir (str): Path to temporary directory for intermediate files.
    """

    def __init__(
        self,
        input_fasta: Path,
        template_fasta: Path,
        output_prefix: str,
        template_pdb: Path,
    ):
        """
        Initialize the Alignment object.

        Args:
            input_fasta (Path): Path to the input FASTA file.
            template_fasta (Path): Path to the template FASTA file.
            output_prefix (str): Prefix for output files.
            template_pdb (Path): Path to the template PDB file.

        Raises:
            FileNotFoundError: If any of the input files are not found.
            ValueError: If the input files are empty or in an incorrect format.
        """
        self.input_fasta = input_fasta
        self.template_fasta = template_fasta
        self.output_prefix = output_prefix
        self.template_pdb = template_pdb

        try:
            self.input_sequences = list(SeqIO.parse(input_fasta, "fasta"))
            self.template_sequences = list(SeqIO.parse(template_fasta, "fasta"))
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Input file not found: {e.filename}") from e
        except ValueError as e:
            raise ValueError(f"Error parsing input files: {str(e)}") from e

        if not self.input_sequences or not self.template_sequences:
            raise ValueError("Input or template sequences are empty")

        self.hydroxyproline_positions: Dict[str, Dict[str, List[int]]] = {
            "template": {},
            "input": {},
        }
        self.temp_dir = tempfile.mkdtemp()

    def add_stagger_to_sequences(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        """
        Add stagger to sequences based on their ID (register configuration). A: 0, B: +1, C: +2

        Args:
            sequences (List[SeqRecord]): List of sequence records to process.

        Returns:
            List[SeqRecord]: List of sequence records with stagger added.

        Raises:
            ValueError: If a sequence ID is in an unexpected format.
        """

        def calculate_stagger(seq_id: str) -> int:
            try:
                return ord(seq_id.split(":")[1]) - 64
            except (IndexError, ValueError):
                raise ValueError(f"Invalid sequence ID format: {seq_id}")

        return [
            SeqRecord(
                Seq("-" * calculate_stagger(seq.id) + str(seq.seq)),
                id=seq.id,
                description=seq.description,
            )
            for seq in sequences
        ]

    def process_hydroxyprolines(
        self, sequence: str, seq_type: str, chain: str, original_id: str
    ) -> Tuple[str, List[int]]:
        """
        Process hydroxyprolines (OLC: O) in a sequence.

        Args:
            sequence (str): The input sequence.
            seq_type (str): The type of sequence ('input' or 'template').
            chain (str): The chain identifier.
            original_id (str): The original sequence ID.

        Returns:
            Tuple[str, List[int]]: Processed sequence and list of hydroxyproline positions.
        """
        seq_list = list(sequence)
        positions = [i for i, char in enumerate(seq_list) if char.upper() == "O"]
        for i in positions:
            seq_list[i] = "P"
        self.hydroxyproline_positions[seq_type][chain] = {
            "positions": positions,
            "original_id": original_id,
        }
        return "".join(seq_list), positions

    def create_position_mapping(
        self, original_seq: str, aligned_seq: str
    ) -> Dict[int, int]:
        """
        Create a position mapping between the original and aligned sequences.

        Args:
            original_seq (str): The original sequence.
            aligned_seq (str): The aligned sequence.

        Returns:
            Dict[int, int]: Mapping of original positions to aligned positions.
        """
        map_original_to_aligned = {}
        aligned_index = 0
        for original_index, char in enumerate(original_seq):
            while (
                aligned_index < len(aligned_seq) and aligned_seq[aligned_index] == "-"
            ):
                aligned_index += 1
            if aligned_index < len(aligned_seq):
                map_original_to_aligned[original_index] = aligned_index
                aligned_index += 1
        return map_original_to_aligned

    def restore_hydroxyproline_with_mapping(
        self,
        aligned_seq: str,
        original_positions: List[int],
        position_mapping: Dict[int, int],
    ) -> str:
        """
        Restore hydroxyprolines in the aligned sequence using the position mapping.

        Args:
            aligned_seq (str): The aligned sequence.
            original_positions (List[int]): Original positions of hydroxyprolines.
            position_mapping (Dict[int, int]): Mapping of original positions to aligned positions.

        Returns:
            str: Aligned sequence with hydroxyprolines restored.
        """
        aligned_seq_list = list(aligned_seq)
        for original_pos in original_positions:
            if original_pos in position_mapping:
                new_pos = position_mapping[original_pos]
                aligned_seq_list[new_pos] = "O"
        return "".join(aligned_seq_list)

    def get_muscle_version(self) -> str:
        """
        Get the version of the MUSCLE tool.

        Returns:
            str: MUSCLE version.

        Raises:
            RuntimeError: If MUSCLE is not installed or not found in PATH.
        """
        try:
            result = subprocess.run(
                ["muscle", "-version"], capture_output=True, text=True, check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            raise RuntimeError("MUSCLE is not installed or not found in PATH.")

    def align_sequences_with_muscle(self, input_path: str, output_path: str) -> None:
        """
        Align sequences using MUSCLE.

        Args:
            input_path (str): Path to the input FASTA file.
            output_path (str): Path to save the aligned output.

        Raises:
            RuntimeError: If MUSCLE alignment fails.
        """
        muscle_version = self.get_muscle_version()
        if "3.8" in muscle_version:
            muscle_command = ["muscle", "-in", input_path, "-out", output_path]
        else:
            muscle_command = ["muscle", "-align", input_path, "-output", output_path]
        try:
            subprocess.run(
                muscle_command,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            raise RuntimeError("MUSCLE alignment failed.")

    def extract_last_atom_serial_number(self, pdb_file: str) -> str:
        """
        Extract the last atom serial number from a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.

        Returns:
            str: Last atom serial number.
        """
        with open(pdb_file, "r") as file:
            for line in reversed(file.readlines()):
                if line.startswith("ATOM"):
                    return line[6:11].strip()
        return "0"

    def write_modeller_formatted_output(
        self, sequences: List[SeqRecord], output_file: str, end_atom: str
    ) -> None:
        """
        Write sequences in a format suitable for MODELLER.

        Args:
            sequences (List[SeqRecord]): List of sequences to write.
            output_file (str): Path to save the output file.
            end_atom (str): Last atom serial number.
        """
        pdb_file_name_without_ext = os.path.splitext(
            os.path.basename(self.template_pdb)
        )[0]

        template_sequences = [
            seq for seq in sequences if seq.id.startswith("template:")
        ]
        input_sequences = [seq for seq in sequences if seq.id.startswith("input:")]

        with open(output_file, "w") as out_file:
            out_file.write(">P1;template\n")
            out_file.write(
                f"structure:{pdb_file_name_without_ext}:1:A:+{end_atom}:C:::-1.00:-1.00\n"
            )
            out_file.write("/".join(str(seq.seq) for seq in template_sequences) + "*\n")

            out_file.write("\n>P1;target\n")
            out_file.write(f"sequence:{self.output_prefix}:: :: :::-1.00:-1.00\n")
            out_file.write("/".join(str(seq.seq) for seq in input_sequences) + "*\n")

    def equalize_sequence_lengths(self, sequences: List[SeqRecord]) -> List[SeqRecord]:
        """
        Equalize the lengths of sequences by padding with '-'.

        Args:
            sequences (List[SeqRecord]): List of sequences to equalize.

        Returns:
            List[SeqRecord]: List of sequences with equalized lengths.
        """
        max_length = max(len(seq.seq) for seq in sequences)
        return [
            SeqRecord(
                Seq(str(seq.seq).ljust(max_length, "-")),
                id=seq.id,
                description=seq.description,
            )
            for seq in sequences
        ]

    @timeit
    def align_sequences(self) -> Tuple[str, str]:
        """
        Align input sequences with template sequences and restore hydroxyprolines.

        Returns:
            Tuple[List[SeqRecord], str, str]: Aligned and restored sequences,
            MSA output path, and MODELLER output path.

        Raises:
            Exception: If any error occurs during the alignment process.
        """
        try:
            # Step 1: Align input sequences amongst themselves
            input_modified_sequences = [
                SeqRecord(
                    Seq(
                        self.process_hydroxyprolines(
                            str(seq.seq), "input", chr(65 + i), seq.id
                        )[0]
                    ),
                    id=f"input:{chr(65 + i)}",
                    description=seq.description,
                )
                for i, seq in enumerate(self.input_sequences)
            ]

            temp_input_path = os.path.join(
                self.temp_dir, f"{self.output_prefix}_input.fasta"
            )
            temp_input_aligned_path = os.path.join(
                self.temp_dir, f"{self.output_prefix}_input_aligned.afa"
            )

            SeqIO.write(input_modified_sequences, temp_input_path, "fasta")
            self.align_sequences_with_muscle(temp_input_path, temp_input_aligned_path)

            input_aligned_sequences = list(
                SeqIO.parse(temp_input_aligned_path, "fasta")
            )

            # Step 2: Align input MSA with template
            template_modified_sequences = [
                SeqRecord(
                    Seq(
                        self.process_hydroxyprolines(
                            str(seq.seq), "template", chr(65 + i), seq.id
                        )[0]
                    ),
                    id=f"template:{chr(65 + i)}",
                    description=seq.description,
                )
                for i, seq in enumerate(self.template_sequences)
            ]

            combined_sequences = template_modified_sequences + input_aligned_sequences
            temp_combined_path = os.path.join(
                self.temp_dir, f"{self.output_prefix}_combined.fasta"
            )
            temp_final_aligned_path = os.path.join(
                self.temp_dir, f"{self.output_prefix}_final_aligned.afa"
            )

            SeqIO.write(combined_sequences, temp_combined_path, "fasta")
            self.align_sequences_with_muscle(
                temp_combined_path, temp_final_aligned_path
            )

            final_aligned_sequences = list(
                SeqIO.parse(temp_final_aligned_path, "fasta")
            )

            original_order = {seq.id: i for i, seq in enumerate(combined_sequences)}
            final_aligned_sequences.sort(key=lambda seq: original_order[seq.id])

            # Step 3: Restore hydroxyprolines
            restored_sequences = []
            for seq in final_aligned_sequences:
                seq_type, chain = seq.id.split(":")
                if (
                    seq_type in self.hydroxyproline_positions
                    and chain in self.hydroxyproline_positions[seq_type]
                ):
                    original_id = self.hydroxyproline_positions[seq_type][chain][
                        "original_id"
                    ]
                    original_seq = next(
                        s
                        for s in (
                            self.input_sequences
                            if seq_type == "input"
                            else self.template_sequences
                        )
                        if s.id == original_id
                    )
                    position_mapping = self.create_position_mapping(
                        str(original_seq.seq), str(seq.seq)
                    )
                    restored_seq = self.restore_hydroxyproline_with_mapping(
                        str(seq.seq),
                        self.hydroxyproline_positions[seq_type][chain]["positions"],
                        position_mapping,
                    )
                    restored_sequences.append(
                        SeqRecord(
                            Seq(restored_seq), id=seq.id, description=seq.description
                        )
                    )
                else:
                    restored_sequences.append(seq)

            # Step 4: Add stagger to template and input sequences
            template_staggered = self.add_stagger_to_sequences(
                [seq for seq in restored_sequences if seq.id.startswith("template:")]
            )
            input_staggered = self.add_stagger_to_sequences(
                [seq for seq in restored_sequences if seq.id.startswith("input:")]
            )

            template_staggered.sort(key=lambda seq: seq.id.split(":")[1])
            input_staggered.sort(key=lambda seq: seq.id.split(":")[1])

            template_staggered = self.equalize_sequence_lengths(template_staggered)
            input_staggered = self.equalize_sequence_lengths(input_staggered)

            staggered_restored_sequences = template_staggered + input_staggered

            # Step 5: Write output files
            # 5.1 Write MSA
            msa_output = f"{self.output_prefix}_alignment.fasta"
            SeqIO.write(staggered_restored_sequences, msa_output, "fasta")

            # 5.2 Write MODELLER-formatted alignment
            end_atom = self.extract_last_atom_serial_number(self.template_pdb)
            modeller_output = f"{self.output_prefix}_modeller.ali"
            self.write_modeller_formatted_output(
                staggered_restored_sequences, modeller_output, end_atom
            )

            return msa_output, modeller_output

        except Exception as e:
            LOG.error(f"Error in align_sequences: {str(e)}", exc_info=True)
            raise
        finally:
            shutil.rmtree(self.temp_dir, ignore_errors=True)


@timeit
def align_sequences(
    input_fasta_path: Path,
    template_fasta_path: Path,
    output_prefix: str,
    template_pdb: Path,
) -> t.Tuple[str, str]:
    """
    Align input sequences with template sequences and restore hydroxyprolines.

    Args:
        input_fasta_path (Path): Path to the input FASTA file.
        template_fasta_path (Path): Path to the template FASTA file.
        output_prefix (str): Prefix for output files.
        template_pdb (Path): Path to the template PDB file.

    Returns:
        Tuple[List[SeqRecord], Path, Path]: Aligned and restored sequences,
        MSA output path, and MODELLER output path.

    Raises:
        Exception: If any error occurs during the alignment process.
    """
    try:
        LOG.debug(f"Starting alignment process with input: {input_fasta_path}")
        alignment = Alignment(
            input_fasta_path, template_fasta_path, output_prefix, template_pdb
        )
        msa_output, modeller_output = alignment.align_sequences()
        LOG.debug("Alignment process completed successfully")
        return msa_output, modeller_output
    except Exception as e:
        LOG.error(f"Error in align_sequences: {str(e)}")
        raise
