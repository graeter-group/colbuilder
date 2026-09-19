"""
This module provides validation utilities for biological file formats, including FASTA and PDB files,
to ensure compatibility with the ColBuilder pipeline. It includes functionality for validating input
files, parsing content, and checking structural integrity.

Key Features:
--------------
1. **FASTA Validation**:
   - Validate the existence and content of FASTA files.
   - Ensure sequences are properly formatted with headers and contain exactly three sequences.
   - Check sequence lengths against configurable minimum and maximum thresholds.

2. **PDB Validation**:
   - Validate the existence and content of PDB files.
   - Ensure the presence of required records (e.g., CRYST1, ATOM, TER).
   - Verify the correct number of chains (exactly three) and proper termination with TER records.

3. **Error Handling**:
   - Raise detailed exceptions for invalid files or structural issues.
   - Provide warnings for non-critical issues, such as sequence length deviations.

Usage:
------
This module is designed to be used as part of the ColBuilder pipeline to validate input files before
processing. The `BioformatValidator` class provides methods for validating FASTA and PDB files.

Example:
--------
```python
from pathlib import Path
from colbuilder.core.utils.validators import BioformatValidator

# Initialize the validator
validator = BioformatValidator(min_sequence_length=950, max_sequence_length=1100)

# Validate input files
fasta_path = Path("input_sequences.fasta")
pdb_path = Path("input_structure.pdb")

try:
    warnings = validator.validate_input_files(fasta_path=fasta_path, pdb_path=pdb_path)
    if warnings:
        print("Warnings:", warnings)
    print("Validation successful!")
except Exception as e:
    print(f"Validation failed: {e}")
```
"""

from dataclasses import dataclass
from typing import List, Optional
from pathlib import Path

from .exceptions import (
    ColbuilderError,
    SequenceGenerationError,
    GeometryGenerationError,
)


class BioformatValidator:
    """Validates biological file formats including FASTA and PDB files."""

    def __init__(self, min_sequence_length: int = 950, max_sequence_length: int = 1100):
        """
        Initialize the validator with customizable parameters.

        Args:
            min_sequence_length: Minimum acceptable sequence length
            max_sequence_length: Maximum acceptable sequence length
        """
        self.min_sequence_length = min_sequence_length
        self.max_sequence_length = max_sequence_length

    def validate_input_files(
        self, fasta_path: Optional[Path] = None, pdb_path: Optional[Path] = None
    ) -> List[str]:
        """
        Validates input files if they are provided.

        Args:
            fasta_path: Optional path to the FASTA file
            pdb_path: Optional path to the PDB file

        Returns:
            List of warning messages

        Raises:
            SequenceGenerationError: If FASTA validation fails
            GeometryGenerationError: If PDB validation fails
        """
        warnings = []

        if fasta_path is not None:
            if not fasta_path.exists():
                raise SequenceGenerationError(
                    message=f"FASTA file not found: {fasta_path}",
                    error_code="SEQ_ERR_001",
                )

            try:
                content = fasta_path.read_text()
                errors = []

                lines = [line.strip() for line in content.split("\n") if line.strip()]
                if not lines:
                    errors.append("Empty FASTA file")

                sequences = self._parse_fasta_content(lines, errors)
                self._validate_sequences(sequences, errors, warnings)

                if errors:
                    raise SequenceGenerationError(
                        message="FASTA validation failed",
                        error_code="SEQ_ERR_002",
                        original_error=Exception("; ".join(errors)),
                    )

            except Exception as e:
                if not isinstance(e, ColbuilderError):
                    raise SequenceGenerationError(
                        message="Error reading FASTA file",
                        error_code="SEQ_ERR_003",
                        original_error=e,
                    )
                raise

        if pdb_path is not None:
            if not pdb_path.exists():
                raise GeometryGenerationError(
                    message=f"PDB file not found: {pdb_path}", error_code="GEO_ERR_001"
                )

            try:
                content = pdb_path.read_text()
                errors = []

                self._validate_pdb_content(content, errors)

                if errors:
                    raise GeometryGenerationError(
                        message="PDB validation failed",
                        error_code="GEO_ERR_002",
                        original_error=Exception("; ".join(errors)),
                    )

            except Exception as e:
                if not isinstance(e, ColbuilderError):
                    raise GeometryGenerationError(
                        message="Error reading PDB file",
                        error_code="GEO_ERR_003",
                        original_error=e,
                    )
                raise

        return warnings

    def _parse_fasta_content(self, lines: List[str], errors: List[str]) -> dict:
        """Parse FASTA content and collect any errors."""
        sequences = {}
        current_header = None
        current_sequence = []
        sequence_count = 0

        for line in lines:
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_sequence)
                    current_sequence = []

                current_header = line
                sequence_count += 1
            else:
                if not current_header:
                    errors.append("Sequence found before header")
                    continue
                current_sequence.append(line)

        if current_header and current_sequence:
            sequences[current_header] = "".join(current_sequence)

        if sequence_count != 3:
            errors.append(f"Expected exactly 3 sequences, found {sequence_count}")
        elif len(sequences) != sequence_count:
            errors.append("Duplicate FASTA headers found; each sequence must have a unique header")

        return sequences

    def _validate_sequences(
        self, sequences: dict, errors: List[str], warnings: List[str]
    ) -> None:
        """Validate sequence lengths and add warnings if outside acceptable range."""
        for header, sequence in sequences.items():
            length = len(sequence)
            if length < self.min_sequence_length or length > self.max_sequence_length:
                warnings.append(
                    f"Sequence {header} length ({length}) is outside recommended range "
                    f"({self.min_sequence_length}-{self.max_sequence_length})"
                )

    def _validate_pdb_content(self, content: str, errors: List[str]) -> None:
        """Validate PDB content and collect errors."""
        if not content.strip():
            raise GeometryGenerationError(
                message="Empty PDB file provided", error_code="GEO_ERR_001"
            )
        lines = content.split("\n")

        cryst_found = False
        chains = set()
        current_chain = None
        ter_count = 0

        for line in lines:
            if not line.strip():
                continue

            if line.startswith("CRYST1"):
                cryst_found = True
            elif line.startswith("ATOM"):
                try:
                    chain_id = line[21]
                    chains.add(chain_id)
                    current_chain = chain_id
                except IndexError:
                    raise GeometryGenerationError(
                        message=f"Invalid ATOM record format found: {line}",
                        error_code="GEO_ERR_005",
                    )
            elif line.startswith("TER"):
                if current_chain:
                    ter_count += 1

        # Check for CRYST1 record
        if not cryst_found:
            raise GeometryGenerationError(
                message="Missing CRYST1 record in PDB file. Crystal contacts information is required.",
                error_code="GEO_ERR_002",
            )

        # Check for correct number of chains
        if len(chains) != 3:
            raise GeometryGenerationError(
                message=f"Invalid number of chains in PDB file. Expected exactly 3 chains, found {len(chains)}. "
                f"Found chains: {', '.join(sorted(chains)) if chains else 'none'}",
                error_code="GEO_ERR_006",
            )

        # Check for correct number of TER records
        if ter_count != 3:
            raise GeometryGenerationError(
                message=f"Invalid number of TER records in PDB file. Expected 3 TER records, found {ter_count}. "
                "Each chain must be properly terminated with a TER record.",
                error_code="GEO_ERR_007",
            )
