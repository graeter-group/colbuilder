"""
Colbuilder Error Codes Module

This module defines all error codes and their associated information for the
Colbuilder system. Each error type includes:
- A unique code
- A descriptive message
- A list of suggestions for resolution
- Optional documentation URL

Error categories:
- System (SYS): System-level errors
- Configuration (CFG): Configuration and setup errors
- Sequence (SEQ): Sequence generation errors
- Geometry (GEO): Geometry generation errors
- Topology (TOP): Topology generation errors
"""

from typing import Dict, List, NamedTuple, Optional

class ErrorInfo(NamedTuple):
    """
    Information associated with an error code.
    
    Attributes:
        code: Unique error identifier
        message: Human-readable error description
        suggestions: List of potential solutions
        docs_url: Optional link to documentation
    """
    code: str
    message: str
    suggestions: List[str]
    docs_url: Optional[str] = None

# System-related errors
SYSTEM_ERRORS: Dict[str, ErrorInfo] = {
    "SYS_ERR_001": ErrorInfo(
        code="SYS_ERR_001",
        message="Unexpected system error occurred during operation",
        suggestions=[
            "Check system resources (CPU, memory, disk space)",
            "Verify all required dependencies are installed",
            "Check logs for detailed error information",
            "Ensure sufficient permissions for file operations"
        ]
    ),
    "SYS_ERR_002": ErrorInfo(
        code="SYS_ERR_002",
        message="Failed to import required module",
        suggestions=[
            "Verify all required packages are installed",
            "Check Python environment configuration",
            "Ensure module path is correct",
            "Check for version conflicts between dependencies"
        ]
    )
}

# Configuration-related errors
CONFIGURATION_ERRORS: Dict[str, ErrorInfo] = {
    "CFG_ERR_001": ErrorInfo(
        code="CFG_ERR_001",
        message="Configuration initialization failed",
        suggestions=[
            "Check if all required fields are provided",
            "Verify the YAML syntax in your configuration file",
            "Ensure all paths exist and are accessible",
            "Validate configuration values against requirements"
        ],
        docs_url="https://colbuilder.readthedocs.io/en/latest/configuration.html"
    ),
    "CFG_ERR_002": ErrorInfo(
        code="CFG_ERR_002",
        message="Required file or directory not found",
        suggestions=[
            "Check if the path exists",
            "Verify file permissions",
            "Ensure the directory structure is correct",
            "Check for case sensitivity in filenames"
        ]
    ),
    "CFG_ERR_003": ErrorInfo(
        code="CFG_ERR_003",
        message="Invalid species configuration",
        suggestions=[
            "Use a supported species or provide custom FASTA file",
            "Check species name spelling",
            "Verify species is in supported list",
            "Ensure species data files are properly installed"
        ],
        docs_url="https://colbuilder.readthedocs.io/en/latest/configuration.html#supported-species"
    ),
    "CFG_ERR_004": ErrorInfo(
        code="CFG_ERR_004",
        message="Invalid mixing configuration",
        suggestions=[
            "Ensure mixing ratios sum to exactly 100%",
            "Use format 'Type:percentage Type:percentage'",
            "Verify all required files are provided",
            "Check that mixing types match available options"
        ]
    ),
    "CFG_ERR_005": ErrorInfo(
        code="CFG_ERR_005",
        message="Invalid file format detected",
        suggestions=[
            "Check file format requirements in documentation",
            "Verify file is not corrupted",
            "Ensure file meets structural requirements",
            "Validate file content against schema"
        ]
    ),
    "CFG_ERR_006": ErrorInfo(
        code="CFG_ERR_006",
        message="Invalid parameter value specified",
        suggestions=[
            "Contact distance must be a positive number",
            "Fibril length must be a positive number less than 334",
            "Terminal combinations must follow format 'ResidueNumber.Chain - ResidueNumber.Chain'",
            "Force field must be either 'amber99' or 'martini'"
        ]
    )
}

# Sequence-related errors
SEQUENCE_ERRORS: Dict[str, ErrorInfo] = {
    "SEQ_ERR_001": ErrorInfo(
        code="SEQ_ERR_001",
        message="Invalid sequence format detected",
        suggestions=[
            "Ensure sequence is in valid FASTA format",
            "Check for invalid characters in sequence",
            "Verify sequence matches collagen requirements",
            "Validate sequence length and composition"
        ],
        docs_url="https://colbuilder.readthedocs.io/en/latest/input_formats.html#sequence-requirements"
    ),
    "SEQ_ERR_002": ErrorInfo(
        code="SEQ_ERR_002",
        message="Invalid species or crosslink specification",
        suggestions=[
            "Verify species name is correct",
            "Check if specified crosslink type is available for this species",
            "Ensure crosslink combination is valid",
            "Validate against supported crosslink types"
        ],
        docs_url="https://colbuilder.readthedocs.io/en/latest/crosslinks.html#available-crosslinks"
    ),
    "SEQ_ERR_003": ErrorInfo(
        code="SEQ_ERR_003",
        message="Structure optimization failed during sequence generation",
        suggestions=[
            "Try running sequence generation again",
            "Check if crosslink specifications are valid",
            "Consider using different crosslink combinations",
            "Verify input structure meets requirements"
        ]
    ),
    "SEQ_ERR_004": ErrorInfo(
        code="SEQ_ERR_004",
        message="Required input file not found for sequence generation",
        suggestions=[
            "Ensure all required files are uploaded",
            "Check if FASTA file was uploaded correctly",
            "Verify file paths in configuration",
            "Check file permissions"
        ]
    )
}

# Geometry-related errors
GEOMETRY_ERRORS: Dict[str, ErrorInfo] = {
    "GEO_ERR_001": ErrorInfo(
        code="GEO_ERR_001",
        message="Invalid input parameters for geometry generation",
        suggestions=[
            "Provide either contact distance or crystal contacts file",
            "Check if PDB file is correctly formatted",
            "Verify input parameters match requirements",
            "Ensure all required fields are provided"
        ],
        docs_url="https://colbuilder.readthedocs.io/en/latest/geometry.html#input-requirements"
    ),
    "GEO_ERR_002": ErrorInfo(
        code="GEO_ERR_002",
        message="Crystal contacts generation failed",
        suggestions=[
            "Ensure CRYST1 record is present at PDB file top",
            "CRYST1 record must contain unit cell parameters (a, b, c, alpha, beta, gamma)",
            "Check if PDB file follows standard format",
            "Verify crystal symmetry information is correct"
        ]
    ),
    "GEO_ERR_003": ErrorInfo(
        code="GEO_ERR_003",
        message="Geometry mixing operation failed",
        suggestions=[
            "Verify all required PDB files are provided",
            "Check if mixing ratios sum to 100%",
            "Ensure crosslink types match available options",
            "Validate mixing parameters"
        ]
    ),
    "GEO_ERR_004": ErrorInfo(
        code="GEO_ERR_004",
        message="Crosslink replacement operation failed",
        suggestions=[
            "Check if replacement ratio is between 0 and 100",
            "Verify structure contains crosslinks to replace",
            "Ensure replacement types are valid",
            "Check replacement file format"
        ]
    ),
    "GEO_ERR_005": ErrorInfo(
        code="GEO_ERR_005",
        message="Invalid ATOM record format in PDB file",
        suggestions=[
            "Check if PDB file follows standard format",
            "Ensure ATOM records contain all required fields",
            "Verify chain identifiers are in column 22",
            "Use PDB validation tool to check format"
        ]
    ),
    "GEO_ERR_006": ErrorInfo(
        code="GEO_ERR_006",
        message="Invalid number of chains in structure",
        suggestions=[
            "Ensure PDB file contains exactly 3 chains",
            "Check chain ID assignments in ATOM records",
            "Verify each chain represents complete collagen strand",
            "Ensure chain identifiers are consistent"
        ]
    ),
    "GEO_ERR_007": ErrorInfo(
        code="GEO_ERR_007",
        message="Invalid termination records in structure",
        suggestions=[
            "Add TER records after each chain",
            "Ensure each chain is properly terminated",
            "Check TER record format",
            "Verify chain terminations match chain count"
        ]
    )
}

# Topology-related errors
TOPOLOGY_ERRORS: Dict[str, ErrorInfo] = {
    "TOP_ERR_001": ErrorInfo(
        code="TOP_ERR_001",
        message="Invalid force field specification",
        suggestions=[
            "Use a supported force field (currently only amber99)",
            "Check force field spelling",
            "Verify force field is properly configured",
            "Ensure force field files are available in the expected location"
        ]
    ),
    "TOP_ERR_002": ErrorInfo(
        code="TOP_ERR_002",
        message="Force field directory setup failed",
        suggestions=[
            "Verify force field directory exists",
            "Check file permissions",
            "Ensure sufficient disk space",
            "Validate force field file structure"
        ]
    ),
    "TOP_ERR_003": ErrorInfo(
        code="TOP_ERR_003",
        message="Missing required force field files",
        suggestions=[
            "Check if all required files are present in force field directory",
            "Verify file permissions",
            "Ensure force field installation is complete",
            "Check for corrupted force field files"
        ]
    ),
    "TOP_ERR_004": ErrorInfo(
        code="TOP_ERR_004",
        message="Invalid or missing merged PDB file",
        suggestions=[
            "Verify PDB file was generated correctly",
            "Check file permissions and path",
            "Ensure sufficient disk space",
            "Validate PDB file format"
        ]
    ),
    "TOP_ERR_005": ErrorInfo(
        code="TOP_ERR_005",
        message="GROMACS pdb2gmx execution failed",
        suggestions=[
            "Check GROMACS installation",
            "Verify PDB file format",
            "Ensure all required force field parameters are available",
            "Check GROMACS error message for specific issues"
        ]
    ),
    "TOP_ERR_006": ErrorInfo(
        code="TOP_ERR_006",
        message="No models successfully processed",
        suggestions=[
            "Check input models for validity",
            "Verify force field compatibility",
            "Review error messages for specific model failures",
            "Ensure proper model format and structure"
        ]
    ),
    "TOP_ERR_007": ErrorInfo(
        code="TOP_ERR_007",
        message="Failed to write final topology files",
        suggestions=[
            "Check write permissions in output directory",
            "Ensure sufficient disk space",
            "Verify topology generation completed successfully",
            "Check for conflicting file names"
        ]
    )
}