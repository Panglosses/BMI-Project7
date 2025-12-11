"""
Validation Utilities
"""

from pathlib import Path

from core.data_models import AnalysisConfig


def validate_pdb_file(filepath: str | Path) -> bool:
    """
    Validates a PDB file.

    Args:
        filepath: Path to the PDB file.

    Returns:
        bool: True if the file is valid.
    """
    path = Path(filepath)

    # Check file existence
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {filepath}")

    # Check file size
    if path.stat().st_size == 0:
        raise ValueError(f"PDB file is empty: {filepath}")

    # Check file extension
    if path.suffix.lower() not in [".pdb", ".ent"]:
        raise ValueError(f"File extension is not a recognized PDB format: {filepath}")

    # Basic content validation
    try:
        with open(path, "r") as f:
            first_line = f.readline().strip()
            # PDB files typically start with specific record types
            if not any(
                first_line.startswith(prefix)
                for prefix in ["HEADER", "ATOM", "HETATM", "MODEL"]
            ):
                # Some PDB files might lack a HEADER
                if len(first_line) > 0 and first_line[0] not in [" ", "\t"]:
                    # Check for a valid PDB record type
                    record_type = first_line[0:6].strip()
                    if record_type not in [
                        "HEADER",
                        "ATOM",
                        "HETATM",
                        "MODEL",
                        "ENDMDL",
                        "TER",
                        "END",
                    ]:
                        raise ValueError(
                            f"Invalid PDB file format. First line: {first_line[:50]}..."
                        )
    except UnicodeDecodeError:
        raise ValueError(f"PDB file is not a text file: {filepath}")

    return True


def validate_config(config: AnalysisConfig) -> bool:
    """
    Validates the analysis configuration.

    Args:
        config: Analysis configuration object.

    Returns:
        bool: True if the configuration is valid.
    """
    config.validate()  # Uses the built-in validate method
    return True


def validate_output_dir(directory: str | Path) -> Path:
    """
    Validates and prepares an output directory.

    Args:
        directory: Path to the output directory.

    Returns:
        Path: A validated, writable Path object for the directory.
    """
    path = Path(directory)

    # Create directory if it doesn't exist
    if not path.exists():
        try:
            path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise ValueError(f"Failed to create output directory {directory}: {e}")

    # Ensure it's a directory
    if not path.is_dir():
        raise ValueError(f"Output path is not a directory: {directory}")

    # Simple write-permission test
    test_file = path / ".write_test"
    try:
        test_file.touch()
        test_file.unlink()
    except Exception as e:
        raise ValueError(f"Output directory is not writable {directory}: {e}")

    return path