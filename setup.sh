#!/bin/bash
# Setup script for Bioinformatics Project - Solvent Accessibility Analysis Toolkit

set -e  # Exit on error

echo "========================================"
echo "Solvent Accessibility Analysis Toolkit"
echo "Setup Script"
echo "========================================"
echo ""

# Check Python version
REQUIRED_PYTHON="3.14"
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' || python --version 2>&1 | awk '{print $2}')

if [[ -z "$PYTHON_VERSION" ]]; then
    echo "ERROR: Python is not installed or not in PATH."
    echo "Please install Python $REQUIRED_PYTHON or higher."
    exit 1
fi

echo "Detected Python version: $PYTHON_VERSION"

# Compare versions
if [[ "$(printf '%s\n' "$REQUIRED_PYTHON" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_PYTHON" ]]; then
    echo "ERROR: Python $REQUIRED_PYTHON or higher is required."
    echo "Current version: $PYTHON_VERSION"
    exit 1
fi

echo "✓ Python version meets requirement."
echo ""

# Check if virtual environment already exists
VENV_DIR="./bmi"
if [[ -d "$VENV_DIR" ]]; then
    echo "Virtual environment already exists at: $VENV_DIR"
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
        echo "Creating new virtual environment..."
        python -m venv "$VENV_DIR"
        echo "✓ Virtual environment created."
    else
        echo "Using existing virtual environment."
    fi
else
    echo "Creating virtual environment at: $VENV_DIR"
    python -m venv "$VENV_DIR"
    echo "✓ Virtual environment created."
fi

echo ""

# Activate virtual environment based on OS
if [[ "$OSTYPE" == "linux-gnu"* ]] || [[ "$OSTYPE" == "darwin"* ]]; then
    # Linux or macOS
    ACTIVATE_SCRIPT="$VENV_DIR/bin/activate"
    if [[ -f "$ACTIVATE_SCRIPT" ]]; then
        source "$ACTIVATE_SCRIPT"
        echo "Virtual environment activated (Linux/macOS)."
    else
        echo "WARNING: Could not find activation script: $ACTIVATE_SCRIPT"
        echo "Please activate manually: source $ACTIVATE_SCRIPT"
    fi
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    # Windows (Git Bash, Cygwin)
    ACTIVATE_SCRIPT="$VENV_DIR/Scripts/activate"
    if [[ -f "$ACTIVATE_SCRIPT" ]]; then
        source "$ACTIVATE_SCRIPT"
        echo "Virtual environment activated (Windows Git Bash)."
    else
        echo "WARNING: Could not find activation script: $ACTIVATE_SCRIPT"
        echo "Please activate manually: source $ACTIVATE_SCRIPT"
    fi
else
    echo "NOTE: Unsupported OS type: $OSTYPE"
    echo "Please activate virtual environment manually:"
    echo "  Linux/macOS: source $VENV_DIR/bin/activate"
    echo "  Windows:     $VENV_DIR\\Scripts\\Activate.ps1 (PowerShell)"
    echo "  Windows:     $VENV_DIR\\Scripts\\activate.bat (CMD)"
fi

echo ""

# Upgrade pip
echo "Upgrading pip..."
python -m pip install --upgrade pip
echo "✓ pip upgraded."

echo ""

# Install dependencies from requirements.txt
REQUIREMENTS_FILE="requirements.txt"
if [[ -f "$REQUIREMENTS_FILE" ]]; then
    echo "Installing dependencies from $REQUIREMENTS_FILE..."
    pip install -r "$REQUIREMENTS_FILE"
    echo "✓ Dependencies installed."
else
    echo "ERROR: $REQUIREMENTS_FILE not found in current directory."
    echo "Please ensure you are in the project root directory."
    exit 1
fi

echo ""
echo "========================================"
echo "Setup completed successfully!"
echo "========================================"
echo ""
echo "To use the project:"
echo ""
echo "1. Activate virtual environment:"
if [[ "$OSTYPE" == "linux-gnu"* ]] || [[ "$OSTYPE" == "darwin"* ]]; then
    echo "   source $VENV_DIR/bin/activate"
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "   source $VENV_DIR/Scripts/activate"
else
    echo "   Linux/macOS: source $VENV_DIR/bin/activate"
    echo "   Windows PowerShell: .\\$VENV_DIR\\Scripts\\Activate.ps1"
    echo "   Windows CMD: $VENV_DIR\\Scripts\\activate.bat"
fi
echo ""
echo "2. Run the analysis:"
echo "   python -m solvent_analysis --help"
echo ""
echo "3. Deactivate virtual environment when done:"
echo "   deactivate"
echo ""
echo "For detailed usage, see README.md"
echo "========================================"