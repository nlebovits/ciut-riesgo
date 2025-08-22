#!/usr/bin/env python3
# numpydoc ignore=GL08
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Match


def process_python_code(
    code: str, ruff_args: List[str], verbose: bool = False
) -> str:  # numpydoc ignore=RT01
    """Process Python code using Ruff with custom arguments."""
    # Use temporary file approach like the Lua extension for better reliability
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False, encoding="utf-8") as temp_file:
        temp_file.write(code)
        temp_filepath = temp_file.name
    
    try:
        # Process the temporary file with the specified ruff command
        cmd = ["ruff"] + ruff_args + [temp_filepath]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        
        # Read the formatted content back
        with open(temp_filepath, 'r', encoding='utf-8') as f:
            formatted_content = f.read()
        
        return formatted_content
        
    except subprocess.CalledProcessError as e:
        if verbose:
            print(f"Ruff command: {' '.join(cmd)}", file=sys.stderr)
            print(f"Code block content preview: {code[:100]}...", file=sys.stderr)
            print(f"Temporary file: {temp_filepath}", file=sys.stderr)
        
        print(
            f"Error: Failed to process Python code with Ruff. "
            f"Exit code: {e.returncode}", 
            file=sys.stderr
        )
        if e.stderr:
            print(f"Ruff stderr: {e.stderr}", file=sys.stderr)
        return code
        
    finally:
        # Clean up temporary file
        try:
            os.unlink(temp_filepath)
        except OSError:
            pass  # File already deleted or doesn't exist


def replace_code_block(
    match: Match[str], ruff_args: List[str], verbose: bool = False
) -> str:  # numpydoc ignore=RT01
    """Replace code block with processed version."""
    return f"{match.group(1)}\n{process_python_code(match.group(2), ruff_args, verbose)}{match.group(3)}"


def process_file(
    filepath: Path, ruff_args: List[str], verbose: bool = False
) -> None:  # numpydoc ignore=RT01
    """Process the given file, formatting Python code blocks."""
    # Check for unsupported operations
    if "check" in ruff_args or "--fix" in ruff_args:
        print("Error: 'ruff check' and '--fix' operations are not supported.", file=sys.stderr)
        print("This tool only supports 'ruff format' operations on individual code blocks.", file=sys.stderr)
        print("Use 'ruff check' directly on extracted Python files for linting.", file=sys.stderr)
        sys.exit(1)
    
    python_code_block_pattern = r"(```\{python\})(.*?)(```)"
    try:
        content = filepath.read_text()
        formatted_content = re.sub(
            python_code_block_pattern,
            lambda m: replace_code_block(m, ruff_args, verbose),
            content,
            flags=re.DOTALL,
        )

        with tempfile.NamedTemporaryFile(mode="w", delete=False, encoding="utf-8") as temp_file:
            temp_file.write(formatted_content)
            temp_filepath = Path(temp_file.name)

        shutil.move(str(temp_filepath), str(filepath))
        
        if verbose:
            print(f"Successfully processed file: {filepath}", file=sys.stderr)
            
    except IOError as e:
        print(f"Error processing file {filepath}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            'Usage: python tools/qmd_ruff.py "RUFF_FORMAT_ARGS" <filename1.qmd> [filename2.qmd ...] [--verbose]'
        )
        print('Examples:')
        print('       python tools/qmd_ruff.py "format" renabap.qmd')
        print('       python tools/qmd_ruff.py "format --diff" renabap.qmd --verbose')
        print('       python tools/qmd_ruff.py "" renabap.qmd  # Basic formatting')
        print('')
        print('Note: Only ruff format operations are supported. Use ruff check directly on extracted Python files for linting.')
        sys.exit(1)

    # Check for verbose flag
    verbose = "--verbose" in sys.argv
    if verbose:
        sys.argv.remove("--verbose")

    ruff_args = sys.argv[1].split()

    missing_files = [file for file in sys.argv[2:] if not Path(file).exists()]
    if missing_files:
        raise FileNotFoundError(
            f"The following file(s) do not exist: {', '.join(missing_files)}."
        )
    
    for filepath in sys.argv[2:]:
        path = Path(filepath)
        if verbose:
            print(f"Processing file: {filepath}", file=sys.stderr)
        process_file(path, ruff_args, verbose)