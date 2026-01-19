#!/usr/bin/env python3
"""Black on save hook - Auto-format Python files with Black after edits."""

import json
import subprocess
import sys
from pathlib import Path


def main():
    """Main hook entry point."""
    # Read input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract file path from tool input
    tool_input = input_data.get("tool_input", {})
    file_path = tool_input.get("file_path", "")

    if not file_path:
        sys.exit(0)

    # Check if it's a Python file
    path = Path(file_path)
    if path.suffix.lower() != ".py":
        sys.exit(0)

    # Check if file exists
    if not path.exists():
        sys.exit(0)

    # Run isort first (if available)
    try:
        subprocess.run(
            ["isort", str(path)],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Run black
    try:
        result = subprocess.run(
            ["black", str(path)],
            capture_output=True,
            text=True,
            timeout=25,
        )

        if result.returncode == 0:
            print(f"Formatted: {file_path}")
        else:
            # Black returns non-zero for "file already formatted" too
            if "unchanged" not in result.stderr.lower():
                print(f"Black warning: {result.stderr}", file=sys.stderr)

    except FileNotFoundError:
        # Black not installed, fail silently
        pass
    except subprocess.TimeoutExpired:
        print("Black timeout", file=sys.stderr)
    except Exception as e:
        print(f"Black error: {e}", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
