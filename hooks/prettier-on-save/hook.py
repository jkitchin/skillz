#!/usr/bin/env python3
"""Prettier on save hook - Auto-format files with Prettier after edits."""

import json
import subprocess
import sys
from pathlib import Path

# File extensions that Prettier can format
PRETTIER_EXTENSIONS = {
    ".js",
    ".jsx",
    ".ts",
    ".tsx",
    ".css",
    ".scss",
    ".less",
    ".json",
    ".md",
    ".mdx",
    ".html",
    ".htm",
    ".yaml",
    ".yml",
    ".graphql",
    ".gql",
    ".vue",
    ".svelte",
}


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

    # Check if file extension is supported
    path = Path(file_path)
    if path.suffix.lower() not in PRETTIER_EXTENSIONS:
        sys.exit(0)

    # Check if file exists
    if not path.exists():
        sys.exit(0)

    # Run prettier
    try:
        result = subprocess.run(
            ["npx", "prettier", "--write", str(path)],
            capture_output=True,
            text=True,
            timeout=25,
        )

        if result.returncode != 0:
            # Try with global prettier
            result = subprocess.run(
                ["prettier", "--write", str(path)],
                capture_output=True,
                text=True,
                timeout=25,
            )

        if result.returncode == 0:
            print(f"Formatted: {file_path}")
        else:
            print(f"Prettier warning: {result.stderr}", file=sys.stderr)

    except FileNotFoundError:
        # Prettier not installed, fail silently
        pass
    except subprocess.TimeoutExpired:
        print("Prettier timeout", file=sys.stderr)
    except Exception as e:
        print(f"Prettier error: {e}", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
