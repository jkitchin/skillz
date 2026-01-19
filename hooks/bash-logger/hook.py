#!/usr/bin/env python3
"""Bash logger hook - Log all bash commands for audit."""

import json
import os
import sys
from datetime import datetime
from pathlib import Path


def get_log_file() -> Path:
    """Get the log file path from environment or default."""
    log_file = os.environ.get("BASH_LOGGER_FILE", "~/.claude/logs/bash-commands.log")
    return Path(os.path.expanduser(log_file))


def rotate_log_if_needed(log_file: Path, max_size_mb: int = 10):
    """Rotate log file if it exceeds max size."""
    if not log_file.exists():
        return

    size_mb = log_file.stat().st_size / (1024 * 1024)
    if size_mb < max_size_mb:
        return

    # Rotate logs (keep 5 backups)
    for i in range(4, 0, -1):
        old = log_file.with_suffix(f".log.{i}")
        new = log_file.with_suffix(f".log.{i + 1}")
        if old.exists():
            old.rename(new)

    # Rotate current log
    log_file.rename(log_file.with_suffix(".log.1"))


def format_entry(
    timestamp: datetime,
    session_id: str,
    cwd: str,
    command: str,
    description: str,
    format_type: str = "full",
) -> str:
    """Format a log entry."""
    ts = timestamp.strftime("%Y-%m-%d %H:%M:%S")

    if format_type == "json":
        entry = {
            "timestamp": timestamp.isoformat(),
            "session_id": session_id,
            "cwd": cwd,
            "command": command,
            "description": description,
        }
        return json.dumps(entry)

    elif format_type == "compact":
        return f"[{ts}] [{session_id[:8]}] {command}"

    else:  # full
        lines = [
            f"[{ts}] [session:{session_id[:12]}] [cwd:{cwd}]",
            f"Command: {command}",
        ]
        if description:
            lines.append(f"Description: {description}")
        lines.append("---")
        return "\n".join(lines)


def should_exclude(command: str) -> bool:
    """Check if command should be excluded from logging."""
    import re

    exclude_pattern = os.environ.get("BASH_LOGGER_EXCLUDE", "")
    if not exclude_pattern:
        return False

    try:
        return bool(re.match(exclude_pattern, command))
    except re.error:
        return False


def main():
    """Main hook entry point."""
    # Read input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract relevant fields
    session_id = input_data.get("session_id", "unknown")
    cwd = input_data.get("cwd", "")
    tool_input = input_data.get("tool_input", {})

    command = tool_input.get("command", "")
    description = tool_input.get("description", "")

    if not command:
        sys.exit(0)

    # Check exclusion pattern
    if should_exclude(command):
        sys.exit(0)

    # Get configuration
    log_file = get_log_file()
    format_type = os.environ.get("BASH_LOGGER_FORMAT", "full").lower()

    # Ensure log directory exists
    log_file.parent.mkdir(parents=True, exist_ok=True)

    # Rotate if needed
    rotate_log_if_needed(log_file)

    # Format entry
    timestamp = datetime.now()
    entry = format_entry(timestamp, session_id, cwd, command, description, format_type)

    # Append to log
    try:
        with open(log_file, "a") as f:
            f.write(entry + "\n")
    except Exception as e:
        print(f"Error writing to log: {e}", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
