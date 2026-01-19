#!/usr/bin/env python3
"""Ralph safety check hook - Block dangerous bash commands during Ralph mode."""

import json
import os
import re
import sys

# Dangerous patterns to block
DANGEROUS_PATTERNS = [
    # Remote code execution
    (r"curl\s+.*\|\s*(ba)?sh", "Remote code execution via curl"),
    (r"wget\s+.*\|\s*(ba)?sh", "Remote code execution via wget"),
    (r"curl\s+-[a-zA-Z]*o\s*-.*\|\s*(ba)?sh", "Remote code execution via curl"),

    # Credential access patterns
    (r"(cat|less|more|head|tail|grep)\s+.*\.ssh/", "Accessing SSH credentials"),
    (r"(cat|less|more|head|tail|grep)\s+.*\.aws/", "Accessing AWS credentials"),
    (r"(cat|less|more|head|tail|grep)\s+.*\.config/gcloud", "Accessing GCloud credentials"),
    (r"(cat|less|more|head|tail|grep)\s+.*/etc/shadow", "Accessing shadow file"),

    # Destructive operations
    (r"rm\s+-[a-zA-Z]*r[a-zA-Z]*f[a-zA-Z]*\s+/\s*$", "Recursive delete of root"),
    (r"rm\s+-[a-zA-Z]*r[a-zA-Z]*f[a-zA-Z]*\s+/\*", "Recursive delete of root contents"),
    (r"rm\s+-[a-zA-Z]*r[a-zA-Z]*f[a-zA-Z]*\s+~/?$", "Recursive delete of home directory"),
    (r"rm\s+-[a-zA-Z]*r[a-zA-Z]*f[a-zA-Z]*\s+~/\*", "Recursive delete of home contents"),
    (r"mkfs\.", "Filesystem creation (destructive)"),
    (r"dd\s+.*of=/dev/", "Direct disk write"),

    # Network reconnaissance
    (r"\bnmap\b", "Network scanning"),
    (r"\bnc\s+-[a-zA-Z]*l", "Netcat listener"),
    (r"\bnetcat\s+-[a-zA-Z]*l", "Netcat listener"),

    # System modification
    (r"chmod\s+777", "Overly permissive chmod"),
    (r"chown\s+root", "Changing ownership to root"),
    (r">\s*/etc/", "Writing to /etc"),
    (r"tee\s+/etc/", "Writing to /etc via tee"),

    # Exfiltration patterns
    (r"curl\s+.*-d\s+.*@", "Data exfiltration via curl POST"),
    (r"curl\s+.*--data.*@", "Data exfiltration via curl POST"),

    # Privilege escalation
    (r"\bsudo\b", "Privilege escalation attempt"),
    (r"\bsu\s+-", "Switching to root user"),
    (r"\bsu\s+root", "Switching to root user"),

    # Persistence mechanisms
    (r"crontab\s+-[a-zA-Z]*e", "Modifying crontab"),
    (r">>\s*~/\.bashrc", "Modifying .bashrc"),
    (r">>\s*~/\.profile", "Modifying .profile"),
    (r">>\s*~/\.zshrc", "Modifying .zshrc"),
]


def is_dangerous(command: str) -> tuple[bool, str]:
    """Check if a command matches dangerous patterns."""
    # Get allowed patterns from environment
    allowed = os.environ.get("RALPH_SAFETY_ALLOW", "")
    allowed_list = [p.strip().lower() for p in allowed.split(",") if p.strip()]

    # Get additional blocked patterns from environment
    extra_blocked = os.environ.get("RALPH_SAFETY_BLOCK", "")
    extra_patterns = [(p.strip(), f"Custom blocked pattern: {p.strip()}")
                      for p in extra_blocked.split(",") if p.strip()]

    all_patterns = DANGEROUS_PATTERNS + extra_patterns

    # Normalize command for checking
    cmd_lower = command.lower()

    for pattern, description in all_patterns:
        # Check if this pattern is allowed
        pattern_name = description.split()[0].lower()
        if pattern_name in allowed_list:
            continue

        if re.search(pattern, command, re.IGNORECASE):
            return True, description

    return False, ""


def main():
    """Main hook entry point."""
    # Check if we're in Ralph mode (optional - can be used standalone too)
    # ralph_mode = os.environ.get("RALPH_MODE", "false").lower() == "true"

    # Read input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract command from tool input
    tool_input = input_data.get("tool_input", {})
    command = tool_input.get("command", "")

    if not command:
        sys.exit(0)

    # Check if command is dangerous
    dangerous, reason = is_dangerous(command)

    if dangerous:
        error_msg = f"""BLOCKED: Dangerous command detected.

Command: {command}
Reason: {reason}

This command was blocked by the ralph-safety-check hook as a safety measure.

If this is a legitimate operation:
1. Ask the user to run it manually, OR
2. Add to RALPH_SAFETY_ALLOW environment variable

Example: RALPH_SAFETY_ALLOW="{reason.split()[0].lower()}" ./ralph.sh"""

        print(error_msg, file=sys.stderr)
        sys.exit(2)  # Exit 2 = block operation

    sys.exit(0)


if __name__ == "__main__":
    main()
