#!/usr/bin/env python3
"""Protect secrets hook - Block writes to sensitive files."""

from __future__ import annotations

import fnmatch
import json
import os
import sys
from pathlib import Path

# Default protected patterns
PROTECTED_PATTERNS = [
    # Environment files
    ".env",
    ".env.*",
    "*.env",
    # Credential files
    "credentials.json",
    "credentials.yaml",
    "credentials.yml",
    "secrets.json",
    "secrets.yaml",
    "secrets.yml",
    "secret.json",
    "secret.yaml",
    "secret.yml",
    # Private keys
    "*.pem",
    "*.key",
    "*.p12",
    "*.pfx",
    "id_rsa",
    "id_rsa.*",
    "id_ed25519",
    "id_ed25519.*",
    "id_ecdsa",
    "id_ecdsa.*",
    "id_dsa",
    "id_dsa.*",
    # AWS
    ".aws/credentials",
    ".aws/config",
    # SSH
    ".ssh/config",
    ".ssh/known_hosts",
    ".ssh/authorized_keys",
    # Git credentials
    ".git-credentials",
    ".gitconfig",
    ".netrc",
    # Docker
    ".docker/config.json",
    # Kubernetes
    "kubeconfig",
    ".kube/config",
    # NPM
    ".npmrc",
    # PyPI
    ".pypirc",
    # GCP
    "gcloud/credentials.db",
    "application_default_credentials.json",
    # Azure
    ".azure/",
    # Terraform
    "*.tfvars",
    "terraform.tfstate",
    "terraform.tfstate.backup",
]


def is_protected(file_path: str) -> tuple[bool, str]:
    """Check if a file path matches protected patterns."""
    path = Path(file_path)
    name = path.name

    # Get extra patterns from environment
    extra_patterns = os.environ.get("PROTECT_SECRETS_EXTRA", "")
    extra_list = [p.strip() for p in extra_patterns.split(",") if p.strip()]

    # Get allowed patterns from environment
    allowed = os.environ.get("PROTECT_SECRETS_ALLOW", "")
    allowed_list = [p.strip() for p in allowed.split(",") if p.strip()]

    # Check allowed list first
    for pattern in allowed_list:
        if fnmatch.fnmatch(name, pattern) or fnmatch.fnmatch(str(path), pattern):
            return False, ""

    # Check against protected patterns
    all_patterns = PROTECTED_PATTERNS + extra_list

    for pattern in all_patterns:
        # Match against filename
        if fnmatch.fnmatch(name, pattern):
            return True, pattern
        # Match against full path
        if fnmatch.fnmatch(str(path), f"*{pattern}") or fnmatch.fnmatch(str(path), f"*/{pattern}"):
            return True, pattern
        # Match against path parts
        for part in path.parts:
            if fnmatch.fnmatch(part, pattern):
                return True, pattern

    return False, ""


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

    # Check if file is protected
    protected, pattern = is_protected(file_path)

    if protected:
        error_msg = f"""BLOCKED: Cannot write to protected file.

File: {file_path}
Matched pattern: {pattern}

This file appears to contain secrets or sensitive configuration.
To proceed, please:
1. Ask the user to manually edit this file, OR
2. Use a template file (e.g., .env.example) instead

Set PROTECT_SECRETS_ALLOW="{Path(file_path).name}" to allow this file."""

        print(error_msg, file=sys.stderr)
        sys.exit(2)  # Exit 2 = block operation

    sys.exit(0)


if __name__ == "__main__":
    main()
