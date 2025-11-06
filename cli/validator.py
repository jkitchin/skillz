"""Validation for skills and commands."""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from cli.utils import validate_description, validate_name


class ValidationError(Exception):
    """Custom exception for validation errors."""

    pass


class SkillValidator:
    """Validator for skill directories and SKILL.md files."""

    REQUIRED_FRONTMATTER_FIELDS = ["name", "description"]
    OPTIONAL_FRONTMATTER_FIELDS = ["allowed-tools"]
    VALID_TOOLS = {
        "Bash",
        "Read",
        "Write",
        "Edit",
        "Glob",
        "Grep",
        "Task",
        "WebFetch",
        "WebSearch",
        "TodoWrite",
        "AskUserQuestion",
        "Skill",
        "SlashCommand",
        "NotebookEdit",
        "BashOutput",
        "KillShell",
    }

    @classmethod
    def validate_skill_directory(cls, skill_path: Path) -> Tuple[bool, List[str]]:
        """
        Validate a skill directory.

        Args:
            skill_path: Path to skill directory

        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        errors = []

        # Check if directory exists
        if not skill_path.exists():
            errors.append(f"Directory does not exist: {skill_path}")
            return False, errors

        if not skill_path.is_dir():
            errors.append(f"Path is not a directory: {skill_path}")
            return False, errors

        # Check for SKILL.md
        skill_file = skill_path / "SKILL.md"
        if not skill_file.exists():
            errors.append(f"Missing SKILL.md in {skill_path}")
            return False, errors

        # Validate SKILL.md content
        valid, file_errors = cls.validate_skill_file(skill_file)
        errors.extend(file_errors)

        return len(errors) == 0, errors

    @classmethod
    def validate_skill_file(cls, skill_file: Path) -> Tuple[bool, List[str]]:
        """
        Validate a SKILL.md file.

        Args:
            skill_file: Path to SKILL.md file

        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        errors = []

        try:
            with open(skill_file, "r") as f:
                content = f.read()
        except Exception as e:
            errors.append(f"Error reading file: {e}")
            return False, errors

        # Parse frontmatter
        frontmatter = cls._parse_frontmatter(content)
        if frontmatter is None:
            errors.append("Missing or invalid YAML frontmatter")
            return False, errors

        # Validate required fields
        for field in cls.REQUIRED_FRONTMATTER_FIELDS:
            if field not in frontmatter:
                errors.append(f"Missing required field: {field}")

        # Validate name
        if "name" in frontmatter:
            name = frontmatter["name"]
            if not validate_name(name):
                errors.append(
                    f"Invalid name '{name}': must be lowercase with hyphens/numbers only, "
                    f"max 64 chars"
                )

        # Validate description
        if "description" in frontmatter:
            desc = frontmatter["description"]
            if not validate_description(desc):
                errors.append(f"Description too long: max 1024 characters")

        # Validate allowed-tools
        if "allowed-tools" in frontmatter:
            tools = frontmatter["allowed-tools"]
            if tools and tools != "*":
                if isinstance(tools, list):
                    for tool in tools:
                        if tool not in cls.VALID_TOOLS:
                            errors.append(f"Unknown tool: {tool}")
                else:
                    errors.append("allowed-tools must be a list or '*'")

        return len(errors) == 0, errors

    @staticmethod
    def _parse_frontmatter(content: str) -> Optional[Dict]:
        """
        Parse YAML frontmatter from markdown content.

        Args:
            content: Markdown content

        Returns:
            Dictionary of frontmatter or None if invalid
        """
        # Match frontmatter between --- delimiters
        pattern = r"^---\s*\n(.*?)\n---\s*\n"
        match = re.match(pattern, content, re.DOTALL)

        if not match:
            return None

        frontmatter_str = match.group(1)
        try:
            return yaml.safe_load(frontmatter_str)
        except yaml.YAMLError:
            return None


class CommandValidator:
    """Validator for command files."""

    OPTIONAL_FRONTMATTER_FIELDS = [
        "description",
        "allowed-tools",
        "model",
        "argument-hint",
        "disable-model-invocation",
    ]

    @classmethod
    def validate_command_file(cls, command_file: Path) -> Tuple[bool, List[str]]:
        """
        Validate a command file.

        Args:
            command_file: Path to command markdown file

        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        errors = []

        if not command_file.exists():
            errors.append(f"File does not exist: {command_file}")
            return False, errors

        if not command_file.suffix == ".md":
            errors.append(f"Command file must be a .md file: {command_file}")
            return False, errors

        try:
            with open(command_file, "r") as f:
                content = f.read()
        except Exception as e:
            errors.append(f"Error reading file: {e}")
            return False, errors

        # Parse optional frontmatter
        frontmatter = cls._parse_frontmatter(content)

        # Validate frontmatter if present
        if frontmatter:
            # Validate description length
            if "description" in frontmatter:
                desc = frontmatter["description"]
                if not validate_description(desc, max_length=256):
                    errors.append("Description too long: max 256 characters for commands")

            # Validate model
            if "model" in frontmatter:
                model = frontmatter["model"]
                valid_models = ["sonnet", "opus", "haiku"]
                if model not in valid_models:
                    errors.append(
                        f"Invalid model '{model}': must be one of {valid_models}"
                    )

            # Validate allowed-tools
            if "allowed-tools" in frontmatter:
                tools = frontmatter["allowed-tools"]
                if tools and tools != "*":
                    if not isinstance(tools, list):
                        errors.append("allowed-tools must be a list or '*'")

        # Check content is not empty (excluding frontmatter)
        content_without_frontmatter = cls._remove_frontmatter(content)
        if not content_without_frontmatter.strip():
            errors.append("Command content is empty")

        return len(errors) == 0, errors

    @staticmethod
    def _parse_frontmatter(content: str) -> Optional[Dict]:
        """Parse YAML frontmatter from markdown content."""
        pattern = r"^---\s*\n(.*?)\n---\s*\n"
        match = re.match(pattern, content, re.DOTALL)

        if not match:
            return None

        frontmatter_str = match.group(1)
        try:
            return yaml.safe_load(frontmatter_str)
        except yaml.YAMLError:
            return None

    @staticmethod
    def _remove_frontmatter(content: str) -> str:
        """Remove frontmatter from content."""
        pattern = r"^---\s*\n.*?\n---\s*\n"
        return re.sub(pattern, "", content, count=1, flags=re.DOTALL)
