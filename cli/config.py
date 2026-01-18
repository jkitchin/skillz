"""Configuration management for skillz."""

import os
from pathlib import Path
from typing import Dict, Optional

import yaml


class Config:
    """Configuration manager for skillz."""

    DEFAULT_CONFIG = {
        "personal_skills_dir": "~/.claude/skills",
        "personal_commands_dir": "~/.claude/commands",
        "project_skills_dir": ".claude/skills",
        "project_commands_dir": ".claude/commands",
        "repository_path": None,  # Path to the local clone of skills repository
        "default_target": "personal",  # personal or project
        "default_platform": "claude",  # Default platform: claude, opencode, codex, gemini
        "platforms": {
            "claude": {
                "skills_dir": "~/.claude/skills",
                "commands_dir": "~/.claude/commands",
            },
            "opencode": {
                "skills_dir": "~/.config/opencode/skills",
                "commands_dir": "~/.config/opencode/command",
            },
            "codex": {
                "skills_dir": "~/.codex/skills",
                "commands_dir": "~/.codex/commands",
            },
            "gemini": {
                "skills_dir": "~/.config/gemini/skills",
                "commands_dir": "~/.config/gemini/commands",
            },
        },
    }

    def __init__(self, config_path: Optional[Path] = None):
        """Initialize configuration."""
        self.config_path = config_path or self._default_config_path()
        self.config = self._load_config()

    @staticmethod
    def _default_config_path() -> Path:
        """Get the default configuration file path."""
        return Path.home() / ".config" / "skillz" / "config.yaml"

    def _load_config(self) -> Dict:
        """Load configuration from file or use defaults."""
        if self.config_path.exists():
            with open(self.config_path) as f:
                user_config = yaml.safe_load(f) or {}
                # Deep merge with defaults
                config = self._deep_merge(self.DEFAULT_CONFIG.copy(), user_config)
                return config
        return self.DEFAULT_CONFIG.copy()

    def _deep_merge(self, base: Dict, override: Dict) -> Dict:
        """Deep merge two dictionaries, with override taking precedence."""
        result = base.copy()
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._deep_merge(result[key], value)
            else:
                result[key] = value
        return result

    def save_config(self) -> None:
        """Save current configuration to file."""
        self.config_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.config_path, "w") as f:
            yaml.safe_dump(self.config, f, default_flow_style=False)

    def get_skills_dir(self, target: str = "personal", platform: str = "claude") -> Path:
        """Get the skills directory for a given target and platform."""
        if target == "personal":
            if platform in self.config["platforms"]:
                path = self.config["platforms"][platform]["skills_dir"]
            else:
                path = self.config["personal_skills_dir"]
        else:  # project
            path = self.config["project_skills_dir"]

        return Path(os.path.expanduser(path))

    def get_commands_dir(self, target: str = "personal", platform: str = "claude") -> Path:
        """Get the commands directory for a given target and platform."""
        if target == "personal":
            if platform in self.config["platforms"]:
                path = self.config["platforms"][platform]["commands_dir"]
            else:
                path = self.config["personal_commands_dir"]
        else:  # project
            path = self.config["project_commands_dir"]

        return Path(os.path.expanduser(path))

    def get_repository_path(self) -> Optional[Path]:
        """Get the repository path."""
        repo_path = self.config.get("repository_path")
        if repo_path:
            return Path(os.path.expanduser(repo_path))
        return None

    def set_repository_path(self, path: Path) -> None:
        """Set the repository path."""
        self.config["repository_path"] = str(path)
        self.save_config()
