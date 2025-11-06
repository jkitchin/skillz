"""Install command for claude-skills."""

from pathlib import Path

import click
from rich.console import Console

from cli.config import Config
from cli.utils import confirm_action, copy_directory, copy_file
from cli.validator import CommandValidator, SkillValidator

console = Console()


@click.command()
@click.argument("name")
@click.option(
    "--target",
    "-t",
    type=click.Choice(["personal", "project"]),
    default="personal",
    help="Installation target (personal or project)",
)
@click.option(
    "--platform",
    "-p",
    default="claude",
    help="Target platform (claude, codex, gemini)",
)
@click.option("--type", "item_type", type=click.Choice(["skill", "command"]), help="Item type")
@click.option("--force", "-f", is_flag=True, help="Overwrite existing files")
@click.option("--dry-run", is_flag=True, help="Preview without making changes")
@click.pass_context
def install(ctx, name, target, platform, item_type, force, dry_run):
    """
    Install a skill or command.

    NAME is the name of the skill or command to install.
    """
    verbose = ctx.obj.get("verbose", False)
    config = Config()

    # Get repository path
    repo_path = config.get_repository_path()
    if not repo_path or not repo_path.exists():
        console.print(
            "[red]Error: Repository path not configured or does not exist.[/red]"
        )
        console.print("Run: claude-skills config set repository <path>")
        raise click.Abort()

    # Determine item type if not specified
    if not item_type:
        item_type = _detect_item_type(repo_path, name)
        if not item_type:
            console.print(f"[red]Error: Could not find skill or command '{name}'[/red]")
            raise click.Abort()

    # Find source
    if item_type == "skill":
        source_path = _find_skill(repo_path, name)
        if not source_path:
            console.print(f"[red]Error: Skill '{name}' not found in repository[/red]")
            raise click.Abort()

        # Validate skill
        valid, errors = SkillValidator.validate_skill_directory(source_path)
        if not valid:
            console.print(f"[red]Error: Invalid skill:[/red]")
            for error in errors:
                console.print(f"  - {error}")
            raise click.Abort()

        # Get destination
        dest_dir = config.get_skills_dir(target, platform)
        dest_path = dest_dir / name

    else:  # command
        source_path = _find_command(repo_path, name)
        if not source_path:
            console.print(f"[red]Error: Command '{name}' not found in repository[/red]")
            raise click.Abort()

        # Validate command
        valid, errors = CommandValidator.validate_command_file(source_path)
        if not valid:
            console.print(f"[red]Error: Invalid command:[/red]")
            for error in errors:
                console.print(f"  - {error}")
            raise click.Abort()

        # Get destination
        dest_dir = config.get_commands_dir(target, platform)
        dest_path = dest_dir / source_path.name

    if verbose:
        console.print(f"Source: {source_path}")
        console.print(f"Destination: {dest_path}")

    # Check if already exists
    if dest_path.exists() and not force:
        if not confirm_action(
            f"{item_type.capitalize()} '{name}' already exists. Overwrite?", default=False
        ):
            console.print("[yellow]Installation cancelled[/yellow]")
            return

    # Dry run
    if dry_run:
        console.print(f"[blue]Would install {item_type} '{name}' to {dest_path}[/blue]")
        return

    # Install
    dest_dir.mkdir(parents=True, exist_ok=True)

    if item_type == "skill":
        success = copy_directory(source_path, dest_path, force=True)
    else:
        success = copy_file(source_path, dest_path, force=True)

    if success:
        console.print(f"[green]Successfully installed {item_type} '{name}'[/green]")
    else:
        console.print(f"[red]Failed to install {item_type} '{name}'[/red]")


def _detect_item_type(repo_path: Path, name: str) -> str:
    """Detect whether name is a skill or command."""
    if _find_skill(repo_path, name):
        return "skill"
    if _find_command(repo_path, name):
        return "command"
    return None


def _find_skill(repo_path: Path, name: str) -> Path:
    """Find a skill by name in repository."""
    skills_dir = repo_path / "skills"
    if not skills_dir.exists():
        return None

    # Search for skill directory
    for skill_path in skills_dir.rglob("SKILL.md"):
        if skill_path.parent.name == name:
            return skill_path.parent

    return None


def _find_command(repo_path: Path, name: str) -> Path:
    """Find a command by name in repository."""
    commands_dir = repo_path / "commands"
    if not commands_dir.exists():
        return None

    # Search for command file
    for cmd_path in commands_dir.rglob("*.md"):
        if cmd_path.stem == name:
            return cmd_path

    return None
