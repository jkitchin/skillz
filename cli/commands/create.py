"""Create command for claude-skills."""

from pathlib import Path

import click
from rich.console import Console
from rich.prompt import Confirm, Prompt

from cli.config import Config
from cli.utils import validate_description, validate_name

console = Console()


@click.command()
@click.option("--type", "item_type", type=click.Choice(["skill", "command"]), required=True)
@click.option("--name", "-n", help="Name of the skill or command")
@click.option("--interactive", "-i", is_flag=True, default=True, help="Interactive mode")
@click.pass_context
def create(ctx, item_type, name, interactive):
    """
    Create a new skill or command from a template.

    Creates a new skill or command with proper structure and frontmatter.
    """
    verbose = ctx.obj.get("verbose", False)
    config = Config()

    console.print(f"\n[bold]Creating a new {item_type}[/bold]\n")

    # Get name
    if not name:
        if interactive:
            name = Prompt.ask(
                "Enter name (lowercase, hyphens only)",
                default="my-" + item_type,
            )
        else:
            console.print("[red]Error: --name is required in non-interactive mode[/red]")
            raise click.Abort()

    # Validate name
    if not validate_name(name):
        console.print(
            "[red]Error: Invalid name. Use lowercase letters, numbers, and hyphens only[/red]"
        )
        raise click.Abort()

    # Get description
    if interactive:
        description = Prompt.ask("Enter description")
        while not validate_description(description):
            console.print("[red]Description too long (max 1024 chars)[/red]")
            description = Prompt.ask("Enter description")
    else:
        description = f"A {item_type} for Claude Code"

    # Get location
    if interactive:
        use_repo = Confirm.ask(
            "Create in repository? (no = create in current directory)", default=True
        )
    else:
        use_repo = True

    if use_repo:
        repo_path = config.get_repository_path()
        if not repo_path or not repo_path.exists():
            console.print("[red]Error: Repository path not configured[/red]")
            raise click.Abort()

        if item_type == "skill":
            base_path = repo_path / "skills"
        else:
            base_path = repo_path / "commands"
    else:
        base_path = Path.cwd()

    # Create the item
    if item_type == "skill":
        _create_skill(base_path, name, description, verbose)
    else:
        _create_command(base_path, name, description, verbose)

    console.print(f"\n[green]Successfully created {item_type} '{name}'![/green]")


def _create_skill(base_path: Path, name: str, description: str, verbose: bool):
    """Create a new skill directory and SKILL.md."""
    skill_path = base_path / name
    skill_path.mkdir(parents=True, exist_ok=True)

    skill_file = skill_path / "SKILL.md"

    content = f"""---
name: {name}
description: {description}
---

# {name.replace('-', ' ').title()}

{description}

## Usage

When to use this skill:
- TODO: Describe when Claude should use this skill

## Instructions

TODO: Add detailed instructions for how to use this skill

## Examples

TODO: Add examples of using this skill
"""

    skill_file.write_text(content)

    if verbose:
        console.print(f"Created: {skill_file}")


def _create_command(base_path: Path, name: str, description: str, verbose: bool):
    """Create a new command file."""
    cmd_file = base_path / f"{name}.md"
    cmd_file.parent.mkdir(parents=True, exist_ok=True)

    content = f"""---
description: {description}
---

# {name.replace('-', ' ').title()} Command

TODO: Add the command prompt content here.

You can use:
- `$ARGUMENTS` - All arguments passed to the command
- `$1`, `$2`, etc. - Individual positional arguments
- `@filename` - Include file content
- `!command` - Execute bash commands (requires Bash in allowed-tools)

## Example Usage

```
/{name} <arguments>
```
"""

    cmd_file.write_text(content)

    if verbose:
        console.print(f"Created: {cmd_file}")
