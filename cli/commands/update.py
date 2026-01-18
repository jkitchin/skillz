"""Update command for skillz."""

import click
from rich.console import Console

console = Console()


@click.command()
@click.argument("name", required=False)
@click.option("--all", "-a", is_flag=True, help="Update all installed items")
@click.option(
    "--target",
    "-t",
    type=click.Choice(["personal", "project"]),
    help="Target location to update",
)
@click.option(
    "--platform", "-p", default="claude", help="Target platform (claude, opencode, codex, gemini)"
)
@click.option("--dry-run", is_flag=True, help="Preview without making changes")
@click.pass_context
def update(ctx, name, all, target, platform, dry_run):
    """
    Update installed skills and commands.

    NAME is the name of a specific item to update, or use --all to update everything.
    """
    _ = ctx.obj.get("verbose", False)  # Reserved for future use

    if not name and not all:
        console.print("[red]Error: Specify a name or use --all[/red]")
        raise click.Abort()

    if all:
        console.print("[yellow]Updating all items...[/yellow]")
        # TODO: Implement bulk update logic
        # This would:
        # 1. List all installed items
        # 2. Check repository for newer versions
        # 3. Update each item that has changes
        console.print("[yellow]Update --all not yet implemented[/yellow]")
    else:
        console.print(f"[yellow]Updating '{name}'...[/yellow]")
        # TODO: Implement single item update logic
        # This would:
        # 1. Find the installed item
        # 2. Compare with repository version
        # 3. Update if different
        console.print(f"[yellow]Update for '{name}' not yet implemented[/yellow]")

    # For now, suggest using uninstall + install
    console.print("\n[dim]Tip: You can update by running:[/dim]")
    if name:
        console.print(f"[dim]  skillz uninstall {name}[/dim]")
        console.print(f"[dim]  skillz install {name}[/dim]")
    else:
        console.print("[dim]  skillz uninstall <name>[/dim]")
        console.print("[dim]  skillz install <name>[/dim]")
