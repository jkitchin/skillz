"""Main CLI entry point for claude-skills."""

import click
from rich.console import Console

from cli.commands import create, info, install, search, uninstall, update
from cli.commands.list import list_skills

console = Console()


@click.group()
@click.version_option(version="0.1.0", prog_name="claude-skills")
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose output")
@click.pass_context
def cli(ctx, verbose):
    """
    Claude Skills - Manage AI assistant skills and slash commands.

    A CLI tool for installing, managing, and creating skills and commands
    for Claude Code and other LLM platforms.
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


# Register commands
cli.add_command(install.install)
cli.add_command(uninstall.uninstall)
cli.add_command(list_skills)
cli.add_command(search.search)
cli.add_command(info.info)
cli.add_command(update.update)
cli.add_command(create.create)


if __name__ == "__main__":
    cli()
