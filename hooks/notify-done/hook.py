#!/usr/bin/env python3
"""Notify done hook - Send desktop notification when Claude needs input."""

import json
import os
import platform
import subprocess
import sys


def send_linux_notification(title: str, message: str, urgency: str = "normal"):
    """Send notification on Linux using notify-send."""
    try:
        cmd = ["notify-send", title, message, f"--urgency={urgency}"]

        # Add sound if configured
        if os.environ.get("NOTIFY_DONE_SOUND", "false").lower() == "true":
            cmd.extend(["--hint", "int:transient:1"])

        subprocess.run(cmd, capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def send_macos_notification(title: str, message: str):
    """Send notification on macOS using osascript."""
    try:
        script = f'display notification "{message}" with title "{title}"'
        subprocess.run(
            ["osascript", "-e", script],
            capture_output=True,
            timeout=5,
        )
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def send_windows_notification(title: str, message: str):
    """Send notification on Windows (WSL) using PowerShell."""
    try:
        # Try BurntToast first (if installed)
        ps_script = f"""
        if (Get-Module -ListAvailable -Name BurntToast) {{
            New-BurntToastNotification -Text "{title}", "{message}"
        }} else {{
            [Windows.UI.Notifications.ToastNotificationManager, Windows.UI.Notifications, ContentType = WindowsRuntime] | Out-Null
            $template = [Windows.UI.Notifications.ToastNotificationManager]::GetTemplateContent([Windows.UI.Notifications.ToastTemplateType]::ToastText02)
            $template.GetElementsByTagName("text")[0].AppendChild($template.CreateTextNode("{title}")) | Out-Null
            $template.GetElementsByTagName("text")[1].AppendChild($template.CreateTextNode("{message}")) | Out-Null
            $notifier = [Windows.UI.Notifications.ToastNotificationManager]::CreateToastNotifier("Claude Code")
            $notifier.Show([Windows.UI.Notifications.ToastNotification]::new($template))
        }}
        """
        subprocess.run(
            ["powershell.exe", "-Command", ps_script],
            capture_output=True,
            timeout=10,
        )
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def send_notification(title: str, message: str):
    """Send notification based on platform."""
    system = platform.system().lower()
    urgency = os.environ.get("NOTIFY_DONE_URGENCY", "normal")

    if system == "linux":
        # Check if running in WSL
        try:
            with open("/proc/version") as f:
                if "microsoft" in f.read().lower():
                    if send_windows_notification(title, message):
                        return True
        except FileNotFoundError:
            pass

        return send_linux_notification(title, message, urgency)

    elif system == "darwin":
        return send_macos_notification(title, message)

    elif system == "windows":
        return send_windows_notification(title, message)

    return False


def main():
    """Main hook entry point."""
    # Read input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON input: {e}", file=sys.stderr)
        sys.exit(1)

    # Get configuration
    title = os.environ.get("NOTIFY_DONE_TITLE", "Claude Code")
    message = "Claude is waiting for your input"

    # Try to extract context from input
    cwd = input_data.get("cwd", "")
    if cwd:
        project_name = os.path.basename(cwd)
        message = f"Ready for input in {project_name}"

    # Send notification
    success = send_notification(title, message)

    if success:
        print(f"Notification sent: {message}")
    else:
        print("Could not send notification (notification tools not available)", file=sys.stderr)

    sys.exit(0)


if __name__ == "__main__":
    main()
