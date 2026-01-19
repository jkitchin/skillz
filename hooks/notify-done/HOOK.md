---
name: notify-done
description: Send desktop notification when Claude needs input or completes a task
event: Stop
type: command
timeout: 10
---

# Notify Done

Sends a desktop notification when Claude Code stops and is waiting for user input. Perfect for long-running tasks where you've switched to another window.

## Purpose

Get notified when:
- Claude completes a task and needs your input
- Claude encounters an error and needs guidance
- A long-running operation finishes

## Supported Platforms

- **Linux**: Uses `notify-send` (libnotify)
- **macOS**: Uses `osascript` (AppleScript)
- **Windows WSL**: Uses `powershell.exe` with BurntToast or native toast

## Requirements

### Linux
```bash
# Debian/Ubuntu
sudo apt install libnotify-bin

# Fedora
sudo dnf install libnotify

# Arch
sudo pacman -S libnotify
```

### macOS
No additional installation required (uses built-in osascript)

### Windows WSL
PowerShell notifications work out of the box

## Configuration

```bash
# Custom notification title
export NOTIFY_DONE_TITLE="Claude Code"

# Custom notification sound (Linux only)
export NOTIFY_DONE_SOUND="true"

# Notification urgency: low, normal, critical (Linux only)
export NOTIFY_DONE_URGENCY="normal"
```

## Behavior

1. Hook triggered when Claude stops (Stop event)
2. Detects operating system
3. Sends appropriate desktop notification
4. Exits 0 (never blocks)

## Notes

- **Non-blocking** - notification is fire-and-forget
- Falls back gracefully if notification tools aren't available
- Can be combined with other Stop hooks
