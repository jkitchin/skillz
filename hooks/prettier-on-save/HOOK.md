---
name: prettier-on-save
description: Auto-format JavaScript, TypeScript, CSS, and JSON files with Prettier after edits
event: PostToolUse
matcher: Edit|Write
type: command
timeout: 30
---

# Prettier on Save

Automatically format JavaScript, TypeScript, CSS, JSON, and other files with Prettier after Claude edits or creates them.

## Purpose

Ensures consistent code formatting without manual intervention. When Claude writes or edits a file that Prettier supports, this hook automatically formats it.

## Supported File Types

- JavaScript (`.js`, `.jsx`)
- TypeScript (`.ts`, `.tsx`)
- CSS/SCSS/Less (`.css`, `.scss`, `.less`)
- JSON (`.json`)
- Markdown (`.md`)
- HTML (`.html`)
- YAML (`.yaml`, `.yml`)
- GraphQL (`.graphql`)

## Requirements

- Node.js installed
- Prettier installed (`npm install -g prettier` or project-local)

## Configuration

The hook respects your project's Prettier configuration:
- `.prettierrc`
- `.prettierrc.json`
- `prettier.config.js`
- `package.json` (prettier key)

## Behavior

1. Hook receives the file path from the Edit/Write tool
2. Checks if the file extension is supported by Prettier
3. Runs `npx prettier --write` on the file
4. Silent success, reports errors to stderr

## Notes

- Does not block the operation (PostToolUse)
- Fails silently if Prettier is not installed
- Respects `.prettierignore`
