---
description: Generate API documentation
allowed-tools: ["Read", "Glob", "Grep"]
argument-hint: <file-or-module-path>
---

# Generate API Documentation

Generate comprehensive API documentation for a module or file.

## Instructions

1. Read the source files specified: $ARGUMENTS
   - If a directory/module, find all public interfaces
   - If a file, document all exports

2. Extract and document:
   - Classes with their methods and properties
   - Functions with signatures
   - Constants and configuration
   - Types and interfaces (for TypeScript)
   - Enums

3. For each item, include:
   - Name and signature
   - Description
   - Parameters with types
   - Return type and description
   - Examples
   - Related items

## Output Format

```markdown
# API Reference: <module-name>

## Overview
<Brief description of the module>

## Classes

### `ClassName`

<Description>

#### Constructor

```<language>
ClassName(param1: Type, param2: Type)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| param1 | `Type` | Description |
| param2 | `Type` | Description |

#### Methods

##### `methodName(params) -> ReturnType`

<Description>

**Parameters:**
- `param`: Description

**Returns:** Description

**Example:**
```<language>
instance.methodName(value)
```

## Functions

### `functionName(params) -> ReturnType`

<Description>

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| param1 | `Type` | - | Description |

**Returns:** `ReturnType` - Description

**Example:**
```<language>
result = functionName(value)
```

## Constants

| Name | Type | Value | Description |
|------|------|-------|-------------|
| CONST_NAME | `Type` | `value` | Description |

## Types

### `TypeName`

```<language>
<type definition>
```

<Description>
```

Module to document: $ARGUMENTS
