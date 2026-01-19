---
description: Generate documentation for code
allowed-tools: ["Read", "Glob", "Grep"]
argument-hint: <file-path or function-name>
---

# Generate Documentation

Generate documentation (docstrings, comments, or markdown) for code.

## Instructions

1. Read the code specified: $ARGUMENTS
   - If a file, document all public functions/classes
   - If a specific function/class, document just that

2. Detect the language and use appropriate documentation format:
   - **Python**: Google-style or NumPy-style docstrings
   - **JavaScript/TypeScript**: JSDoc comments
   - **Go**: GoDoc comments
   - **Rust**: Rustdoc comments
   - **Other**: Language-appropriate format

3. Document:
   - Purpose/description
   - Parameters with types and descriptions
   - Return values
   - Exceptions/errors that may be raised
   - Usage examples where helpful

## Output Format

For Python:
```python
def function_name(param1: Type, param2: Type) -> ReturnType:
    """Short description of function.

    Longer description if needed, explaining the purpose
    and behavior of the function.

    Args:
        param1: Description of param1.
        param2: Description of param2.

    Returns:
        Description of return value.

    Raises:
        ExceptionType: When this error occurs.

    Example:
        >>> function_name(value1, value2)
        expected_result
    """
```

For JavaScript/TypeScript:
```javascript
/**
 * Short description of function.
 *
 * @param {Type} param1 - Description of param1
 * @param {Type} param2 - Description of param2
 * @returns {ReturnType} Description of return value
 * @throws {ErrorType} When this error occurs
 * @example
 * functionName(value1, value2)
 */
```

Code to document: $ARGUMENTS
