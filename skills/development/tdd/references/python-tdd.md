# Python TDD Guide

Comprehensive guide to Test-Driven Development in Python using pytest and unittest.

## Python Testing Ecosystem

### pytest (Recommended)
Modern, powerful, widely adopted testing framework.

**Advantages:**
- Simple assert statements (no need for self.assertEqual)
- Automatic test discovery
- Fixtures for setup/teardown
- Parametrized tests
- Rich plugin ecosystem
- Better failure reporting

**Installation:**
```bash
pip install pytest pytest-cov
```

### unittest
Python standard library testing framework.

**Advantages:**
- No installation needed (standard library)
- Similar to xUnit frameworks (JUnit, NUnit)
- Class-based structure

**Use when:**
- Can't add dependencies
- Working with existing unittest codebase

### doctest
Tests embedded in docstrings.

**Use for:**
- Documentation examples
- Simple function tests
- Not recommended as primary testing approach

## pytest Basics

### Simple Test Example

```python
# test_calculator.py
from calculator import Calculator

def test_add_positive_numbers():
    calc = Calculator()
    assert calc.add(2, 3) == 5

def test_add_with_zero():
    calc = Calculator()
    assert calc.add(5, 0) == 5
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific file
pytest test_calculator.py

# Run specific test
pytest test_calculator.py::test_add_positive_numbers

# Run with verbose output
pytest -v

# Run with coverage
pytest --cov=src --cov-report=html
```

## Assertions

### pytest Assertions

```python
# Equality
assert result == expected

# Inequality
assert result != unexpected

# Boolean
assert is_valid
assert not is_empty

# Containment
assert item in collection
assert key in dictionary

# Exceptions
import pytest

with pytest.raises(ValueError):
    function_that_raises()

with pytest.raises(ValueError, match="invalid input"):
    function_with_specific_message()

# Approximate equality (floats)
assert result == pytest.approx(0.333, rel=1e-3)
```

### unittest Assertions

```python
import unittest

class TestCalculator(unittest.TestCase):
    def test_add(self):
        calc = Calculator()
        self.assertEqual(calc.add(2, 3), 5)

    def test_raises_error(self):
        with self.assertRaises(ValueError):
            function_that_raises()
```

## Fixtures (pytest)

### Basic Fixture

```python
import pytest

@pytest.fixture
def calculator():
    """Provide a Calculator instance for tests."""
    return Calculator()

def test_add(calculator):
    assert calculator.add(2, 3) == 5

def test_subtract(calculator):
    assert calculator.subtract(5, 3) == 2
```

### Setup and Teardown

```python
@pytest.fixture
def database():
    """Create database connection, yield it, then clean up."""
    db = Database.connect()
    yield db
    db.close()

def test_insert(database):
    database.insert("test data")
    assert database.count() == 1
```

### Fixture Scopes

```python
@pytest.fixture(scope="function")  # Default: new instance per test
def func_scope():
    pass

@pytest.fixture(scope="class")  # Shared across test class
def class_scope():
    pass

@pytest.fixture(scope="module")  # Shared across test module
def module_scope():
    pass

@pytest.fixture(scope="session")  # Shared across entire test session
def session_scope():
    pass
```

### conftest.py

Share fixtures across multiple test files:

```python
# tests/conftest.py
import pytest

@pytest.fixture
def database():
    """Available to all test files in this directory."""
    db = Database.connect()
    yield db
    db.close()
```

## Parametrized Tests

Test same function with multiple inputs:

```python
import pytest

@pytest.mark.parametrize("a,b,expected", [
    (2, 3, 5),
    (0, 0, 0),
    (-1, 1, 0),
    (100, 200, 300),
])
def test_add(a, b, expected):
    assert add(a, b) == expected

@pytest.mark.parametrize("value", [
    "",
    None,
    [],
    {},
])
def test_empty_values(value):
    assert is_empty(value)
```

## Mocking

### unittest.mock

```python
from unittest.mock import Mock, patch, MagicMock

# Mock object
mock_db = Mock()
mock_db.query.return_value = ["result1", "result2"]

# Patch function
with patch('mymodule.expensive_function') as mock_func:
    mock_func.return_value = 42
    result = code_that_calls_expensive_function()
    assert result == 42

# Patch as decorator
@patch('mymodule.external_api_call')
def test_api_integration(mock_api):
    mock_api.return_value = {"status": "ok"}
    result = function_using_api()
    assert result["status"] == "ok"

# Mock class
with patch('mymodule.Database') as MockDatabase:
    mock_instance = MockDatabase.return_value
    mock_instance.query.return_value = []
    # Use mock instance
```

### pytest-mock

```bash
pip install pytest-mock
```

```python
def test_with_mocker(mocker):
    mock_func = mocker.patch('mymodule.function')
    mock_func.return_value = 42
    assert code_using_function() == 42
```

## Testing Async Code

### pytest-asyncio

```bash
pip install pytest-asyncio
```

```python
import pytest

@pytest.mark.asyncio
async def test_async_function():
    result = await async_function()
    assert result == expected_value

@pytest.mark.asyncio
async def test_async_with_fixture(async_database):
    await async_database.insert("data")
    count = await async_database.count()
    assert count == 1
```

## Test Organization

### Directory Structure

```
project/
├── src/
│   ├── __init__.py
│   ├── calculator.py
│   └── database.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_calculator.py
│   └── test_database.py
├── pytest.ini
└── setup.py
```

### pytest.ini Configuration

```ini
[pytest]
testpaths = tests
python_files = test_*.py *_test.py
python_functions = test_*
addopts = -v --strict-markers --cov=src
markers =
    slow: marks tests as slow
    integration: marks tests as integration tests
    unit: marks tests as unit tests
```

### Running Specific Test Types

```bash
# Run only unit tests
pytest -m unit

# Run all except slow tests
pytest -m "not slow"

# Run integration tests
pytest -m integration
```

## Coverage

### pytest-cov

```bash
# Run with coverage
pytest --cov=src

# Generate HTML report
pytest --cov=src --cov-report=html

# Show missing lines
pytest --cov=src --cov-report=term-missing

# Fail if coverage below threshold
pytest --cov=src --cov-fail-under=80
```

### .coveragerc Configuration

```ini
[run]
source = src
omit = */tests/*

[report]
exclude_lines =
    pragma: no cover
    def __repr__
    raise AssertionError
    raise NotImplementedError
    if __name__ == .__main__.:
```

## TDD Workflow Example

### RED Phase - Write Failing Test

```python
# test_string_utils.py
def test_reverse_string():
    assert reverse("hello") == "olleh"
```

Run: `pytest test_string_utils.py`
Result: ❌ `NameError: name 'reverse' is not defined`

### GREEN Phase - Minimal Implementation

```python
# string_utils.py
def reverse(s):
    return s[::-1]
```

Run: `pytest test_string_utils.py`
Result: ✅ Test passes

### REFACTOR Phase

```python
# string_utils.py
def reverse(s: str) -> str:
    """Reverse a string.

    Args:
        s: String to reverse

    Returns:
        Reversed string
    """
    return s[::-1]
```

Run: `pytest test_string_utils.py`
Result: ✅ Still passes

## Best Practices

### Test Naming
```python
# Good: Descriptive, says what's being tested
def test_add_two_positive_numbers_returns_sum():
    pass

def test_divide_by_zero_raises_value_error():
    pass

# Bad: Vague, doesn't describe what's tested
def test_add():
    pass

def test_1():
    pass
```

### Test Independence
```python
# Good: Each test is independent
def test_add():
    calc = Calculator()
    assert calc.add(2, 3) == 5

def test_subtract():
    calc = Calculator()  # Fresh instance
    assert calc.subtract(5, 3) == 2

# Bad: Tests depend on order
calc = Calculator()

def test_add():
    result = calc.add(2, 3)
    assert result == 5

def test_subtract():
    # Assumes test_add ran first!
    result = calc.subtract(10, 5)
    assert result == 5
```

### One Assertion Per Test (Guideline)
```python
# Good: One logical assertion
def test_add_returns_correct_sum():
    assert add(2, 3) == 5

def test_add_is_commutative():
    assert add(2, 3) == add(3, 2)

# Acceptable: Multiple assertions testing same behavior
def test_parse_config():
    config = parse("config.yaml")
    assert config["host"] == "localhost"
    assert config["port"] == 5432
    assert config["database"] == "testdb"

# Bad: Testing multiple unrelated behaviors
def test_everything():
    assert add(2, 3) == 5
    assert subtract(5, 3) == 2
    assert multiply(2, 3) == 6
    # If one fails, all fail
```

## Common Patterns

### Testing Exceptions
```python
import pytest

def test_divide_by_zero():
    with pytest.raises(ZeroDivisionError):
        divide(1, 0)

def test_invalid_input():
    with pytest.raises(ValueError, match="must be positive"):
        process(-1)
```

### Testing File Operations
```python
import pytest
from pathlib import Path

@pytest.fixture
def temp_file(tmp_path):
    """pytest provides tmp_path fixture."""
    file = tmp_path / "test.txt"
    file.write_text("content")
    return file

def test_read_file(temp_file):
    content = read_file(temp_file)
    assert content == "content"
```

### Testing with Database
```python
@pytest.fixture
def db():
    database = Database(":memory:")  # SQLite in-memory
    database.create_tables()
    yield database
    database.close()

def test_insert(db):
    db.insert(User(name="Alice"))
    users = db.query_all()
    assert len(users) == 1
    assert users[0].name == "Alice"
```

## Resources

**pytest documentation:** https://docs.pytest.org/
**unittest documentation:** https://docs.python.org/3/library/unittest.html
**pytest-cov:** https://pytest-cov.readthedocs.io/
**pytest-mock:** https://pytest-mock.readthedocs.io/
