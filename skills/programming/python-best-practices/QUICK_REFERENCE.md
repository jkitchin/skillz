# Python Best Practices Quick Reference

A concise reference guide for Python best practices, tooling commands, and common patterns.

## Quick Setup

### Initialize New Project with uv

```bash
# Create virtual environment
uv venv

# Activate environment
source .venv/bin/activate  # Unix/Mac
.venv\Scripts\activate     # Windows

# Install dependencies
uv pip install pandas numpy pytest black ruff mypy

# Create requirements file
uv pip freeze > requirements.txt
```

### Project Structure

```
my_project/
├── pyproject.toml          # Project metadata and tool config
├── README.md               # Project documentation
├── .gitignore             # Git ignore patterns
├── src/
│   └── my_project/        # Source code package
│       ├── __init__.py
│       ├── core.py
│       └── utils.py
├── tests/                 # Test files
│   ├── __init__.py
│   ├── test_core.py
│   └── test_utils.py
└── docs/                  # Documentation
```

## Tooling Commands

### Ruff (Linting and Formatting)

```bash
# Check code for issues
ruff check .

# Auto-fix issues
ruff check --fix .

# Format code
ruff format .

# Check formatting without changes
ruff format --check .

# Run on specific file
ruff check src/my_module.py
```

### Black (Formatting)

```bash
# Format all files
black .

# Check without modifying
black --check .

# Format specific file
black src/my_module.py

# Show diff
black --diff src/my_module.py
```

### Mypy (Type Checking)

```bash
# Check all files
mypy src/

# Check specific file
mypy src/my_module.py

# Ignore missing imports
mypy --ignore-missing-imports src/

# Generate report
mypy --html-report mypy_report src/
```

### Pytest (Testing)

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=src --cov-report=html

# Run specific test file
pytest tests/test_core.py

# Run specific test function
pytest tests/test_core.py::test_function_name

# Run with verbose output
pytest -v

# Run only failed tests
pytest --lf

# Run and stop at first failure
pytest -x
```

## PEP 8 Quick Reference

### Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Function | `snake_case` | `def calculate_total():` |
| Variable | `snake_case` | `user_count = 10` |
| Constant | `UPPER_CASE` | `MAX_SIZE = 100` |
| Class | `PascalCase` | `class UserAccount:` |
| Private | `_leading_underscore` | `def _internal_method():` |
| Module | `snake_case` | `data_processor.py` |
| Package | `snake_case` | `mypackage/` |

### Import Organization

```python
# 1. Standard library imports
import os
import sys
from pathlib import Path

# 2. Third-party imports
import numpy as np
import pandas as pd
import requests

# 3. Local application imports
from myapp.core import process_data
from myapp.utils import validate_input
```

### Line Length and Formatting

```python
# Maximum 79 characters per line (or 99 for code)

# Good: break long function calls
result = some_function(
    argument_one,
    argument_two,
    argument_three,
    keyword_arg=value
)

# Good: break long conditionals
if (condition_one and condition_two
        and condition_three):
    do_something()

# Good: break long strings
message = (
    "This is a very long message that "
    "spans multiple lines for better "
    "readability."
)
```

## Type Hints Cheat Sheet

### Basic Types

```python
from typing import Optional

# Simple types
name: str = "Alice"
age: int = 30
price: float = 19.99
is_active: bool = True

# Collections
numbers: list[int] = [1, 2, 3]
scores: dict[str, int] = {"alice": 90, "bob": 85}
coordinates: tuple[int, int] = (10, 20)
unique_ids: set[int] = {1, 2, 3}

# Optional (can be None)
middle_name: str | None = None  # Python 3.10+
middle_name: Optional[str] = None  # Older syntax
```

### Function Signatures

```python
from collections.abc import Sequence, Callable

# Basic function
def greet(name: str) -> str:
    return f"Hello, {name}"

# Multiple parameters
def add(a: int, b: int) -> int:
    return a + b

# Default values
def repeat(text: str, times: int = 3) -> str:
    return text * times

# No return value
def log_message(message: str) -> None:
    print(message)

# Multiple return types
def divide(a: float, b: float) -> float | None:
    return a / b if b != 0 else None

# Sequence (list, tuple, etc.)
def sum_values(values: Sequence[int]) -> int:
    return sum(values)

# Callable (function as parameter)
def apply(func: Callable[[int], int], value: int) -> int:
    return func(value)
```

### Advanced Types

```python
from typing import TypeVar, Generic, Protocol
from collections.abc import Iterator

# Generic type variable
T = TypeVar('T')

def first_or_none(items: list[T]) -> T | None:
    return items[0] if items else None

# Generic class
class Stack(Generic[T]):
    def __init__(self) -> None:
        self._items: list[T] = []

    def push(self, item: T) -> None:
        self._items.append(item)

    def pop(self) -> T:
        return self._items.pop()

# Protocol (structural typing)
class Drawable(Protocol):
    def draw(self) -> None: ...

def render(obj: Drawable) -> None:
    obj.draw()

# Iterator
def count_up(n: int) -> Iterator[int]:
    for i in range(n):
        yield i
```

## Docstring Templates

### Function Docstring (Google Style)

```python
def calculate_discount(price: float, discount_percent: float) -> float:
    """Calculate discounted price.

    Args:
        price: Original price in dollars
        discount_percent: Discount percentage (0-100)

    Returns:
        Final price after discount

    Raises:
        ValueError: If discount_percent is not in range [0, 100]

    Example:
        >>> calculate_discount(100.0, 20.0)
        80.0
    """
    if not 0 <= discount_percent <= 100:
        raise ValueError("Discount must be between 0 and 100")
    return price * (1 - discount_percent / 100)
```

### Class Docstring

```python
class CustomerAccount:
    """Manage customer account information.

    This class handles customer data including personal information,
    transaction history, and account status.

    Attributes:
        account_id: Unique identifier for the account
        name: Customer full name
        balance: Current account balance in dollars
        is_active: Whether the account is currently active

    Example:
        >>> account = CustomerAccount("A001", "John Doe")
        >>> account.deposit(100.0)
        >>> print(account.balance)
        100.0
    """

    def __init__(self, account_id: str, name: str):
        """Initialize customer account.

        Args:
            account_id: Unique account identifier
            name: Customer's full name
        """
        self.account_id = account_id
        self.name = name
        self.balance = 0.0
        self.is_active = True
```

## Common Patterns

### Context Manager

```python
from contextlib import contextmanager
from typing import Iterator

@contextmanager
def temporary_setting(config: dict, key: str, value: any) -> Iterator[None]:
    """Temporarily change a configuration setting.

    Args:
        config: Configuration dictionary
        key: Setting key to change
        value: Temporary value

    Yields:
        None
    """
    original_value = config.get(key)
    config[key] = value
    try:
        yield
    finally:
        config[key] = original_value

# Usage
with temporary_setting(app_config, 'debug', True):
    run_tests()
```

### Dataclass

```python
from dataclasses import dataclass, field
from datetime import datetime

@dataclass
class User:
    """Represent a user in the system."""
    username: str
    email: str
    created_at: datetime = field(default_factory=datetime.now)
    is_admin: bool = False
    tags: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate after initialization."""
        if '@' not in self.email:
            raise ValueError("Invalid email format")

# Usage
user = User(username="alice", email="alice@example.com")
```

### Property Decorator

```python
class Temperature:
    """Temperature with Celsius/Fahrenheit conversion."""

    def __init__(self, celsius: float):
        self._celsius = celsius

    @property
    def celsius(self) -> float:
        """Get temperature in Celsius."""
        return self._celsius

    @celsius.setter
    def celsius(self, value: float) -> None:
        """Set temperature in Celsius."""
        if value < -273.15:
            raise ValueError("Temperature below absolute zero")
        self._celsius = value

    @property
    def fahrenheit(self) -> float:
        """Get temperature in Fahrenheit."""
        return self._celsius * 9/5 + 32

# Usage
temp = Temperature(25)
print(temp.fahrenheit)  # 77.0
```

### Enum

```python
from enum import Enum, auto

class Status(Enum):
    """Order status enumeration."""
    PENDING = auto()
    PROCESSING = auto()
    SHIPPED = auto()
    DELIVERED = auto()
    CANCELLED = auto()

# Usage
order_status = Status.PENDING
if order_status == Status.PENDING:
    process_order()
```

## Error Handling Patterns

### Custom Exceptions

```python
class ApplicationError(Exception):
    """Base exception for application errors."""
    pass

class ValidationError(ApplicationError):
    """Raised when validation fails."""
    pass

class ResourceNotFoundError(ApplicationError):
    """Raised when requested resource doesn't exist."""
    pass

# Usage
def get_user(user_id: int) -> dict:
    """Retrieve user by ID."""
    user = database.find(user_id)
    if user is None:
        raise ResourceNotFoundError(f"User {user_id} not found")
    return user
```

### Try-Except Patterns

```python
import logging

logger = logging.getLogger(__name__)

def safe_divide(a: float, b: float) -> float | None:
    """Safely divide two numbers."""
    try:
        result = a / b
    except ZeroDivisionError:
        logger.warning(f"Attempted division by zero: {a} / {b}")
        return None
    except TypeError as e:
        logger.error(f"Type error in division: {e}")
        raise
    else:
        logger.debug(f"Division successful: {a} / {b} = {result}")
        return result
    finally:
        logger.debug("Division operation completed")
```

## Testing Patterns

### Basic Test Structure

```python
import pytest
from mymodule import Calculator

class TestCalculator:
    """Test suite for Calculator class."""

    @pytest.fixture
    def calc(self):
        """Provide calculator instance for tests."""
        return Calculator()

    def test_addition(self, calc):
        """Test addition operation."""
        assert calc.add(2, 3) == 5

    def test_division_by_zero(self, calc):
        """Test that division by zero raises exception."""
        with pytest.raises(ZeroDivisionError):
            calc.divide(10, 0)

    @pytest.mark.parametrize("a,b,expected", [
        (2, 3, 5),
        (-1, 1, 0),
        (0, 0, 0),
    ])
    def test_add_parameterized(self, calc, a, b, expected):
        """Test addition with multiple inputs."""
        assert calc.add(a, b) == expected
```

### Mocking

```python
from unittest.mock import Mock, patch
import pytest

def test_api_call():
    """Test function that makes API call."""
    with patch('requests.get') as mock_get:
        # Setup mock response
        mock_get.return_value.json.return_value = {'status': 'ok'}
        mock_get.return_value.status_code = 200

        # Call function under test
        result = fetch_data('https://api.example.com/data')

        # Assertions
        assert result == {'status': 'ok'}
        mock_get.assert_called_once()
```

## Logging Configuration

```python
import logging
from pathlib import Path

def setup_logging(
    level: str = "INFO",
    log_file: Path | None = None,
    format_string: str | None = None
) -> None:
    """Configure application logging.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional file path for log output
        format_string: Custom format string
    """
    if format_string is None:
        format_string = (
            '%(asctime)s - %(name)s - %(levelname)s - '
            '%(filename)s:%(lineno)d - %(message)s'
        )

    handlers: list[logging.Handler] = [logging.StreamHandler()]

    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format=format_string,
        handlers=handlers
    )

# Usage
setup_logging(level="DEBUG", log_file=Path("app.log"))
logger = logging.getLogger(__name__)
logger.info("Application started")
```

## Configuration Files

### pyproject.toml Template

```toml
[project]
name = "my-project"
version = "0.1.0"
description = "Project description"
authors = [{name = "Your Name", email = "you@example.com"}]
requires-python = ">=3.11"
dependencies = [
    "requests>=2.28.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "black>=23.0.0",
    "ruff>=0.1.0",
    "mypy>=1.0.0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.ruff]
line-length = 99
target-version = "py311"
select = ["E", "F", "I", "N", "W", "B", "C4", "UP"]
ignore = ["E501"]

[tool.ruff.per-file-ignores]
"__init__.py" = ["F401"]  # Allow unused imports in __init__.py

[tool.black]
line-length = 99
target-version = ['py311']

[tool.mypy]
python_version = "3.11"
strict = true
warn_return_any = true
warn_unused_configs = true

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = "test_*.py"
python_functions = "test_*"
addopts = "--cov=src --cov-report=html --cov-report=term"
```

## Pre-commit Hooks

### .pre-commit-config.yaml

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.6
    hooks:
      - id: ruff
        args: [--fix]
      - id: ruff-format

  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.5.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files

  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.7.0
    hooks:
      - id: mypy
        additional_dependencies: [types-requests]
```

### Install pre-commit

```bash
uv pip install pre-commit
pre-commit install
pre-commit run --all-files
```

## Checklist for Code Quality

- [ ] All functions have type hints
- [ ] All public functions have docstrings
- [ ] Code passes `ruff check` with no errors
- [ ] Code is formatted with `ruff format` or `black`
- [ ] Types check with `mypy --strict`
- [ ] All tests pass (`pytest`)
- [ ] Test coverage > 80% (`pytest --cov`)
- [ ] No hardcoded values (use constants)
- [ ] No code duplication (DRY)
- [ ] Error handling is appropriate
- [ ] Logging is configured
- [ ] Dependencies are in pyproject.toml
- [ ] README.md is up to date
