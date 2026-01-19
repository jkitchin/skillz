"""Example test suite demonstrating Python testing best practices.

This module shows:
- pytest usage and fixtures
- Test organization and naming
- Parametrized tests
- Mocking and patching
- Test coverage patterns
- TDD workflow
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from datetime import datetime
from typing import Any


# Sample code under test (normally in separate module)
class Calculator:
    """Simple calculator for demonstration."""

    def add(self, a: float, b: float) -> float:
        """Add two numbers."""
        return a + b

    def subtract(self, a: float, b: float) -> float:
        """Subtract b from a."""
        return a - b

    def divide(self, a: float, b: float) -> float:
        """Divide a by b."""
        if b == 0:
            raise ZeroDivisionError("Cannot divide by zero")
        return a / b

    def multiply(self, a: float, b: float) -> float:
        """Multiply two numbers."""
        return a * b


class UserManager:
    """Manage users (example with external dependencies)."""

    def __init__(self, database):
        self.database = database

    def create_user(self, username: str, email: str) -> dict[str, Any]:
        """Create new user."""
        if self.database.user_exists(username):
            raise ValueError(f"User {username} already exists")

        user_data = {"username": username, "email": email, "created_at": datetime.now().isoformat()}

        self.database.save_user(user_data)
        return user_data

    def get_user(self, username: str) -> dict[str, Any] | None:
        """Get user by username."""
        return self.database.find_user(username)


# Test Suite for Calculator
class TestCalculator:
    """Test suite for Calculator class.

    Demonstrates:
    - Fixture usage
    - Basic assertions
    - Exception testing
    - Parametrized tests
    """

    @pytest.fixture
    def calc(self) -> Calculator:
        """Provide Calculator instance for tests.

        This fixture is automatically injected into test methods
        that have 'calc' parameter.
        """
        return Calculator()

    def test_addition(self, calc: Calculator):
        """Test basic addition."""
        result = calc.add(2, 3)
        assert result == 5

    def test_addition_negative_numbers(self, calc: Calculator):
        """Test addition with negative numbers."""
        result = calc.add(-1, -1)
        assert result == -2

    def test_addition_floats(self, calc: Calculator):
        """Test addition with floating point numbers."""
        result = calc.add(1.5, 2.5)
        assert result == pytest.approx(4.0)

    def test_subtraction(self, calc: Calculator):
        """Test basic subtraction."""
        result = calc.subtract(5, 3)
        assert result == 2

    def test_multiplication(self, calc: Calculator):
        """Test basic multiplication."""
        result = calc.multiply(4, 3)
        assert result == 12

    def test_division(self, calc: Calculator):
        """Test basic division."""
        result = calc.divide(10, 2)
        assert result == 5.0

    def test_division_by_zero_raises_error(self, calc: Calculator):
        """Test that division by zero raises ZeroDivisionError."""
        with pytest.raises(ZeroDivisionError) as exc_info:
            calc.divide(10, 0)

        assert "Cannot divide by zero" in str(exc_info.value)

    # Parametrized tests - run same test with different inputs
    @pytest.mark.parametrize(
        "a,b,expected",
        [
            (2, 3, 5),
            (0, 0, 0),
            (-1, 1, 0),
            (100, -50, 50),
            (1.5, 2.5, 4.0),
        ],
    )
    def test_add_parametrized(self, calc: Calculator, a: float, b: float, expected: float):
        """Test addition with multiple parameter sets."""
        result = calc.add(a, b)
        assert result == pytest.approx(expected)

    @pytest.mark.parametrize(
        "a,b,expected",
        [
            (10, 2, 5.0),
            (9, 3, 3.0),
            (1, 2, 0.5),
            (-10, 2, -5.0),
        ],
    )
    def test_divide_parametrized(self, calc: Calculator, a: float, b: float, expected: float):
        """Test division with multiple parameter sets."""
        result = calc.divide(a, b)
        assert result == pytest.approx(expected)

    @pytest.mark.parametrize(
        "dividend,divisor",
        [
            (10, 0),
            (0, 0),
            (-5, 0),
        ],
    )
    def test_divide_by_zero_parametrized(self, calc: Calculator, dividend: float, divisor: float):
        """Test division by zero with multiple values."""
        with pytest.raises(ZeroDivisionError):
            calc.divide(dividend, divisor)


# Test Suite with Mocking
class TestUserManager:
    """Test suite for UserManager class.

    Demonstrates:
    - Mocking external dependencies
    - Patching
    - Assertion on mock calls
    """

    @pytest.fixture
    def mock_database(self) -> Mock:
        """Provide mock database for testing."""
        mock_db = Mock()
        mock_db.user_exists.return_value = False
        mock_db.find_user.return_value = None
        return mock_db

    @pytest.fixture
    def user_manager(self, mock_database: Mock) -> UserManager:
        """Provide UserManager with mock database."""
        return UserManager(mock_database)

    def test_create_user_success(self, user_manager: UserManager, mock_database: Mock):
        """Test successful user creation."""
        result = user_manager.create_user("alice", "alice@example.com")

        # Assert return value
        assert result["username"] == "alice"
        assert result["email"] == "alice@example.com"
        assert "created_at" in result

        # Assert database methods were called correctly
        mock_database.user_exists.assert_called_once_with("alice")
        mock_database.save_user.assert_called_once()

    def test_create_user_already_exists(self, user_manager: UserManager, mock_database: Mock):
        """Test that creating existing user raises error."""
        # Configure mock to indicate user exists
        mock_database.user_exists.return_value = True

        with pytest.raises(ValueError) as exc_info:
            user_manager.create_user("alice", "alice@example.com")

        assert "already exists" in str(exc_info.value)

        # Verify save was never called
        mock_database.save_user.assert_not_called()

    def test_get_user_found(self, user_manager: UserManager, mock_database: Mock):
        """Test getting existing user."""
        # Configure mock return value
        expected_user = {"username": "alice", "email": "alice@example.com"}
        mock_database.find_user.return_value = expected_user

        result = user_manager.get_user("alice")

        assert result == expected_user
        mock_database.find_user.assert_called_once_with("alice")

    def test_get_user_not_found(self, user_manager: UserManager, mock_database: Mock):
        """Test getting non-existent user."""
        mock_database.find_user.return_value = None

        result = user_manager.get_user("nonexistent")

        assert result is None


# Advanced Fixtures
@pytest.fixture(scope="session")
def test_data_dir(tmp_path_factory) -> Path:
    """Create temporary directory for test data (session scope).

    This fixture is created once per test session and shared
    across all tests.
    """
    data_dir = tmp_path_factory.mktemp("test_data")
    return data_dir


@pytest.fixture
def sample_file(tmp_path: Path) -> Path:
    """Create sample file for testing (function scope).

    This fixture creates a new file for each test function.
    """
    file_path = tmp_path / "sample.txt"
    file_path.write_text("Sample content\nLine 2\nLine 3")
    return file_path


@pytest.fixture(autouse=True)
def reset_state():
    """Fixture that runs automatically before each test.

    Use autouse=True for setup/teardown that should always run.
    """
    # Setup code here
    yield
    # Teardown code here (runs after test)


# Test using file fixtures
def test_read_file(sample_file: Path):
    """Test reading file content."""
    content = sample_file.read_text()
    assert "Sample content" in content
    assert len(content.splitlines()) == 3


# Patching examples
@patch("datetime.datetime")
def test_with_patched_datetime(mock_datetime):
    """Test with patched datetime."""
    # Configure mock
    mock_datetime.now.return_value = datetime(2025, 1, 1, 12, 0, 0)

    # Code under test would use datetime.now()
    from datetime import datetime as dt
    # Note: In real code, you'd patch where it's used, not where it's defined

    mock_datetime.now.assert_called()


# Test with context manager patching
def test_with_context_manager_patch():
    """Test using patch as context manager."""
    with patch("builtins.open", create=True) as mock_open:
        mock_open.return_value.__enter__.return_value.read.return_value = "mocked content"

        # Code that uses open() would get mocked version
        # with open('file.txt') as f:
        #     content = f.read()

        mock_open.assert_called()


# Markers for test organization
@pytest.mark.slow
def test_slow_operation():
    """Test marked as slow (can be skipped with -m "not slow")."""
    # Simulate slow operation
    import time

    time.sleep(0.1)
    assert True


@pytest.mark.integration
def test_integration():
    """Integration test (can run separately with -m integration)."""
    assert True


@pytest.mark.skip(reason="Not implemented yet")
def test_future_feature():
    """Test that is currently skipped."""
    pass


@pytest.mark.skipif(
    condition=True,  # Replace with actual condition
    reason="Skipped on certain conditions",
)
def test_conditional_skip():
    """Test skipped based on condition."""
    pass


@pytest.mark.xfail(reason="Known bug")
def test_known_bug():
    """Test expected to fail (xfail)."""
    assert False


# Custom assertions and helpers
def assert_valid_email(email: str) -> None:
    """Custom assertion for email validation."""
    assert "@" in email, f"Invalid email: {email}"
    assert "." in email.split("@")[1], f"Invalid email domain: {email}"


def test_custom_assertion():
    """Test using custom assertion."""
    assert_valid_email("user@example.com")


# TDD Example: Write test first
def test_parse_config_from_dict():
    """Test configuration parsing (TDD - test first)."""
    from dataclasses import dataclass

    @dataclass
    class Config:
        host: str
        port: int
        debug: bool

    def parse_config(data: dict) -> Config:
        return Config(host=data["host"], port=int(data["port"]), debug=data.get("debug", False))

    # Now test the implementation
    result = parse_config({"host": "localhost", "port": "8080"})
    assert result.host == "localhost"
    assert result.port == 8080
    assert result.debug is False


# Fixture combinations
@pytest.fixture
def user_data() -> dict[str, str]:
    """Provide test user data."""
    return {"username": "testuser", "email": "test@example.com", "password": "secure123"}


@pytest.fixture
def admin_user_data(user_data: dict[str, str]) -> dict[str, str]:
    """Provide admin user data (depends on user_data fixture)."""
    data = user_data.copy()
    data["is_admin"] = True
    return data


def test_fixture_dependencies(admin_user_data: dict[str, str]):
    """Test using fixture that depends on another fixture."""
    assert admin_user_data["username"] == "testuser"
    assert admin_user_data["is_admin"] is True


# Conftest.py patterns (typically in conftest.py file)
def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line("markers", "slow: mark test as slow")
    config.addinivalue_line("markers", "integration: mark test as integration test")


# Running tests:
# pytest                           # Run all tests
# pytest -v                        # Verbose output
# pytest -k test_addition          # Run tests matching pattern
# pytest -m slow                   # Run tests with 'slow' marker
# pytest -m "not slow"             # Skip slow tests
# pytest --cov=module              # Run with coverage
# pytest --cov-report=html         # Generate HTML coverage report
# pytest -x                        # Stop at first failure
# pytest --lf                      # Run last failed tests
# pytest --pdb                     # Drop into debugger on failure
