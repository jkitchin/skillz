# Universal TDD Principles

Language-agnostic TDD guide applicable to any programming language.

## Finding Testing Tools for Your Language

**Search patterns:**
- "[language] testing framework"
- "[language] unit testing"
- "[language] xUnit" (JUnit, NUnit, etc.)

**Common frameworks by language:**
- **JavaScript:** Jest, Mocha, Jasmine
- **Java:** JUnit, TestNG
- **C#:** NUnit, xUnit, MSTest
- **Ruby:** RSpec, Minitest
- **Go:** testing (standard library), Testify
- **Rust:** cargo test (built-in)
- **PHP:** PHPUnit
- **Swift:** XCTest

## Universal Test Structure

### Setup-Exercise-Verify-Teardown

**Setup (Arrange):**
- Create test data
- Configure environment
- Initialize objects

**Exercise (Act):**
- Execute code under test
- Single operation being tested

**Verify (Assert):**
- Check results
- Verify expectations
- Confirm behavior

**Teardown (Cleanup):**
- Release resources
- Reset state
- Close connections

## Common Testing Concepts

### Test Suites
Group of related tests that run together.

### Test Cases
Individual test functions/methods.

### Assertions
Statements that verify expected outcomes.

### Fixtures
Reusable test setup and teardown code.

### Mocking
Replacing real dependencies with controlled substitutes.

### Parametrized Tests
Running same test with multiple inputs.

### Test Coverage
Measuring which code is executed by tests.

## Test Naming Strategies

**Pattern 1:** `test_<function>_<scenario>_<expected>`
```
test_add_positive_numbers_returns_sum
test_divide_by_zero_raises_exception
```

**Pattern 2:** `should_<behavior>_when_<condition>`
```
should_return_true_when_valid_input
should_throw_error_when_null
```

**Pattern 3:** BDD style - descriptive sentences
```
"it returns the sum when adding positive numbers"
"it raises an error when dividing by zero"
```

## Language-Agnostic TDD Example

### Iteration 1: Basic functionality

**RED:** Write test for simplest behavior
```
test: fibonacci(0) should return 0
status: FAIL (function doesn't exist)
```

**GREEN:** Implement minimal code
```
fibonacci(n):
    return 0
status: PASS
```

**REFACTOR:** Nothing to improve yet
```
status: PASS
```

### Iteration 2: Next case

**RED:** Test next behavior
```
test: fibonacci(1) should return 1
status: FAIL (returns 0, expected 1)
```

**GREEN:** Make it pass
```
fibonacci(n):
    if n == 0: return 0
    if n == 1: return 1
status: PASS (both tests)
```

**REFACTOR:** Good enough
```
status: PASS
```

### Iteration 3: General case

**RED:** Test recursive case
```
test: fibonacci(5) should return 5
status: FAIL (returns 1, expected 5)
```

**GREEN:** Implement recursion
```
fibonacci(n):
    if n <= 1: return n
    return fibonacci(n-1) + fibonacci(n-2)
status: PASS (all tests)
```

**REFACTOR:** Add memoization for performance
```
# Add caching logic
status: PASS (all tests still pass)
```

## Build Tool Integration

Most languages have build tools that run tests:

- **JavaScript:** `npm test`, `yarn test`
- **Python:** `pytest`, `python -m unittest`
- **Java:** `mvn test`, `gradle test`
- **Ruby:** `rake test`, `bundle exec rspec`
- **Go:** `go test ./...`
- **Rust:** `cargo test`
- **C#:** `dotnet test`

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Setup [Language]
        uses: [language-setup-action]
      - name: Install dependencies
        run: [install-command]
      - name: Run tests
        run: [test-command]
```

## Universal TDD Principles

1. **Test First** - Always write test before implementation
2. **Minimal Implementation** - Simplest code to pass
3. **Refactor When Green** - Only improve when tests pass
4. **Fast Tests** - Tests should run quickly
5. **Independent Tests** - No dependencies between tests
6. **Readable Tests** - Tests are documentation
7. **One Assertion Focus** - Test one behavior
8. **Avoid Testing Privates** - Test public interface
9. **Mock External Dependencies** - Don't hit real APIs/databases
10. **Keep Test Suite Fast** - Run frequently

## Common Anti-Patterns (Universal)

❌ Implementation before tests
❌ Testing private methods
❌ Tests depending on execution order
❌ Slow tests
❌ Complex test setup
❌ Mock excessive
❌ Brittle tests (break when implementation changes)
❌ Not running tests frequently
❌ Skipping refactor step
❌ Testing framework code

## Resources by Language

Search for "[language] TDD tutorial" or "[language] testing best practices" for language-specific guidance.
