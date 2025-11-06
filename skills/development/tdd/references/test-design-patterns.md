# Test Design Patterns

Comprehensive guide to test design, organization, and anti-patterns.

## What to Test

### Public Interface
**Test:** Methods/functions users call directly
**Why:** This is the contract you're providing

### Business Logic
**Test:** Algorithms, calculations, decision-making
**Why:** Core functionality that drives value

### Edge Cases
**Test:** Boundaries, empty inputs, max values, null/undefined
**Why:** Where bugs hide

### Error Conditions
**Test:** Invalid inputs, exceptions, error handling
**Why:** Verify graceful failure

### Integration Points
**Test:** Where components/services interact
**Why:** Integration failures are common

## What NOT to Test

### Private Implementation
**Don't test:** Private methods, internal state
**Why:** Tests break when refactoring, coupling too tight
**Instead:** Test through public interface

### Third-Party Code
**Don't test:** Library internals
**Why:** Not your responsibility, assumed to work
**Instead:** Test your usage of the library

### Simple Getters/Setters
**Don't test:** Basic property access (unless logic involved)
**Why:** No behavior to verify
**Instead:** Test complex properties

### Framework Code
**Don't test:** Framework internals
**Why:** Framework is tested by its maintainers
**Instead:** Test your code that uses framework

## Test Organization

### One Test Per Behavior
```python
# Good: Focused tests
def test_add_returns_sum_of_positive_numbers():
    assert add(2, 3) == 5

def test_add_returns_negative_when_adding_negatives():
    assert add(-2, -3) == -5

# Bad: Multiple behaviors
def test_add():
    assert add(2, 3) == 5
    assert add(-2, -3) == -5
    assert add(0, 0) == 0
```

### Descriptive Names
Test names should describe: what's tested, scenario, expected result

### Group Related Tests
```python
# Good: Grouped by feature
class TestUserAuthentication:
    def test_valid_credentials_returns_token()
    def test_invalid_credentials_raises_error()
    def test_expired_credentials_refreshes_token()
```

### Separate Test Types
```
tests/
├── unit/          # Fast, isolated
├── integration/   # Components together
└── e2e/           # Full workflows
```

## Test Smells and Anti-Patterns

### Test Smell: Tests Test Multiple Things
**Problem:** If one fails, others don't run
**Fix:** Split into separate tests

### Test Smell: Brittle Tests
**Problem:** Break when implementation changes (even though behavior same)
**Fix:** Test behavior, not implementation

### Test Smell: Slow Tests
**Problem:** Won't be run frequently
**Fix:** Mock slow dependencies, parallel execution

### Test Smell: Tests Depend on Order
**Problem:** Fail when run individually
**Fix:** Make each test independent

### Test Smell: Hidden Dependencies
**Problem:** Tests fail mysteriously
**Fix:** Explicit setup in each test

### Test Smell: Complex Setup
**Problem:** Hard to understand what's being tested
**Fix:** Extract fixtures, simplify

### Test Smell: Testing Implementation
**Problem:** Tests coupled to internals
**Fix:** Test through public API

## Mocking Strategies

### When to Mock

**External Dependencies:**
- Databases
- APIs
- File systems
- Network calls
- Time-dependent code

**Slow Operations:**
- Complex calculations
- Large data processing

**Unpredictable Behavior:**
- Random number generation
- Current time/date
- External service responses

### When NOT to Mock

**Your Own Code:**
Prefer real objects for your own classes

**Simple Objects:**
Value objects, data classes need no mocking

**What You're Testing:**
Never mock the system under test

### Mock vs Stub vs Fake vs Spy

**Mock:** Verifies interactions (method was called with specific args)
**Stub:** Provides responses (returns canned data)
**Fake:** Working implementation (in-memory database)
**Spy:** Records calls for later verification

## Test Case Generation

### Equivalence Partitioning
Group inputs with similar behavior, test one from each group

**Example:** Age validation (0-17: minor, 18-64: adult, 65+: senior)
Test: -1 (invalid), 5 (minor), 30 (adult), 70 (senior)

### Boundary Value Analysis
Test at boundaries where behavior changes

**Example:** Function accepts 1-100
Test: 0, 1, 50, 100, 101

### Error Guessing
Based on experience, guess likely errors

**Common errors:**
- Null/undefined
- Empty collections
- Off-by-one
- Divide by zero
- Missing required fields

## Test Coverage

### Coverage Types

**Line Coverage:** Which lines executed
**Branch Coverage:** Which paths taken
**Function Coverage:** Which functions called

### Coverage Goals

**Target:** 80-90% for critical code
**Don't aim for:** 100% (diminishing returns)
**Focus on:** Business logic, complex code

### Coverage Gaps

**Low coverage areas:**
- Error handling
- Edge cases
- Rarely-used features

**Prioritize:**
- Critical paths first
- Business logic second
- Edge cases third

## Best Practices

### Keep Tests Fast
- Mock slow dependencies
- Use in-memory databases
- Parallel execution
- Only test what's necessary

### Make Tests Readable
- Clear names
- Minimal setup
- Single focus
- Good structure (AAA)

### Independent Tests
- No shared state
- Any order
- Can run individually
- Isolated

### Maintainable Tests
- DRY (extract common setup)
- But not too DRY (prefer clarity)
- Update tests with code
- Delete obsolete tests
