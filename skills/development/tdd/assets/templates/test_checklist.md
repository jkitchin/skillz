# Test Coverage Checklist

Comprehensive checklist for ensuring thorough test coverage of your code.

## Feature: [Feature Name]

**Module/Class:** `[module.ClassName]`
**Date:** [YYYY-MM-DD]
**Developer:** [Name]

---

## 1. Happy Path Tests

Tests that verify normal, expected behavior.

- [ ] **Basic functionality** - Core feature works with valid input
- [ ] **Typical use cases** - Common real-world scenarios
- [ ] **Expected outputs** - Returns correct values/types
- [ ] **State changes** - Objects/systems reach correct state
- [ ] **Success conditions** - Operations complete successfully

**Notes:**
```
[Document any important happy path scenarios]
```

---

## 2. Edge Cases

Tests for boundary conditions and unusual inputs.

- [ ] **Empty inputs** - Empty strings, lists, dictionaries
- [ ] **Null/None values** - Missing or undefined values
- [ ] **Minimum values** - Smallest valid input (0, empty, etc.)
- [ ] **Maximum values** - Largest valid input
- [ ] **Boundary values** - Just inside/outside valid range
- [ ] **Single item** - Collections with one element
- [ ] **Large inputs** - Performance with big data
- [ ] **Special characters** - Unicode, escape sequences, etc.
- [ ] **Whitespace** - Leading/trailing spaces, tabs, newlines

**Notes:**
```
[Document edge cases specific to this feature]
```

---

## 3. Error Conditions

Tests that verify proper error handling.

- [ ] **Invalid input types** - Wrong data types
- [ ] **Invalid values** - Out of range, malformed data
- [ ] **Missing required data** - Null references, missing parameters
- [ ] **Constraint violations** - Business rule violations
- [ ] **External failures** - Database, network, file system errors
- [ ] **Error messages** - Clear, helpful error messages
- [ ] **Error types** - Correct exception types raised
- [ ] **Graceful degradation** - System remains stable after errors

**Error scenarios tested:**
```
1.
2.
3.
```

---

## 4. Integration Points

Tests for interactions with other components.

- [ ] **Dependencies** - Calls to external services/modules
- [ ] **API contracts** - Input/output formats match
- [ ] **Data flow** - Data passes correctly between components
- [ ] **State synchronization** - Shared state remains consistent
- [ ] **Event handling** - Events fired/received correctly
- [ ] **Callbacks** - Callback functions invoked properly

**Integration scenarios:**
```
Component A → Component B:
Component A → External Service:
```

---

## 5. State Management

Tests for state transitions and persistence.

- [ ] **Initial state** - Correct starting state
- [ ] **State transitions** - Valid state changes
- [ ] **Invalid transitions** - Prevented/handled correctly
- [ ] **State persistence** - State saved/loaded correctly
- [ ] **Concurrent access** - Thread-safe (if applicable)
- [ ] **State cleanup** - Resources released properly

**State diagram tested:**
```
[Initial] → [State A] → [State B] → [Final]
```

---

## 6. Performance

Tests for performance characteristics (if applicable).

- [ ] **Response time** - Meets performance targets
- [ ] **Memory usage** - No memory leaks
- [ ] **Scalability** - Handles expected load
- [ ] **Resource cleanup** - Files/connections closed
- [ ] **Timeout handling** - Long operations timeout correctly

**Performance requirements:**
```
- Response time: < [X] ms
- Max memory: < [X] MB
- Throughput: > [X] req/sec
```

---

## 7. Security

Security-related tests (if applicable).

- [ ] **Input validation** - Malicious input rejected
- [ ] **SQL injection** - Parameterized queries used
- [ ] **XSS prevention** - Output escaped properly
- [ ] **Authentication** - Auth required where needed
- [ ] **Authorization** - Permission checks enforced
- [ ] **Data sanitization** - Sensitive data handled properly
- [ ] **Rate limiting** - Abuse prevention in place

**Security considerations:**
```
[Note any security-critical code paths]
```

---

## 8. Regression Tests

Tests that prevent previously fixed bugs from recurring.

- [ ] **Bug #[ID]** - [Description]
- [ ] **Bug #[ID]** - [Description]
- [ ] **Bug #[ID]** - [Description]

**Previous issues:**
```
1. Issue: [Description]
   Test: test_bug_123_description()

2. Issue: [Description]
   Test: test_bug_456_description()
```

---

## 9. Test Quality

Quality checks for the tests themselves.

- [ ] **Test independence** - Tests can run in any order
- [ ] **No shared state** - Tests don't affect each other
- [ ] **Fast execution** - Tests run quickly
- [ ] **Clear names** - Test names describe what's tested
- [ ] **Good assertions** - Assertions are specific and clear
- [ ] **Minimal setup** - Tests are easy to understand
- [ ] **Proper cleanup** - Resources cleaned up after tests
- [ ] **No test duplication** - Tests don't redundantly test same thing

---

## 10. Documentation

Test documentation and maintainability.

- [ ] **Docstrings** - All tests have clear descriptions
- [ ] **Comments** - Complex test logic explained
- [ ] **Test data** - Test data is well-documented
- [ ] **Coverage report** - Coverage metrics captured
- [ ] **README** - Test suite documented

---

## Coverage Metrics

**Current Coverage:**
- Line Coverage: [X]%
- Branch Coverage: [X]%
- Function Coverage: [X]%

**Target Coverage:** [X]%

**Coverage Gaps:**
```
File: [filename]
Uncovered lines: [line numbers]
Reason: [why not covered / plan to cover]
```

---

## Test Execution

**How to run tests:**
```bash
# All tests
[command to run all tests]

# Specific test file
[command to run this test file]

# With coverage
[command to run with coverage]
```

**CI/CD Status:**
- [ ] Tests run on every commit
- [ ] Tests run on pull requests
- [ ] Coverage reported automatically
- [ ] Failing tests block merges

---

## Review Sign-off

**Developer:** [Name] - [Date]
**Reviewer:** [Name] - [Date]
**QA:** [Name] - [Date]

**Review Notes:**
```
[Any feedback or additional testing recommendations]
```

---

## Next Steps

- [ ] Add missing tests identified in review
- [ ] Improve coverage in [specific area]
- [ ] Refactor test suite for clarity
- [ ] Update test documentation
- [ ] Add performance benchmarks

**Priority improvements:**
1.
2.
3.
