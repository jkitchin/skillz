# TDD Session Log

**Feature:** [Feature Name]
**Date:** [YYYY-MM-DD]
**Developer:** [Name]
**Duration:** [Start Time] - [End Time]

---

## Session Goals

**Primary Goal:**
```
[What are you building/fixing with TDD today?]
```

**Acceptance Criteria:**
- [ ] Criterion 1
- [ ] Criterion 2
- [ ] Criterion 3

---

## TDD Cycles

### Cycle 1: [Brief description]

**🔴 RED - Write Failing Test**

```
Test: test_[descriptive_name]
File: [test_file.py:line_number]
Status: ❌ FAIL
```

**Test Code:**
```python
def test_descriptive_name():
    """Test description."""
    # Arrange

    # Act

    # Assert
```

**Failure Message:**
```
[Exact failure output from test runner]
```

**Time:** [HH:MM]

---

**🟢 GREEN - Make It Pass**

```
Implementation: [function/class/method]
File: [source_file.py:line_number]
Status: ✅ PASS
```

**Implementation Code:**
```python
def function_name():
    """Implementation."""
    # Minimal code to pass test
```

**All Tests Passing:** ✅ Yes / ❌ No

**Time:** [HH:MM]

---

**🔵 REFACTOR - Improve Code**

**Changes Made:**
- [ ] Improved naming
- [ ] Extracted function
- [ ] Removed duplication
- [ ] Simplified logic
- [ ] Added documentation
- [ ] No refactoring needed

**Refactoring Details:**
```
[Describe what you refactored and why]
```

**All Tests Still Passing:** ✅ Yes / ❌ No

**Commit:** `[commit hash or message]`

**Time:** [HH:MM]

---

### Cycle 2: [Brief description]

**🔴 RED - Write Failing Test**

```
Test: test_[descriptive_name]
File: [test_file.py:line_number]
Status: ❌ FAIL
```

**Test Code:**
```python
def test_descriptive_name():
    """Test description."""
    # Test code here
```

**Failure Message:**
```
[Exact failure output]
```

**Time:** [HH:MM]

---

**🟢 GREEN - Make It Pass**

```
Implementation: [function/class/method]
File: [source_file.py:line_number]
Status: ✅ PASS
```

**Implementation Code:**
```python
# Implementation here
```

**All Tests Passing:** ✅ Yes / ❌ No

**Time:** [HH:MM]

---

**🔵 REFACTOR - Improve Code**

**Changes Made:**
```
[Refactoring details]
```

**All Tests Still Passing:** ✅ Yes / ❌ No

**Commit:** `[commit hash or message]`

**Time:** [HH:MM]

---

### Cycle 3: [Brief description]

[Repeat pattern for each cycle...]

---

## Session Summary

### Tests Written
- Total new tests: [N]
- Total test assertions: [N]
- Test file(s): `[filename(s)]`

### Code Written
- New functions/methods: [N]
- Lines of code: [approx. N]
- Source file(s): `[filename(s)]`

### Coverage
- Starting coverage: [X]%
- Ending coverage: [X]%
- Improvement: [+X]%

### Commits
1. `[hash]` - [commit message]
2. `[hash]` - [commit message]
3. `[hash]` - [commit message]

---

## Insights & Learnings

### What Went Well
```
1. [Something that went well]
2. [Something that went well]
```

### Challenges Faced
```
1. [Challenge] - [How you solved it]
2. [Challenge] - [How you solved it]
```

### TDD Discipline
- [ ] Wrote test first every time
- [ ] Made minimal implementation
- [ ] Refactored when green
- [ ] Committed after each cycle
- [ ] Ran tests frequently

**Deviations from TDD:**
```
[Note any times you deviated and why]
```

---

## Technical Decisions

### Design Choices
```
1. Choice: [What you decided]
   Reason: [Why you decided it]
   Alternative considered: [Other option]

2. Choice: [What you decided]
   Reason: [Why you decided it]
```

### Dependencies Added
- `[package-name]` - [reason for adding]

### Patterns Used
- [ ] Factory pattern
- [ ] Strategy pattern
- [ ] Dependency injection
- [ ] Observer pattern
- [ ] Other: [pattern name]

---

## Code Quality

### Refactorings Performed
1. **[Refactoring type]** in `[file:line]`
   - Before: [description]
   - After: [description]
   - Reason: [why]

### Code Smells Addressed
- [ ] Long function → Extracted methods
- [ ] Duplication → DRY principle applied
- [ ] Magic numbers → Named constants
- [ ] Complex conditional → Simplified
- [ ] Other: [description]

### Tech Debt Created
```
[Note any shortcuts or TODOs for later]
1. TODO: [what needs to be done]
2. TODO: [what needs to be done]
```

---

## Next Session

### Remaining Work
- [ ] [Task to complete]
- [ ] [Task to complete]
- [ ] [Task to complete]

### Next Tests to Write
1. `test_[name]` - [what it will test]
2. `test_[name]` - [what it will test]
3. `test_[name]` - [what it will test]

### Notes for Next Time
```
[Anything you want to remember for the next session]
```

---

## Metrics

| Metric | Value |
|--------|-------|
| Session duration | [X] hours |
| Cycles completed | [N] |
| Tests written | [N] |
| Tests passing | [N] |
| Lines of code | [N] |
| Lines of tests | [N] |
| Test/code ratio | [X:1] |
| Commits | [N] |
| Coverage gained | [+X]% |

---

## Reflection

**Energy Level:** [1-5] ⭐⭐⭐⭐⭐

**Flow State:** [1-5] ⭐⭐⭐⭐⭐

**Overall Session Quality:**
```
[Brief reflection on how the session went, what you learned,
and how you can improve your TDD practice]
```

**TDD Benefits Observed:**
- [ ] Caught bugs early
- [ ] Better design emerged
- [ ] Felt more confident
- [ ] Easy to refactor
- [ ] Clear progress
- [ ] Living documentation

**Improvement Areas:**
```
[What could you do better in your TDD practice?]
```

---

## Resources Referenced

- [Documentation/article link] - [Topic]
- [Stack Overflow link] - [Problem]
- [Reference material] - [Concept]

---

**Session End Time:** [HH:MM]
**Status:** ✅ Complete / ⏸️ Paused / 🔄 In Progress
