# Refactoring with Tests

Safe refactoring process enabled by comprehensive test suite.

## Refactoring Rules

### 1. All Tests Must Be Green
Never start refactoring with failing tests.

### 2. Make Small Changes
One refactoring at a time, not multiple.

### 3. Run Tests After Each Change
Catch breaks immediately.

### 4. Don't Add Features
Refactoring changes structure, not behavior.

### 5. Commit Frequently
Commit after each successful refactoring.

## Common Refactorings

### Extract Function
**Before:**
```python
def process_order(order):
    # Calculate total
    total = 0
    for item in order.items:
        total += item.price * item.quantity
    
    # Apply discount
    if order.customer.is_premium:
        total *= 0.9
    
    return total
```

**After:**
```python
def process_order(order):
    total = calculate_total(order.items)
    total = apply_discount(total, order.customer)
    return total

def calculate_total(items):
    return sum(item.price * item.quantity for item in items)

def apply_discount(total, customer):
    return total * 0.9 if customer.is_premium else total
```

### Rename for Clarity
```python
# Before
def proc(d):
    return d * 2

# After
def double_value(number):
    return number * 2
```

### Remove Duplication
```python
# Before
def calculate_area_rectangle(width, height):
    return width * height

def calculate_area_square(side):
    return side * side

# After
def calculate_area_rectangle(width, height):
    return width * height

def calculate_area_square(side):
    return calculate_area_rectangle(side, side)
```

### Simplify Conditional
```python
# Before
def get_discount(customer):
    if customer.orders > 10:
        if customer.total_spent > 1000:
            return 0.2
        else:
            return 0.1
    else:
        return 0

# After
def get_discount(customer):
    if customer.orders <= 10:
        return 0
    if customer.total_spent > 1000:
        return 0.2
    return 0.1
```

### Extract Variable
```python
# Before
if order.total * 0.1 > order.customer.credit_limit - order.customer.balance:
    raise CreditLimitError()

# After
tax = order.total * 0.1
available_credit = order.customer.credit_limit - order.customer.balance
if tax > available_credit:
    raise CreditLimitError()
```

## Safe Refactoring Process

### Step 1: Ensure Green
```bash
$ run_tests
All tests pass ✅
```

### Step 2: Make One Change
Extract function, rename, remove duplication - pick ONE.

### Step 3: Run Tests
```bash
$ run_tests
All tests pass ✅
```

### Step 4: Commit
```bash
$ git commit -m "Extract calculate_discount function"
```

### Step 5: Repeat
Next refactoring.

## When Tests Break During Refactoring

### Test is Testing Implementation Detail
**Symptom:** Test breaks but behavior unchanged
**Fix:** Update test to test behavior, not implementation

**Example:**
```python
# Bad test (tests implementation)
def test_uses_quicksort():
    assert "quicksort" in str(sort_function.__code__)

# Good test (tests behavior)
def test_returns_sorted_list():
    assert sort_function([3, 1, 2]) == [1, 2, 3]
```

### Test is Testing Behavior
**Symptom:** Test breaks and behavior changed
**Fix:** Fix the code, not the test (you broke behavior)

### Too Many Tests Break
**Symptom:** Many tests fail
**Fix:** Change is too big, revert and make smaller changes

## Red-Green-Refactor Discipline

### In RED Phase
- Write test
- See it fail
- **DON'T** refactor existing code

### In GREEN Phase
- Make test pass
- Minimal code
- **DON'T** refactor yet

### In REFACTOR Phase
- Improve code quality
- Remove duplication
- Better names
- **DON'T** add features
- **DON'T** write new tests

## Refactoring Patterns

### Replace Magic Numbers
```python
# Before
if age >= 18:
    ...

# After
ADULT_AGE = 18
if age >= ADULT_AGE:
    ...
```

### Introduce Parameter Object
```python
# Before
def create_user(first_name, last_name, email, age):
    ...

# After
def create_user(user_data):
    ...
```

### Replace Conditional with Polymorphism
```python
# Before
def calculate_pay(employee):
    if employee.type == "manager":
        return employee.salary * 1.5
    elif employee.type == "engineer":
        return employee.salary * 1.2
    else:
        return employee.salary

# After (with polymorphism)
class Employee:
    def calculate_pay(self):
        return self.salary

class Manager(Employee):
    def calculate_pay(self):
        return self.salary * 1.5

class Engineer(Employee):
    def calculate_pay(self):
        return self.salary * 1.2
```

## Tools for Safe Refactoring

- **IDE Refactoring:** Automated, safe renames and extractions
- **Static Analysis:** Catch issues before running
- **Test Coverage:** Know what's protected
- **Version Control:** Easy to revert
- **CI/CD:** Run full suite automatically

## When to Refactor

### During TDD Cycle
After making test pass, before next test

### When Adding Features
Make the change easy, then make the easy change

### When You See Duplication
Third time you write similar code, refactor

### When Code Smells
Long functions, complex conditionals, unclear names

### When Reading Code
Leave code better than you found it

## When NOT to Refactor

### Tests Are Red
Fix tests first

### Close to Deadline
Stability over perfection

### Code Isn't Tested
Add tests first, then refactor

### You Don't Understand It
Learn it first

## Refactoring Checklist

- [ ] All tests green before starting
- [ ] One small change at a time
- [ ] Run tests after each change
- [ ] All tests still green
- [ ] Commit
- [ ] Next refactoring

Remember: Tests are your safety net. Trust them. Run them frequently. Keep them green.
