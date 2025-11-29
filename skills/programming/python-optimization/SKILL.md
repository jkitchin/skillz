---
name: python-optimization
description: Expert guidance for mathematical optimization in Python - systematic problem classification, library selection (scipy, pyomo, cvxpy, GEKKO), solver configuration, and implementation patterns for LP, QP, NLP, MIP, convex, and global optimization problems
allowed-tools: ["*"]
---

# Python Optimization Expert

## Overview

Master mathematical optimization in Python through systematic problem analysis, appropriate library selection, and robust implementation. This skill guides you through classifying optimization problems, choosing the right solver, implementing efficient solutions, and troubleshooting convergence issues across linear programming, nonlinear programming, mixed-integer problems, and global optimization.

**Core value:** Transform optimization problems from mathematical formulation to working Python code with the most appropriate library and solver for your specific problem type.

## When to Use

Use this skill when:
- Minimizing or maximizing objective functions
- Parameter estimation and curve fitting
- Portfolio optimization and resource allocation
- Production planning and scheduling
- Engineering design optimization
- Control and dynamic optimization
- Solving constrained or unconstrained optimization problems
- Choosing between scipy, pyomo, cvxpy, or other optimization libraries

**Don't use when:**
- Simple root finding (use scipy.optimize.root instead)
- Linear algebra only (use numpy/scipy.linalg)
- Symbolic mathematics needed (use sympy)

## Problem Analysis Framework

### Step 1: Extract Problem Components

**Always identify:**

1. **Objective function**: What to minimize or maximize?
   - Example: minimize cost, maximize profit, minimize error

2. **Decision variables**: What can be changed?
   - Continuous (real numbers) or discrete (integers)?
   - How many variables?

3. **Constraints**: What restrictions exist?
   - Equality constraints: h(x) = 0
   - Inequality constraints: g(x) ≥ 0
   - Bounds: lower ≤ x ≤ upper

4. **Data/parameters**: What values are given?
   - Constants, measurements, input data

### Step 2: Classify Problem Type

**Linear Programming (LP)**
- Linear objective function
- Linear constraints only
- Example: max c'x subject to Ax ≤ b

**Quadratic Programming (QP)**
- Quadratic objective function
- Linear constraints
- Example: min (1/2)x'Qx + c'x subject to Ax ≤ b

**Nonlinear Programming (NLP)**
- Nonlinear objective or constraints
- Smooth, differentiable functions
- Example: min f(x) subject to g(x) ≥ 0

**Mixed-Integer Programming (MIP/MILP/MINLP)**
- Some or all variables must be integers
- Combinatorial optimization
- Example: scheduling, assignment problems

**Convex Optimization**
- Convex objective function
- Convex feasible region
- Global optimum guaranteed

**Global Optimization**
- Multiple local minima suspected
- Non-convex problems
- Requires global search strategies

**Least Squares**
- Minimize sum of squared residuals
- Parameter estimation, curve fitting
- Special structure exploited by dedicated solvers

**Dynamic Optimization / Optimal Control**
- Time-dependent variables
- Differential equations as constraints
- Trajectory optimization

### Step 3: Characterize Problem Scale

**Small problems**: < 100 variables
- Most methods work well
- scipy.optimize is sufficient

**Medium problems**: 100-1,000 variables
- Need efficient algorithms
- Consider specialized solvers

**Large problems**: > 1,000 variables
- Requires sparse matrix techniques
- Use pyomo or specialized solvers
- Exploit problem structure

**Other characteristics:**
- Are derivatives available? Analytical or numerical?
- Smooth or non-smooth?
- Sparse or dense?
- Special structure (separable, etc.)?

## Library Selection Guide

### scipy.optimize - Default Choice

**Use when:**
- Small to medium problems (<1,000 variables)
- Unconstrained or smooth constrained NLP
- Quick prototyping and development
- Least squares / curve fitting
- Python-native solution preferred

**Strengths:**
- Part of scipy (widely available)
- Easy to use, good documentation
- Multiple algorithms available
- No additional dependencies

**Limitations:**
- Not ideal for large-scale problems
- Limited support for mixed-integer
- Less sophisticated than commercial solvers

### pyomo - Algebraic Modeling

**Use when:**
- Large-scale problems (>1,000 variables)
- Mixed-integer programming
- Need multiple solver backends
- Algebraic modeling preferred
- Production environments

**Strengths:**
- Scales to very large problems
- Supports many external solvers
- Clean mathematical notation
- Active development and community

**Limitations:**
- Requires external solvers for many problem types
- Steeper learning curve
- Additional installation complexity

### cvxpy - Convex Optimization

**Use when:**
- Problem is (or might be) convex
- Want automatic convexity verification
- Prefer disciplined convex programming
- Portfolio optimization, control

**Strengths:**
- Automatic convexity detection
- Mathematical notation syntax
- Multiple solver backends
- Excellent for convex problems

**Limitations:**
- Only for convex problems
- Less flexible for non-convex
- Learning curve for DCP rules

### GEKKO - Dynamic Optimization

**Use when:**
- Black-box functions (no derivatives)
- Dynamic optimization / optimal control
- Differential equations involved
- Process control, trajectory optimization

**Strengths:**
- Handles differential equations naturally
- Derivative-free optimization
- Good for dynamic systems
- Time-dependent problems

**Limitations:**
- Specialized use case
- Less common than scipy
- Requires GEKKO server

### Comparison Table

| Library | Best For | Problem Types | Scale |
|---------|----------|---------------|-------|
| scipy.optimize | General NLP, prototyping | LP, QP, NLP, LS | Small-Medium |
| pyomo | Large-scale, MIP | LP, MILP, NLP, MINLP | Any |
| cvxpy | Convex problems | LP, QP, SOCP, SDP | Medium-Large |
| GEKKO | Dynamic optimization | NLP, DAE, MPC | Small-Medium |

## scipy.optimize - Method Selection

### Unconstrained Optimization

**`minimize(f, x0, method='BFGS', jac=gradient)`**
- **Use for:** Smooth unconstrained problems
- **Requires:** Gradient (or uses finite differences)
- **Pros:** Fast convergence, robust
- **Cons:** Can get stuck in local minima

**`minimize(f, x0, method='L-BFGS-B', bounds=bounds)`**
- **Use for:** Large-scale with box constraints
- **Pros:** Memory efficient, handles bounds
- **Cons:** Only box constraints (no general constraints)

**`minimize(f, x0, method='Nelder-Mead')`**
- **Use for:** Derivative-free, non-smooth
- **Pros:** No derivatives needed
- **Cons:** Slow convergence, unreliable for large problems

**`minimize(f, x0, method='Powell')`**
- **Use for:** Derivative-free, better than Nelder-Mead
- **Pros:** More robust than Nelder-Mead
- **Cons:** Still derivative-free (slower)

### Constrained Optimization

**`minimize(f, x0, method='SLSQP', constraints=cons)`**
- **Use for:** General constraints (equality + inequality)
- **Pros:** Handles all constraint types
- **Cons:** Less robust than trust-constr

**`minimize(f, x0, method='trust-constr', constraints=cons)`**
- **Use for:** Difficult constrained problems
- **Pros:** Most robust, can use Hessian
- **Cons:** Slower than SLSQP

**`minimize(f, x0, method='COBYLA', constraints=cons)`**
- **Use for:** Derivative-free with constraints
- **Pros:** No derivatives needed
- **Cons:** Slow, less accurate

### Special-Purpose Functions

**`least_squares(residual, x0, bounds=bounds)`**
- **Use for:** Nonlinear least squares
- **Pros:** Exploits sum-of-squares structure
- **Cons:** Only for least squares problems

**`curve_fit(model, xdata, ydata, p0=initial)`**
- **Use for:** Parameter estimation, curve fitting
- **Pros:** Simple interface for fitting
- **Cons:** Just a wrapper for least_squares

**`differential_evolution(f, bounds)`**
- **Use for:** Global optimization, multiple local minima
- **Pros:** Finds global optimum, derivative-free
- **Cons:** Slow, must specify bounds

**`linprog(c, A_ub, b_ub)`**
- **Use for:** Linear programming only
- **Pros:** Specialized for LP, efficient
- **Cons:** Only linear problems

## Implementation Patterns

### Pattern 1: Unconstrained Optimization

```python
from scipy.optimize import minimize
import numpy as np

def objective(x):
    """
    Objective function to minimize.

    Example: Rosenbrock function
    f(x) = sum((1 - x[i])^2 + 100*(x[i+1] - x[i]^2)^2)
    """
    return sum((1 - x[:-1])**2 + 100*(x[1:] - x[:-1]**2)**2)

def gradient(x):
    """
    Analytical gradient - ALWAYS provide if possible!
    Dramatically improves convergence speed and reliability.
    """
    grad = np.zeros_like(x)
    grad[:-1] = -2*(1 - x[:-1]) - 400*x[:-1]*(x[1:] - x[:-1]**2)
    grad[1:] += 200*(x[1:] - x[:-1]**2)
    return grad

# Initial guess
x0 = np.array([0.5, 0.5])

# Optimize
result = minimize(objective, x0, method='BFGS', jac=gradient)

# Check results
print(f"Success: {result.success}")
print(f"Message: {result.message}")
print(f"Optimal x: {result.x}")
print(f"Optimal value: {result.fun}")
print(f"Iterations: {result.nit}")
print(f"Function evaluations: {result.nfev}")
```

### Pattern 2: Constrained Optimization

```python
from scipy.optimize import minimize
import numpy as np

def objective(x):
    """Objective: minimize x[0]^2 + x[1]^2"""
    return x[0]**2 + x[1]**2

def constraint_ineq(x):
    """Inequality constraint: x[0] + x[1] - 1 >= 0"""
    return x[0] + x[1] - 1

def constraint_eq(x):
    """Equality constraint: x[0]^2 + x[1]^2 - 1 = 0"""
    return x[0]**2 + x[1]**2 - 1

# Define constraints (can have multiple)
constraints = [
    {'type': 'ineq', 'fun': constraint_ineq},  # g(x) >= 0
    {'type': 'eq', 'fun': constraint_eq}       # h(x) = 0
]

# Bounds: [(lower, upper), ...], use None for no bound
bounds = [(-2, 2), (-2, 2)]

# Initial guess (should be feasible if possible)
x0 = np.array([0.5, 0.5])

# Optimize
result = minimize(
    objective,
    x0,
    method='SLSQP',
    constraints=constraints,
    bounds=bounds,
    options={'ftol': 1e-9, 'maxiter': 1000}
)

print(f"Optimal x: {result.x}")
print(f"Optimal value: {result.fun}")

# Verify constraints are satisfied
print(f"\nConstraint verification:")
print(f"Inequality (should be >= 0): {constraint_ineq(result.x):.6f}")
print(f"Equality (should be ~0): {constraint_eq(result.x):.6e}")
```

### Pattern 3: Parameter Estimation (Curve Fitting)

```python
from scipy.optimize import curve_fit
import numpy as np

def model(x, a, b, c):
    """
    Model function: f(x; parameters)

    Example: Exponential decay + offset
    """
    return a * np.exp(-b * x) + c

# Generate synthetic data (in practice, this is your measured data)
x_data = np.linspace(0, 10, 50)
y_true = 2.5 * np.exp(-0.3 * x_data) + 0.5
y_data = y_true + 0.1 * np.random.randn(len(x_data))  # Add noise

# Fit model to data
# p0 = initial parameter guess (optional but recommended)
params, covariance = curve_fit(
    model,
    x_data,
    y_data,
    p0=[2, 0.5, 0]  # Initial guess: [a, b, c]
)

# Extract parameters and uncertainties
a_opt, b_opt, c_opt = params
uncertainties = np.sqrt(np.diag(covariance))

print(f"Optimal parameters:")
print(f"  a = {a_opt:.3f} ± {uncertainties[0]:.3f}")
print(f"  b = {b_opt:.3f} ± {uncertainties[1]:.3f}")
print(f"  c = {c_opt:.3f} ± {uncertainties[2]:.3f}")

# Calculate R-squared
y_pred = model(x_data, *params)
ss_res = np.sum((y_data - y_pred)**2)
ss_tot = np.sum((y_data - np.mean(y_data))**2)
r_squared = 1 - (ss_res / ss_tot)
print(f"\nR² = {r_squared:.4f}")
```

### Pattern 4: Least Squares (Residuals)

```python
from scipy.optimize import least_squares
import numpy as np

def residuals(params, x, y):
    """
    Residual function for least squares.

    Returns: array of residuals (model - data)
    """
    a, b, c = params
    model = a * np.exp(-b * x) + c
    return model - y

# Data
x_data = np.linspace(0, 10, 50)
y_data = 2.5 * np.exp(-0.3 * x_data) + 0.5 + 0.1*np.random.randn(50)

# Initial guess
p0 = [2, 0.5, 0]

# Optimize with bounds
result = least_squares(
    residuals,
    p0,
    args=(x_data, y_data),
    bounds=([0, 0, -np.inf], [10, 5, np.inf])  # Lower and upper bounds
)

print(f"Success: {result.success}")
print(f"Optimal parameters: {result.x}")
print(f"Final cost (sum of squared residuals): {result.cost}")
print(f"Optimality: {result.optimality}")  # Should be small
```

### Pattern 5: Linear Programming

```python
from scipy.optimize import linprog
import numpy as np

"""
Linear Programming Problem:

Maximize:   3*x1 + 2*x2
Subject to: x1 + x2 <= 4
            2*x1 + x2 <= 5
            x1, x2 >= 0

Note: linprog minimizes, so negate objective for maximization
"""

# Objective coefficients (negate for maximization)
c = [-3, -2]

# Inequality constraints: A_ub @ x <= b_ub
A_ub = np.array([
    [1, 1],   # x1 + x2 <= 4
    [2, 1]    # 2*x1 + x2 <= 5
])
b_ub = np.array([4, 5])

# Bounds: x1, x2 >= 0
bounds = [(0, None), (0, None)]

# Solve
result = linprog(
    c,
    A_ub=A_ub,
    b_ub=b_ub,
    bounds=bounds,
    method='highs'  # Default, most efficient
)

print(f"Success: {result.success}")
print(f"Optimal x: {result.x}")
print(f"Optimal value (minimization): {result.fun}")
print(f"Optimal value (maximization): {-result.fun}")
```

### Pattern 6: Global Optimization

```python
from scipy.optimize import differential_evolution
import numpy as np

def objective(x):
    """
    Function with multiple local minima.

    Example: Rastrigin function
    """
    A = 10
    return A*len(x) + sum(x**2 - A*np.cos(2*np.pi*x))

# MUST specify bounds for each variable
bounds = [(-5.12, 5.12), (-5.12, 5.12)]

# Global optimization
result = differential_evolution(
    objective,
    bounds,
    strategy='best1bin',
    maxiter=1000,
    popsize=15,
    tol=1e-7,
    atol=1e-9,
    workers=-1,  # Use all CPU cores
    polish=True,  # Local refinement at the end
    seed=42  # For reproducibility
)

print(f"Success: {result.success}")
print(f"Global minimum x: {result.x}")
print(f"Global minimum value: {result.fun}")
print(f"Function evaluations: {result.nfev}")

# For comparison, try local optimization from random start
from scipy.optimize import minimize
x0 = np.random.uniform(-5.12, 5.12, 2)
local_result = minimize(objective, x0, method='BFGS')
print(f"\nLocal optimization from random start:")
print(f"Local minimum value: {local_result.fun}")
print(f"(Global is better: {result.fun < local_result.fun})")
```

### Pattern 7: pyomo (Large-Scale / Mixed-Integer)

```python
import pyomo.environ as pyo

# Create model
model = pyo.ConcreteModel()

# Decision variables
model.x = pyo.Var(domain=pyo.NonNegativeReals)
model.y = pyo.Var(domain=pyo.Integers, bounds=(0, 10))
model.z = pyo.Var([1, 2, 3], domain=pyo.Binary)  # Indexed binary variables

# Parameters (optional)
model.a = pyo.Param(initialize=2.5)
model.b = pyo.Param([1, 2, 3], initialize={1: 1.0, 2: 2.0, 3: 3.0})

# Objective function
model.obj = pyo.Objective(
    expr=model.x**2 + model.y + sum(model.z[i] for i in [1,2,3]),
    sense=pyo.minimize  # or pyo.maximize
)

# Constraints
model.con1 = pyo.Constraint(expr=model.x + model.y >= 5)
model.con2 = pyo.Constraint(expr=2*model.x - model.y <= 10)
model.con3 = pyo.Constraint(expr=sum(model.z[i] for i in [1,2,3]) <= 2)

# Solve
# For continuous: 'ipopt' (nonlinear) or 'glpk' (linear)
# For integer: 'cbc', 'glpk', 'gurobi' (commercial)
solver = pyo.SolverFactory('ipopt')
results = solver.solve(model, tee=True)  # tee=True shows solver output

# Check if optimal
if results.solver.termination_condition == pyo.TerminationCondition.optimal:
    print("\nOptimal solution found!")
    print(f"x = {pyo.value(model.x):.4f}")
    print(f"y = {pyo.value(model.y):.4f}")
    print(f"z = {[pyo.value(model.z[i]) for i in [1,2,3]]}")
    print(f"Objective value = {pyo.value(model.obj):.4f}")
else:
    print(f"Solver status: {results.solver.termination_condition}")
```

### Pattern 8: cvxpy (Convex Optimization)

```python
import cvxpy as cp
import numpy as np

"""
Portfolio optimization example:
Minimize variance subject to minimum expected return
"""

# Problem data
n = 4  # Number of assets
mu = np.array([0.12, 0.10, 0.07, 0.03])  # Expected returns
Sigma = np.array([  # Covariance matrix
    [0.10, 0.03, 0.01, 0.00],
    [0.03, 0.08, 0.02, 0.01],
    [0.01, 0.02, 0.05, 0.01],
    [0.00, 0.01, 0.01, 0.02]
])

# Variables
w = cp.Variable(n)  # Portfolio weights

# Objective: minimize variance
objective = cp.Minimize(cp.quad_form(w, Sigma))

# Constraints
constraints = [
    cp.sum(w) == 1,           # Weights sum to 1
    w >= 0,                   # No short selling
    mu @ w >= 0.10            # Minimum expected return
]

# Solve
problem = cp.Problem(objective, constraints)
problem.solve()

# Results
print(f"Status: {problem.status}")
if problem.status == 'optimal':
    print(f"Optimal portfolio weights: {w.value}")
    print(f"Portfolio variance: {problem.value:.6f}")
    print(f"Portfolio std dev: {np.sqrt(problem.value):.4f}")
    print(f"Expected return: {(mu @ w.value):.4f}")
else:
    print(f"Problem could not be solved: {problem.status}")
```

## Best Practices

### 1. Always Provide Analytical Gradients

**Why:** 10-100x faster convergence, more reliable

```python
def objective(x):
    return x[0]**2 + x[1]**2

def gradient(x):
    """Analytical gradient - much better than finite differences"""
    return np.array([2*x[0], 2*x[1]])

# Good: Uses analytical gradient
result = minimize(objective, x0, method='BFGS', jac=gradient)

# Okay but slower: Uses finite differences
result = minimize(objective, x0, method='BFGS')
```

### 2. Verify Gradients

```python
from scipy.optimize import check_grad

# Check gradient accuracy
error = check_grad(objective, gradient, x0)
print(f"Gradient error: {error}")  # Should be < 1e-5

# For constraints
for con in constraints:
    if 'jac' in con:
        error = check_grad(con['fun'], con['jac'], x0)
        print(f"Constraint gradient error: {error}")
```

### 3. Scale Variables

**Problem:** Variables with very different magnitudes
```python
# Bad: x1 ~ 1e-6, x2 ~ 1e6
x = [1e-6, 1e6]
```

**Solution:** Normalize to similar scales
```python
# Good: Both variables ~ O(1)
x_scaled = [x[0] * 1e6, x[1] / 1e6]

def objective_scaled(x_scaled):
    x = [x_scaled[0] / 1e6, x_scaled[1] * 1e6]  # Convert back
    return original_objective(x)
```

### 4. Provide Good Initial Guess

```python
# Bad: Random or arbitrary starting point
x0 = np.zeros(n)

# Good: Physically reasonable or feasible starting point
x0 = np.array([1.0, 2.0, 0.5])  # Based on problem knowledge

# Better: Try multiple starting points for non-convex problems
initial_guesses = [
    np.array([1.0, 2.0]),
    np.array([-1.0, -2.0]),
    np.random.randn(2)
]

results = [minimize(objective, x0) for x0 in initial_guesses]
best_result = min(results, key=lambda r: r.fun)
```

### 5. Validate Solutions

```python
def validate_solution(result, constraints, bounds):
    """Comprehensive solution validation"""

    # Check optimization succeeded
    assert result.success, f"Optimization failed: {result.message}"

    # Check inequality constraints: g(x) >= 0
    for i, con in enumerate([c for c in constraints if c['type'] == 'ineq']):
        val = con['fun'](result.x)
        assert val >= -1e-6, f"Inequality constraint {i} violated: {val:.2e}"

    # Check equality constraints: h(x) = 0
    for i, con in enumerate([c for c in constraints if c['type'] == 'eq']):
        val = con['fun'](result.x)
        assert abs(val) < 1e-6, f"Equality constraint {i} violated: {val:.2e}"

    # Check bounds
    if bounds is not None:
        for i, (lb, ub) in enumerate(bounds):
            if lb is not None:
                assert result.x[i] >= lb - 1e-6, f"Lower bound {i} violated"
            if ub is not None:
                assert result.x[i] <= ub + 1e-6, f"Upper bound {i} violated"

    print("✓ All constraints satisfied")
    return True

# Use it
validate_solution(result, constraints, bounds)
```

### 6. Handle Convergence Issues

**Common problems and solutions:**

```python
# Problem: "Desired error not necessarily achieved"
# Solution 1: Relax tolerances
result = minimize(
    objective, x0,
    method='BFGS',
    options={'ftol': 1e-6, 'gtol': 1e-5}  # Less strict
)

# Solution 2: Scale variables
x_scaled = x / typical_scale
result = minimize(objective_scaled, x_scaled)

# Solution 3: Try different method
result = minimize(objective, x0, method='trust-constr')

# Solution 4: Better initial guess
x0 = feasible_starting_point()

# Problem: Stuck in local minimum
# Solution: Global optimization
result = differential_evolution(objective, bounds)

# Problem: Constraints incompatible
# Solution: Check constraint feasibility
from scipy.optimize import fsolve
# Try to find any feasible point
def feasibility(x):
    violations = []
    for con in constraints:
        if con['type'] == 'eq':
            violations.append(con['fun'](x))
    return violations

x_feasible = fsolve(feasibility, x0)
```

## Common Problem Types - Quick Reference

### 1. Unconstrained: minimize f(x)

```python
minimize(f, x0, method='BFGS', jac=grad)
```

### 2. Box constraints: minimize f(x) subject to lb ≤ x ≤ ub

```python
minimize(f, x0, method='L-BFGS-B', bounds=bounds)
```

### 3. General constraints: minimize f(x) s.t. g(x) ≥ 0, h(x) = 0

```python
cons = [{'type': 'ineq', 'fun': g}, {'type': 'eq', 'fun': h}]
minimize(f, x0, method='SLSQP', constraints=cons)
```

### 4. Least squares: minimize ||r(x)||²

```python
least_squares(residual, x0, bounds=(lb, ub))
```

### 5. Curve fitting: fit model to data

```python
params, cov = curve_fit(model, xdata, ydata, p0=initial)
```

### 6. Linear program: minimize c'x s.t. Ax ≤ b

```python
linprog(c, A_ub=A, b_ub=b, bounds=bounds)
```

### 7. Mixed-integer: some variables must be integers

```python
# Use pyomo
model.x = pyo.Var(domain=pyo.Integers)
solver = pyo.SolverFactory('cbc')
solver.solve(model)
```

### 8. Global: multiple local minima

```python
differential_evolution(f, bounds)
```

### 9. Portfolio optimization

```python
# cvxpy
w = cp.Variable(n)
objective = cp.Minimize(cp.quad_form(w, Sigma))
constraints = [cp.sum(w) == 1, w >= 0, mu @ w >= target_return]
problem = cp.Problem(objective, constraints)
problem.solve()
```

### 10. Parameter estimation with uncertainties

```python
params, cov = curve_fit(model, xdata, ydata)
uncertainties = np.sqrt(np.diag(cov))
```

## Response Format

When solving an optimization problem, always:

### 1. Restate the Problem

"Based on your description, this is a **[problem type]** problem:

**Objective**: minimize/maximize f(x) = ...
**Variables**: x ∈ R^n (or x ∈ Z^n for integers)
**Constraints**:
- [List inequality constraints]
- [List equality constraints]
- [List bounds]"

### 2. Explain Library/Method Choice

"Since this problem has [characteristics], I recommend **[library]** with **[method]**:
- [Reason 1]
- [Reason 2]
- [Alternative if this doesn't work]"

### 3. Provide Complete Implementation

```python
# Complete, runnable code with comments
# Include imports, data, functions, optimization, validation
```

### 4. Display and Interpret Results

"**Results**:
- Optimal x: [values]
- Optimal value: [value]
- Constraints satisfied: ✓ / ✗
- Convergence: [iterations, function evaluations]

**Interpretation**: [What the solution means in context]"

### 5. Suggest Improvements if Needed

"**If convergence issues arise:**
- Try [alternative method]
- Adjust [parameters]
- Consider [problem reformulation]"

## Advanced Topics

### Multi-Objective Optimization

For problems with **multiple conflicting objectives** (minimize cost AND maximize performance), use the dedicated **`python-multiobjective-optimization`** skill which covers:

- Pareto optimality and Pareto fronts
- Evolutionary algorithms (NSGA-II, NSGA-III, MOEA/D)
- Scalarization methods (weighted sum, ε-constraint)
- Libraries: pymoo, platypus, DEAP
- Pareto analysis and decision making
- Visualization techniques

**Quick example** - Weighted sum scalarization with scipy:

```python
from scipy.optimize import minimize

def scalarize_objectives(x, weights):
    """Weighted sum approach for multi-objective"""
    f1 = objective1(x)
    f2 = objective2(x)
    return weights[0] * f1 + weights[1] * f2

# Solve for different weight combinations to trace Pareto front
weights_list = [(1, 0), (0.75, 0.25), (0.5, 0.5), (0.25, 0.75), (0, 1)]
pareto_front = []

for weights in weights_list:
    result = minimize(
        lambda x: scalarize_objectives(x, weights),
        x0,
        method='SLSQP',
        constraints=constraints
    )
    if result.success:
        pareto_front.append({
            'x': result.x,
            'f1': objective1(result.x),
            'f2': objective2(result.x)
        })
```

**Note:** Weighted sum only finds convex portions of Pareto front. For comprehensive multiobjective optimization, see `python-multiobjective-optimization` skill.

### Robust Optimization

```python
def worst_case_objective(x, uncertainty_scenarios):
    """Minimize worst-case performance"""
    return max(objective(x, scenario) for scenario in uncertainty_scenarios)

result = minimize(worst_case_objective, x0, args=(scenarios,))
```

### Stochastic Optimization

```python
def expected_objective(x, n_samples=1000):
    """Minimize expected value over random parameters"""
    samples = np.random.normal(mu, sigma, (n_samples, n_params))
    return np.mean([objective(x, sample) for sample in samples])

result = minimize(expected_objective, x0)
```

## Troubleshooting Guide

| Problem | Symptoms | Solutions |
|---------|----------|-----------|
| Slow convergence | Many iterations, slow progress | Provide analytical gradient, scale variables, better x0 |
| Local minimum | Different x0 gives different solutions | Use `differential_evolution`, try multiple starts |
| Constraint violations | Solution doesn't satisfy constraints | Check constraint formulation, ensure feasible x0, tighten tolerances |
| "Desired error not achieved" | Warning in result | Relax tolerances, scale variables, try different method |
| Numerical errors | NaN or Inf values | Scale variables, avoid extreme values in objective |
| Memory errors | Large problems crash | Use sparse matrices, L-BFGS-B, or pyomo with sparse solvers |

## Installation

```bash
# Core optimization
pip install scipy

# For large-scale and mixed-integer
pip install pyomo

# Install solvers for pyomo
conda install -c conda-forge ipopt  # Nonlinear
conda install -c conda-forge glpk   # Linear/Integer
pip install cplex  # Commercial, requires license

# For convex optimization
pip install cvxpy

# For dynamic optimization
pip install gekko
```

## Additional Resources

- **SciPy Optimization Tutorial:** https://docs.scipy.org/doc/scipy/tutorial/optimize.html
- **Pyomo Documentation:** https://pyomo.readthedocs.io/
- **CVXPY Tutorial:** https://www.cvxpy.org/tutorial/
- **Convex Optimization (Boyd & Vandenberghe):** https://web.stanford.edu/~boyd/cvxbook/
- **Numerical Optimization (Nocedal & Wright):** Classic textbook

## Related Skills

- `python-multiobjective-optimization` - Multiobjective optimization with Pareto fronts (NSGA-II, NSGA-III, pymoo)
- `pycse` - Regression analysis with confidence intervals, ODE solving
- `python-best-practices` - Code quality and testing for optimization code
- `python-ase` - Atomic structure optimization in materials science
