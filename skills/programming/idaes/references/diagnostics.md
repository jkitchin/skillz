# Model Diagnostics and Debugging

Comprehensive guide to diagnosing and fixing IDAES model issues.

## Overview

The IDAES DiagnosticsToolbox provides powerful tools for identifying and resolving model problems. This guide covers how to use diagnostics to debug structural issues, numerical problems, and solver failures.

## DiagnosticsToolbox

### Basic Usage

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

# Create diagnostics for your model
dt = DiagnosticsToolbox(m)

# Run structural diagnostics
dt.report_structural_issues()

# Run numerical diagnostics
dt.report_numerical_issues()
```

## Structural Diagnostics

Structural issues are problems with the model formulation itself, independent of variable values.

### Degrees of Freedom Analysis

```python
# Check DOF for entire model
from idaes.core.util.model_statistics import degrees_of_freedom

dof = degrees_of_freedom(m)
print(f"Model DOF: {dof}")

# Should be 0 for well-posed problem
# DOF > 0: under-specified (need more equations or fix more variables)
# DOF < 0: over-specified (conflicting constraints)

# Check individual units
dof_heater = degrees_of_freedom(m.fs.heater)
print(f"Heater DOF: {dof_heater}")
```

### Finding Unfixed Variables

```python
# List all unfixed variables
unfixed_vars = []
for v in m.component_data_objects(Var, active=True, descend_into=True):
    if not v.fixed and v.value is not None:
        unfixed_vars.append(v.name)

print(f"Unfixed variables: {len(unfixed_vars)}")
```

### Finding Fixed Variables

```python
# List all fixed variables
fixed_vars = []
for v in m.component_data_objects(Var, active=True, descend_into=True):
    if v.fixed:
        fixed_vars.append((v.name, v.value))

print(f"Fixed variables: {len(fixed_vars)}")
for name, value in fixed_vars[:10]:  # Show first 10
    print(f"  {name} = {value}")
```

### Detecting Structural Singularities

Structural singularities indicate fundamental problems in model structure:

```python
dt = DiagnosticsToolbox(m)

# Get structural issues
dt.report_structural_issues()
```

**Common structural issues:**
- **Redundant constraints:** Multiple constraints that are linearly dependent
- **Unconstrained variables:** Variables that don't appear in any active constraint
- **Unused variables:** Variables that are present but not affecting the solution

**Example output:**
```
====================================================================================
Model Statistics

    Activated Blocks: 5 (Deactivated: 0)
    Free Variables in Activated Constraints: 150 (External: 0)
        Free Variables with only lower bounds: 0
        Free Variables with only upper bounds: 0
        Free Variables with upper and lower bounds: 50

    Activated Equality Constraints: 150 (Deactivated: 0)
    Activated Inequality Constraints: 0 (Deactivated: 0)
    Activated Objectives: 0 (Deactivated: 0)

WARNING: 1 Constraint(s) not included in structural analysis

Next steps:

Call report_structural_issues() to identify structural issues in the model
Call display_constraints_with_large_residuals() to identify poorly satisfied constraints
```

### Constraint Analysis

```python
# Display constraints with residuals
dt.display_constraints_with_large_residuals(tol=1e-4)

# Get list of constraints
constraints = [c.name for c in m.component_data_objects(Constraint, active=True)]
print(f"Total constraints: {len(constraints)}")
```

## Numerical Diagnostics

Numerical issues arise from poor variable scaling or ill-conditioning.

### Jacobian Analysis

```python
# Compute Jacobian and analyze conditioning
dt.report_numerical_issues()
```

**What this checks:**
- Badly scaled variables (very large or very small values)
- Badly scaled constraints (Jacobian entries with large range)
- Near-parallel constraints (linearly dependent rows)
- Variables with zero derivatives (potential errors)

**Example output:**
```
====================================================================================
Numerical Diagnostics

Jacobian Condition Number: 1.23e+12  # High = poorly conditioned

Variables with extreme values:
    fs.heater.control_volume.properties_out[0].pressure: 1.01e+05
    fs.heater.control_volume.properties_out[0].flow_mol: 1.00e+02

Badly scaled variables (Jacobian):
    fs.flash.liq_outlet.flow_mol[0]: scale = 1.0, magnitude = 1.5e-08
    fs.flash.vap_outlet.temperature[0]: scale = 1.0, magnitude = 3.7e+02

Constraints with large residuals:
    fs.flash.eq_phase_equilibrium[0]: residual = 1.2e-03

Next steps:
Call display_variables_with_extreme_jacobian_values() for more details
Apply scaling using idaes.core.util.scaling
```

### Checking for NaN/Inf Values

```python
# Check for non-finite values
for v in m.component_data_objects(Var, active=True, descend_into=True):
    if v.value is not None:
        if not np.isfinite(v.value):
            print(f"Non-finite value: {v.name} = {v.value}")
```

### Variable Bounds

```python
# Check variables at bounds
at_bounds = []
for v in m.component_data_objects(Var, active=True, descend_into=True):
    if not v.fixed and v.value is not None:
        if v.lb is not None and abs(v.value - v.lb) < 1e-6:
            at_bounds.append((v.name, "lower", v.value))
        if v.ub is not None and abs(v.value - v.ub) < 1e-6:
            at_bounds.append((v.name, "upper", v.value))

print(f"Variables at bounds: {len(at_bounds)}")
for name, bound_type, value in at_bounds[:10]:
    print(f"  {name} at {bound_type} bound: {value}")
```

## SVD Analysis

Singular Value Decomposition can identify numerical rank deficiency.

```python
from idaes.core.util.model_diagnostics import svd_analysis

# Perform SVD on Jacobian
svd_analysis(m)
```

**Interpretation:**
- **Small singular values:** Indicate near-singularity
- **Singular value ratio (max/min):** Condition number estimate
- **Associated constraints/variables:** Show what causes ill-conditioning

## Degeneracy Hunter

Specialized tool for finding degenerate constraints and variables.

```python
from idaes.core.util.model_diagnostics import DegeneracyHunter

# Create hunter
dh = DegeneracyHunter(m)

# Find issues
dh.check_residuals(tol=1e-5)
dh.check_variable_bounds(tol=1e-6)
```

## Common Problems and Solutions

### Problem 1: DOF ≠ 0

**Diagnosis:**
```python
dof = degrees_of_freedom(m)
print(f"DOF = {dof}")
```

**If DOF > 0 (under-specified):**
- Not enough equations or fixed variables
- Solution: Fix more input variables or add constraints

```python
# Find unfixed variables that should be fixed
m.fs.heater.inlet.flow_mol[0].fix(100)
m.fs.heater.heat_duty[0].fix(10000)
```

**If DOF < 0 (over-specified):**
- Too many equations or fixed variables
- Solution: Unfix some variables or remove redundant constraints

```python
# Unfix a variable
m.fs.heater.outlet.temperature[0].unfix()

# Or deactivate a constraint
m.fs.custom_constraint.deactivate()
```

### Problem 2: Solver Fails with "Infeasible Problem"

**Diagnosis:**
```python
dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
dt.display_constraints_with_large_residuals(tol=1e-3)
```

**Possible causes:**
- Conflicting specifications (e.g., fixing both heat duty and outlet temperature)
- Physical impossibility (e.g., cooling below ambient)
- Incorrect property bounds

**Solutions:**
1. Check for over-specification
2. Verify physical feasibility
3. Relax some constraints
4. Check variable bounds

### Problem 3: Poor Numerical Conditioning

**Diagnosis:**
```python
dt.report_numerical_issues()

# Check Jacobian condition number
# Values > 1e10 indicate poor conditioning
```

**Solutions:**
```python
from idaes.core.util import scaling as iscale

# Apply scaling
iscale.calculate_scaling_factors(m)

# Check scaling was applied
badly_scaled = list(iscale.badly_scaled_var_generator(m))
print(f"Badly scaled variables: {len(badly_scaled)}")
```

See `scaling.md` for detailed scaling guidance.

### Problem 4: Constraint Violations

**Diagnosis:**
```python
# After solve, check constraint satisfaction
dt.display_constraints_with_large_residuals(tol=1e-5)
```

**Solutions:**
1. Tighten solver tolerances
```python
solver = get_solver('ipopt')
solver.options['tol'] = 1e-8
solver.options['constr_viol_tol'] = 1e-8
```

2. Improve initialization
3. Apply better scaling

### Problem 5: Variables at Bounds

**Diagnosis:**
```python
# Check which variables hit bounds
for v in m.fs.heater.component_data_objects(Var, descend_into=True):
    if not v.fixed and v.value is not None:
        if v.lb is not None and abs(v.value - v.lb) < 1e-6:
            print(f"{v.name} at lower bound: {v.value}")
        if v.ub is not None and abs(v.value - v.ub) < 1e-6:
            print(f"{v.name} at upper bound: {v.value}")
```

**Solutions:**
1. Check if bounds are too restrictive
2. Verify physical correctness of bounds
3. May indicate model needs reformulation

## Advanced Diagnostics

### Constraint by Constraint Analysis

```python
# Evaluate each constraint's residual
for c in m.component_data_objects(Constraint, active=True, descend_into=True):
    if c.body is not None:
        try:
            residual = value(c.body) - value(c.upper)
            if abs(residual) > 1e-4:
                print(f"{c.name}: residual = {residual:.2e}")
        except:
            print(f"{c.name}: cannot evaluate")
```

### Variable Scaling Analysis

```python
from idaes.core.util import scaling as iscale

# Check scaling factors
for v in m.component_data_objects(Var, descend_into=True):
    sf = iscale.get_scaling_factor(v)
    if sf is not None and v.value is not None:
        scaled_value = v.value * sf
        if abs(scaled_value) > 100 or abs(scaled_value) < 0.01:
            print(f"{v.name}: value={v.value}, scale={sf}, scaled={scaled_value}")
```

### Jacobian Sparsity

```python
# Visualize Jacobian structure
import matplotlib.pyplot as plt
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

nlp = PyomoNLP(m)
jac = nlp.evaluate_jacobian()

plt.figure(figsize=(10, 10))
plt.spy(jac, markersize=1)
plt.title("Jacobian Sparsity Pattern")
plt.show()
```

## Debugging Workflow

### Step-by-Step Debugging

```python
# 1. Check model structure
dt = DiagnosticsToolbox(m)
dt.report_structural_issues()

# 2. Verify DOF
from idaes.core.util.model_statistics import degrees_of_freedom
dof = degrees_of_freedom(m)
print(f"DOF = {dof}")
assert dof == 0, "Model not fully specified"

# 3. Check for variable issues
# - Unfixed variables that should be fixed
# - Fixed variables that should be unfixed
# - Variables with unreasonable values

# 4. Initialize carefully
# See initialization.md for strategies

# 5. Check numerical conditioning
dt.report_numerical_issues()

# 6. Apply scaling if needed
from idaes.core.util import scaling as iscale
iscale.calculate_scaling_factors(m)

# 7. Attempt solve with verbose output
from idaes.core.solvers import get_solver
solver = get_solver('ipopt')
results = solver.solve(m, tee=True)  # tee=True shows solver output

# 8. If solve fails, check residuals
dt.display_constraints_with_large_residuals(tol=1e-4)

# 9. Examine solver status
from pyomo.opt import TerminationCondition, SolverStatus
print(f"Termination: {results.solver.termination_condition}")
print(f"Status: {results.solver.status}")

# 10. Post-solve diagnostics
if results.solver.termination_condition != TerminationCondition.optimal:
    dt.report_numerical_issues()
```

## Model Statistics

### Comprehensive Statistics

```python
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_activated_constraints,
    number_activated_blocks,
    number_unused_variables
)

print("=" * 50)
print("Model Statistics")
print("=" * 50)
print(f"Degrees of Freedom: {degrees_of_freedom(m)}")
print(f"Variables: {number_variables(m)}")
print(f"Constraints: {number_activated_constraints(m)}")
print(f"Blocks: {number_activated_blocks(m)}")
print(f"Unused Variables: {number_unused_variables(m)}")
```

### Unit-by-Unit Statistics

```python
# Check each unit model
for unit_name in ['heater', 'flash', 'compressor']:
    if hasattr(m.fs, unit_name):
        unit = getattr(m.fs, unit_name)
        dof = degrees_of_freedom(unit)
        n_vars = number_variables(unit)
        n_cons = number_activated_constraints(unit)
        print(f"{unit_name}: DOF={dof}, Vars={n_vars}, Constraints={n_cons}")
```

## Best Practices

### 1. Check Early and Often

```python
# After each major model change
dt = DiagnosticsToolbox(m)
dt.report_structural_issues()

# Before solving
assert degrees_of_freedom(m) == 0
```

### 2. Use Verbose Solver Output

```python
# Always use tee=True when debugging
solver.solve(m, tee=True)

# Shows:
# - Iteration progress
# - Constraint violations
# - Optimization metrics
# - Error messages
```

### 3. Incremental Building

```python
# Test units individually before connecting
m.fs.heater.initialize()
solver.solve(m.fs.heater)  # Should work

# Then connect
m.fs.stream = Arc(source=m.fs.heater.outlet, destination=m.fs.flash.inlet)
expand_arcs(m)
```

### 4. Save Working States

```python
from idaes.core.util import to_json

# After successful solve
to_json(m, fname='working_model.json')

# Can reload if later changes break model
```

### 5. Systematic Scaling

```python
# Apply scaling before attempting solve
from idaes.core.util import scaling as iscale

iscale.calculate_scaling_factors(m.fs.heater)
iscale.calculate_scaling_factors(m.fs.flash)
# ... for all units

iscale.calculate_scaling_factors(m)
```

## Solver Output Interpretation

### IPOPT Messages

**Optimal Solution Found:**
```
EXIT: Optimal Solution Found.
```
✓ Model solved successfully

**Infeasible Problem:**
```
EXIT: Converged to a point of local infeasibility.
```
- Conflicting constraints
- Physically impossible specifications
- Need to relax constraints or fix different variables

**Maximum Iterations:**
```
EXIT: Maximum Number of Iterations Exceeded.
```
- Poor initialization
- Bad scaling
- Need better starting point or more iterations

**Restoration Failed:**
```
EXIT: Restoration Failed!
```
- Severe numerical difficulties
- Check for NaN/Inf values
- Review scaling and bounds

### Checking Results Object

```python
from pyomo.opt import TerminationCondition

if results.solver.termination_condition == TerminationCondition.optimal:
    print("Success!")
elif results.solver.termination_condition == TerminationCondition.infeasible:
    print("Infeasible problem - check constraints")
    dt.display_constraints_with_large_residuals()
elif results.solver.termination_condition == TerminationCondition.maxIterations:
    print("Hit iteration limit - try better initialization")
else:
    print(f"Solver terminated with: {results.solver.termination_condition}")
```

## Quick Reference

### Run Diagnostics

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
dt.report_numerical_issues()
dt.display_constraints_with_large_residuals(tol=1e-5)
```

### Check DOF

```python
from idaes.core.util.model_statistics import degrees_of_freedom
print(f"DOF: {degrees_of_freedom(m)}")
```

### Find Violated Constraints

```python
for c in m.component_data_objects(Constraint, active=True):
    if c.body is not None:
        residual = abs(value(c.body) - value(c.upper))
        if residual > 1e-5:
            print(f"{c.name}: {residual:.2e}")
```

### Check Scaling

```python
from idaes.core.util import scaling as iscale

badly_scaled = list(iscale.badly_scaled_var_generator(m))
print(f"Badly scaled: {len(badly_scaled)}")
```
