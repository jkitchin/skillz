# Initialization Strategies

Comprehensive guide to initializing IDAES models for successful solving.

## Overview

Initialization is the process of finding reasonable starting values for all variables in a model. Good initialization is critical for:
- Solving nonlinear equation systems
- Achieving convergence in complex flowsheets
- Reducing solve time
- Avoiding local optima

## Why Initialization Matters

IDAES models are systems of nonlinear algebraic equations (or DAEs for dynamic models). Solvers like IPOPT use iterative methods that require good starting points. Poor initialization can lead to:
- Solver failures
- Convergence to incorrect solutions
- Very slow convergence
- Numerical errors

## Basic Initialization

### Single Unit Initialization

Most IDAES unit models have an `initialize()` method:

```python
from idaes.models.unit_models import Heater

# Create and specify unit
m.fs.heater = Heater(property_package=m.fs.properties)

# Fix inputs
m.fs.heater.inlet.flow_mol[0].fix(100)
m.fs.heater.inlet.pressure[0].fix(101325)
m.fs.heater.inlet.enth_mol[0].fix(5000)
m.fs.heater.heat_duty[0].fix(10000)

# Initialize
m.fs.heater.initialize()
```

**What happens during initialize():**
1. Unfixed variables are given initial guesses
2. A simplified problem may be solved first
3. Property calculations are performed
4. Mass and energy balances are satisfied
5. The full problem is solved to tight tolerances

### Initialization with Options

```python
# Control solver behavior during initialization
m.fs.heater.initialize(
    optarg={'tol': 1e-6, 'max_iter': 50},
    outlvl=idaes.core.util.initialization.InitializationLogLevel.INFO
)
```

**Common options:**
- `optarg`: Dictionary of solver options
- `outlvl`: Logging verbosity (QUIET, INFO, DEBUG, etc.)
- `solver`: Specify solver to use (default is 'ipopt')

## Sequential Initialization

For connected flowsheets, initialize units in process order:

### Manual Sequential Initialization

```python
# Flowsheet: Feed → Heater → Flash → Separator

# 1. Initialize heater (inlet specified)
m.fs.heater.initialize()

# 2. Propagate outlet to flash inlet
from idaes.core.util import propagate_state
propagate_state(arc=m.fs.stream_01)

# 3. Initialize flash
m.fs.flash.initialize()

# 4. Propagate to separator
propagate_state(arc=m.fs.stream_02)

# 5. Initialize separator
m.fs.separator.initialize()
```

### Using propagate_state

The `propagate_state` utility copies values from one port to another:

```python
from idaes.core.util import propagate_state

# After heater is initialized
propagate_state(
    arc=m.fs.stream_01,  # The arc connecting units
    overwrite_fixed=False  # Don't overwrite fixed variables
)
```

## Advanced Initialization Strategies

### BlockTriangularizationInitializer

IDAES provides advanced initialization tools that automatically determine solve order:

```python
from idaes.core.util.initializer import BlockTriangularizationInitializer

# Create initializer
initializer = BlockTriangularizationInitializer()

# Initialize entire flowsheet
initializer.initialize(m.fs)
```

**Advantages:**
- Automatically determines initialization order
- Handles complex flowsheet topology
- Works with recycle loops
- More robust than manual initialization

**Configuration:**
```python
from idaes.core.util.initialization import InitializationLogLevel

initializer = BlockTriangularizationInitializer(
    constraint_tolerance=1e-6,
    output_level=InitializationLogLevel.INFO
)

initializer.initialize(m.fs)
```

### Hierarchical Initialization

For complex flowsheets, initialize subsystems separately:

```python
# Initialize subsections
initializer = BlockTriangularizationInitializer()

# Initialize feed section
initializer.initialize(m.fs.feed_section)

# Initialize reaction section
initializer.initialize(m.fs.reaction_section)

# Initialize separation section
initializer.initialize(m.fs.separation_section)

# Finally, initialize full flowsheet
initializer.initialize(m.fs)
```

## Handling Recycle Loops

Recycle loops are challenging because the inlet depends on the outlet.

### Tear Stream Approach

```python
# Flowsheet with recycle: Feed → Reactor → Separator → Product
#                                  ↑              ↓
#                                  └─── Recycle ──┘

# 1. Don't create the recycle arc initially
# Build flowsheet without recycle arc

# 2. Initialize forward path
m.fs.reactor.initialize()
propagate_state(m.fs.s01)
m.fs.separator.initialize()

# 3. Make initial guess for recycle
propagate_state(m.fs.separator.recycle_outlet, m.fs.mixer.recycle_inlet)

# 4. Initialize mixer with guessed recycle
m.fs.mixer.initialize()

# 5. Now create recycle arc
m.fs.recycle = Arc(
    source=m.fs.separator.recycle_outlet,
    destination=m.fs.mixer.recycle_inlet
)
expand_arcs(m)

# 6. Solve full flowsheet
solver = get_solver()
solver.solve(m)
```

### Fixed Point Iteration

```python
# Alternative: iterate to convergence
from idaes.core.util import copy_port_values

max_iter = 10
tolerance = 1e-4

for i in range(max_iter):
    # Save current recycle values
    recycle_vals_old = {
        'flow': m.fs.separator.recycle_outlet.flow_mol[0].value,
        'temp': m.fs.separator.recycle_outlet.temperature[0].value,
        'pres': m.fs.separator.recycle_outlet.pressure[0].value
    }

    # Solve flowsheet
    results = solver.solve(m)

    # Check convergence
    recycle_vals_new = {
        'flow': m.fs.separator.recycle_outlet.flow_mol[0].value,
        'temp': m.fs.separator.recycle_outlet.temperature[0].value,
        'pres': m.fs.separator.recycle_outlet.pressure[0].value
    }

    # Calculate change
    max_change = max(
        abs(recycle_vals_new[k] - recycle_vals_old[k]) / recycle_vals_old[k]
        for k in recycle_vals_old
    )

    print(f"Iteration {i}: max change = {max_change}")

    if max_change < tolerance:
        print("Converged!")
        break
```

## Initialization from Previous Solutions

### Saving Initialized States

```python
from idaes.core.util import to_json, from_json

# After successful initialization/solve
to_json(m, fname='initialized_model.json')

# Later, load saved state
from_json(m, fname='initialized_model.json')

# Can now solve without re-initializing
solver.solve(m)
```

### Using Similar Cases

```python
# Solve base case
m.fs.heater.heat_duty[0].fix(10000)
m.fs.heater.initialize()
solver.solve(m)

# Save state
base_state = to_json(m, return_dict=True)

# Solve similar case with different heat duty
m.fs.heater.heat_duty[0].fix(15000)
from_json(m, sd=base_state)  # Use base case as initialization
solver.solve(m)  # Should converge quickly
```

## Initialization Troubleshooting

### Diagnosis

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

# Create diagnostics
dt = DiagnosticsToolbox(m)

# Check before initialization
dt.report_structural_issues()

# If initialization fails, check numerical issues
try:
    m.fs.unit.initialize()
except:
    dt.report_numerical_issues()
```

### Common Issues and Solutions

#### Issue 1: Degrees of Freedom Not Zero

```python
from idaes.core.util.model_statistics import degrees_of_freedom

dof = degrees_of_freedom(m.fs.heater)
if dof != 0:
    print(f"DOF = {dof}, unit is not properly specified")

    # Check what's fixed
    for v in m.fs.heater.component_data_objects(Var):
        if v.fixed:
            print(f"Fixed: {v.name} = {v.value}")
```

**Solution:** Fix or unfix variables until DOF = 0

#### Issue 2: Property Package Errors

```python
# Check state variable bounds
m.fs.heater.inlet.display()

# Ensure values are physical
# For IAPWS: T > 273.15 K, P > 611 Pa
# For ideal gas: T > 0 K, P > 0 Pa
```

**Solution:** Verify specifications are within valid ranges

#### Issue 3: Poor Initial Guesses

```python
# Set better initial guesses before initialization
m.fs.heater.outlet.temperature[0].value = 350  # K (reasonable guess)
m.fs.heater.outlet.pressure[0].value = 101325  # Pa
m.fs.heater.outlet.flow_mol[0].value = 100  # mol/s

m.fs.heater.initialize()
```

**Solution:** Provide reasonable guesses for unfixed variables

#### Issue 4: Scaling Issues

```python
from idaes.core.util import scaling as iscale

# Calculate and set scaling factors before initialization
iscale.calculate_scaling_factors(m.fs.heater)

# Then initialize
m.fs.heater.initialize()
```

**Solution:** Apply proper scaling (see `scaling.md`)

### Relaxed Initialization

For difficult cases, use relaxed tolerances:

```python
# Initialize with loose tolerances first
m.fs.unit.initialize(optarg={'tol': 1e-3, 'constr_viol_tol': 1e-3})

# Then solve with tight tolerances
solver = get_solver('ipopt')
solver.options['tol'] = 1e-6
solver.solve(m)
```

## Best Practices

### 1. Initialize Units Before Connecting

```python
# Good: Initialize individually, then connect
m.fs.heater.initialize()
m.fs.flash.initialize()

# Connect
m.fs.stream = Arc(source=m.fs.heater.outlet, destination=m.fs.flash.inlet)
expand_arcs(m)

# Bad: Connect first, then try to initialize
# m.fs.stream = Arc(...)
# expand_arcs(m)
# m.fs.heater.initialize()  # May fail due to connection
```

### 2. Follow Process Flow

Initialize in the direction of material flow:
- Upstream units first
- Use propagate_state between units
- Downstream units last

### 3. Use Hierarchical Approach

For large flowsheets:
1. Initialize individual units
2. Initialize subsections
3. Initialize full flowsheet

### 4. Save Successful Initializations

```python
# After successful initialization
to_json(m, fname='good_init.json')

# Reuse for similar problems
from_json(m, fname='good_init.json')
```

### 5. Check DOF at Every Step

```python
from idaes.core.util.model_statistics import degrees_of_freedom

# Before initialization
assert degrees_of_freedom(m.fs.heater) == 0

# After adding more units
assert degrees_of_freedom(m) == 0
```

## Initialization Workflow

### Recommended Sequence

```python
# 1. Build flowsheet
m = build_flowsheet()

# 2. Add property packages
add_property_packages(m)

# 3. Add unit models
add_units(m)

# 4. DON'T expand arcs yet
# (or expand, then unfix outlets for initialization)

# 5. Specify inputs (fix inlet conditions, heat duties, etc.)
specify_inputs(m)

# 6. Check DOF for each unit
for unit in [m.fs.unit1, m.fs.unit2, m.fs.unit3]:
    assert degrees_of_freedom(unit) == 0

# 7. Initialize sequentially
initialize_sequential(m)

# 8. Expand arcs (if not done earlier)
expand_arcs(m)

# 9. Verify full flowsheet DOF
assert degrees_of_freedom(m) == 0

# 10. Solve
results = solver.solve(m)
```

## Quick Reference

### Basic Initialization

```python
# Single unit
m.fs.unit.initialize()

# With options
m.fs.unit.initialize(
    optarg={'tol': 1e-6},
    outlvl=idaes.core.util.initialization.InitializationLogLevel.INFO
)
```

### Sequential with propagate_state

```python
from idaes.core.util import propagate_state

m.fs.unit1.initialize()
propagate_state(arc=m.fs.stream1)
m.fs.unit2.initialize()
```

### Automatic Initialization

```python
from idaes.core.util.initializer import BlockTriangularizationInitializer

initializer = BlockTriangularizationInitializer()
initializer.initialize(m.fs)
```

### Save/Load State

```python
from idaes.core.util import to_json, from_json

# Save
to_json(m, fname='model_state.json')

# Load
from_json(m, fname='model_state.json')
```

### Diagnostics

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
dt.report_numerical_issues()
```
