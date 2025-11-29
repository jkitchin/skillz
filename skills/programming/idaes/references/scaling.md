# Scaling

Guide to applying proper scaling for numerical robustness in IDAES models.

## Why Scaling Matters

Numerical solvers work best when:
- All variables are of similar magnitude (ideally O(1))
- Jacobian entries are well-balanced
- Condition numbers are reasonable (<1e8)

Poor scaling causes:
- Solver convergence failures
- Slow convergence
- Numerical precision loss
- Incorrect solutions

## IDAES Scaling Framework

IDAES provides the `iscale` module for systematic scaling.

```python
from idaes.core.util import scaling as iscale
```

## Basic Scaling Workflow

### 1. Set Scaling Factors for Variables

```python
# Set scaling factor for a variable
# Scaling factor = 1 / nominal_value

# For a pressure around 100000 Pa
iscale.set_scaling_factor(m.fs.heater.inlet.pressure, 1e-5)

# For flow around 100 mol/s
iscale.set_scaling_factor(m.fs.heater.inlet.flow_mol, 1e-2)

# For temperature around 300 K
iscale.set_scaling_factor(m.fs.heater.inlet.temperature, 1e-2)
```

**Rule of thumb:** Scaling factor ≈ 1/nominal value
- If variable is typically 100, use scaling factor 0.01
- If variable is typically 1e6, use scaling factor 1e-6

### 2. Calculate Constraint Scaling

```python
# After setting variable scaling, calculate constraint scaling
iscale.calculate_scaling_factors(m)
```

This automatically propagates variable scaling to constraints.

### 3. Verify Scaling

```python
# Check for badly scaled variables
badly_scaled = list(iscale.badly_scaled_var_generator(m))
print(f"Badly scaled variables: {len(badly_scaled)}")

for var, scale_factor in badly_scaled[:10]:  # Show first 10
    print(f"{var.name}: value={var.value}, scale={scale_factor}")
```

## Automatic Scaling

### Property Package Scaling

Many property packages support automatic scaling:

```python
from idaes.models.properties import iapws95

# Property package with auto-scaling
m.fs.properties = iapws95.Iapws95ParameterBlock(
    temperature_scale=1e-2,  # T ~ 300 K
    pressure_scale=1e-5,     # P ~ 100000 Pa
    enthalpy_scale=1e-5      # h ~ 100000 J/mol
)
```

### Unit Model Scaling

```python
# Set scaling on unit model creation
m.fs.heater = Heater(property_package=m.fs.properties)

# Calculate scaling for the unit
iscale.calculate_scaling_factors(m.fs.heater)
```

## Manual Scaling Strategy

### For Custom Models

```python
def set_scaling(m):
    """Apply scaling factors to model variables."""

    # Flows (typical: 10-1000 mol/s)
    for port in [m.fs.unit.inlet, m.fs.unit.outlet]:
        iscale.set_scaling_factor(port.flow_mol, 1e-2)

    # Pressures (typical: 1e5 - 1e7 Pa)
    for port in [m.fs.unit.inlet, m.fs.unit.outlet]:
        iscale.set_scaling_factor(port.pressure, 1e-5)

    # Temperatures (typical: 200-600 K)
    for port in [m.fs.unit.inlet, m.fs.unit.outlet]:
        iscale.set_scaling_factor(port.temperature, 1e-2)

    # Enthalpies (typical: 1e4 - 1e6 J/mol)
    for port in [m.fs.unit.inlet, m.fs.unit.outlet]:
        iscale.set_scaling_factor(port.enth_mol, 1e-5)

    # Heat duties (typical: 1e3 - 1e8 W)
    iscale.set_scaling_factor(m.fs.heater.heat_duty, 1e-5)

    # Calculate constraint scaling
    iscale.calculate_scaling_factors(m)
```

### Time-Indexed Variables

```python
# For dynamic models
for t in m.fs.time:
    iscale.set_scaling_factor(m.fs.heater.inlet.flow_mol[t], 1e-2)
    iscale.set_scaling_factor(m.fs.heater.inlet.pressure[t], 1e-5)
```

## Constraint Scaling

### Automatic from Variables

Most constraints get scaling from their variables:

```python
# After setting variable scaling
iscale.calculate_scaling_factors(m)

# Constraints are automatically scaled based on:
# - Variables that appear in them
# - Their coefficients
```

### Manual Constraint Scaling

For custom constraints:

```python
# Define a constraint
@m.fs.Constraint(m.fs.time)
def custom_constraint(b, t):
    return b.var1[t] == 1000 * b.var2[t] + b.var3[t]

# Set scaling factor
iscale.constraint_scaling_transform(m.fs.custom_constraint, 1e-3)
```

## Scaling Workflow for Flowsheets

### Complete Example

```python
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater, Flash
from idaes.core.util import scaling as iscale
from pyomo.network import Arc, expand_arcs

# 1. Build model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# 2. Add property package with scaling hints
m.fs.properties = iapws95.Iapws95ParameterBlock()

# 3. Add units
m.fs.heater = Heater(property_package=m.fs.properties)
m.fs.flash = Flash(property_package=m.fs.properties)

# 4. Connect
m.fs.s01 = Arc(source=m.fs.heater.outlet, destination=m.fs.flash.inlet)
expand_arcs(m)

# 5. Specify inputs
m.fs.heater.inlet.flow_mol[0].fix(100)      # mol/s
m.fs.heater.inlet.pressure[0].fix(101325)   # Pa
m.fs.heater.inlet.enth_mol[0].fix(5000)     # J/mol
m.fs.heater.heat_duty[0].fix(100000)        # W

m.fs.flash.heat_duty[0].fix(0)
m.fs.flash.deltaP[0].fix(0)

# 6. Set scaling factors
# Heater
iscale.set_scaling_factor(m.fs.heater.inlet.flow_mol, 1e-2)
iscale.set_scaling_factor(m.fs.heater.inlet.pressure, 1e-5)
iscale.set_scaling_factor(m.fs.heater.inlet.enth_mol, 1e-4)
iscale.set_scaling_factor(m.fs.heater.heat_duty, 1e-5)

# Flash
iscale.set_scaling_factor(m.fs.flash.inlet.flow_mol, 1e-2)
iscale.set_scaling_factor(m.fs.flash.inlet.pressure, 1e-5)

# 7. Calculate all constraint scaling
iscale.calculate_scaling_factors(m)

# 8. Check scaling quality
badly_scaled = list(iscale.badly_scaled_var_generator(m))
print(f"Badly scaled variables: {len(badly_scaled)}")

# 9. Initialize and solve
m.fs.heater.initialize()
iscale.propagate_state(m.fs.s01)
m.fs.flash.initialize()

from idaes.core.solvers import get_solver
solver = get_solver()
results = solver.solve(m, tee=True)
```

## Checking Scaling Quality

### Scaled Variable Values

Ideally, scaled variables should be O(1):

```python
# Check scaled variable magnitudes
for v in m.component_data_objects(Var, descend_into=True):
    if not v.fixed and v.value is not None:
        sf = iscale.get_scaling_factor(v)
        if sf is not None:
            scaled_value = v.value * sf
            # Good: 0.01 < |scaled_value| < 100
            if abs(scaled_value) > 100 or abs(scaled_value) < 0.01:
                print(f"Poorly scaled: {v.name}")
                print(f"  Value: {v.value}")
                print(f"  Scale factor: {sf}")
                print(f"  Scaled value: {scaled_value}")
```

### Jacobian Condition Number

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)
dt.report_numerical_issues()

# Look for "Jacobian Condition Number"
# Good: < 1e8
# Acceptable: 1e8 - 1e12
# Poor: > 1e12
```

### Badly Scaled Variables

```python
# IDAES provides a generator for finding badly scaled variables
for var, sf in iscale.badly_scaled_var_generator(m, large=100, small=0.01):
    scaled_val = var.value * sf
    print(f"{var.name}: scaled value = {scaled_val}")
```

## Common Scaling Patterns

### Physical Quantities

Typical scales for common process variables:

```python
# Flows
# Molar: 1-1000 mol/s → scale = 1e-2 to 1e-3
# Mass: 1-1000 kg/s → scale = 1e-2 to 1e-3
iscale.set_scaling_factor(flow_mol, 1e-2)
iscale.set_scaling_factor(flow_mass, 1e-2)

# Pressure
# 1 bar = 100000 Pa → scale = 1e-5
# 10 bar = 1000000 Pa → scale = 1e-6
iscale.set_scaling_factor(pressure, 1e-5)

# Temperature
# 300 K → scale = 1/300 ≈ 3e-3
# For convenience, often use 1e-2
iscale.set_scaling_factor(temperature, 1e-2)

# Enthalpy (molar)
# Typical: 1e4 - 1e6 J/mol → scale = 1e-4 to 1e-6
iscale.set_scaling_factor(enth_mol, 1e-5)

# Heat duties
# kW to MW range: 1e3 - 1e6 W → scale = 1e-3 to 1e-6
iscale.set_scaling_factor(heat_duty, 1e-5)

# Mole fractions
# 0-1 range → scale = 1
iscale.set_scaling_factor(mole_frac, 1)

# Density
# Water: ~1000 kg/m³ → scale = 1e-3
# Gas: 1-10 kg/m³ → scale = 1e-1 to 1
iscale.set_scaling_factor(density, 1e-3)
```

### Component Indexing

For multicomponent systems:

```python
# Set scaling for each component
components = ['H2O', 'CO2', 'N2', 'O2']
for comp in components:
    iscale.set_scaling_factor(m.fs.state.mole_frac_comp[comp], 1)
    iscale.set_scaling_factor(m.fs.state.flow_mol_comp[comp], 1e-2)
```

### Time-Indexed Scaling

```python
# Dynamic models
for t in m.fs.time:
    iscale.set_scaling_factor(m.fs.heater.inlet.flow_mol[t], 1e-2)
    iscale.set_scaling_factor(m.fs.heater.inlet.pressure[t], 1e-5)
```

## Advanced Scaling

### Constraint-Specific Scaling

```python
# For a specific constraint that's problematic
constraint = m.fs.custom_energy_balance

# Set scaling factor directly
iscale.constraint_scaling_transform(constraint, 1e-5)
```

### Scaling State Blocks

```python
# State blocks contain many variables
# Use calculate_scaling_factors to handle them

# For inlet state
iscale.calculate_scaling_factors(m.fs.heater.control_volume.properties_in)

# For outlet state
iscale.calculate_scaling_factors(m.fs.heater.control_volume.properties_out)
```

### Hierarchical Scaling

```python
# For complex flowsheets, scale hierarchically

# 1. Scale property packages
iscale.calculate_scaling_factors(m.fs.properties)

# 2. Scale individual units
for unit in [m.fs.heater, m.fs.flash, m.fs.compressor]:
    iscale.calculate_scaling_factors(unit)

# 3. Scale entire flowsheet
iscale.calculate_scaling_factors(m.fs)

# 4. Scale full model
iscale.calculate_scaling_factors(m)
```

## Debugging Scaling Issues

### Identify Problem Variables

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)

# This will show badly scaled variables
dt.report_numerical_issues()
```

### Check Jacobian Entries

```python
# Variables with extreme Jacobian values
dt.display_variables_with_extreme_jacobian_values()
```

### Iterative Refinement

```python
# 1. Apply initial scaling
apply_scaling(m)

# 2. Check what's still bad
badly_scaled = list(iscale.badly_scaled_var_generator(m))

# 3. Fix worst offenders
for var, sf in badly_scaled[:5]:  # Top 5
    nominal_value = var.value if var.value is not None else 1.0
    new_sf = 1.0 / abs(nominal_value) if nominal_value != 0 else 1.0
    iscale.set_scaling_factor(var, new_sf)

# 4. Recalculate
iscale.calculate_scaling_factors(m)

# 5. Repeat as needed
```

## Best Practices

### 1. Scale Early

Apply scaling before initialization:

```python
# Good
build_model(m)
set_scaling(m)
m.fs.heater.initialize()

# Avoid
build_model(m)
m.fs.heater.initialize()  # May fail without scaling
set_scaling(m)  # Too late
```

### 2. Use Physical Insight

Base scaling on expected variable magnitudes:

```python
# If you know pressure will be ~10 bar (1e6 Pa)
iscale.set_scaling_factor(pressure, 1e-6)

# If flow will be ~50 mol/s
iscale.set_scaling_factor(flow_mol, 0.02)
```

### 3. Be Consistent

Use consistent scaling across similar variables:

```python
# All pressures
for var in [m.fs.unit1.pressure, m.fs.unit2.pressure, m.fs.unit3.pressure]:
    iscale.set_scaling_factor(var, 1e-5)
```

### 4. Document Scaling Assumptions

```python
def apply_model_scaling(m):
    """
    Apply scaling factors to model.

    Assumptions:
    - Pressure: ~1-10 bar (1e5 - 1e6 Pa)
    - Flow: ~10-100 mol/s
    - Temperature: ~300-600 K
    - Heat duties: ~10 kW - 1 MW
    """
    # Implementation
    pass
```

### 5. Verify After Solving

```python
# After solving, check if scaling was adequate
results = solver.solve(m)

if results.solver.termination_condition == TerminationCondition.optimal:
    # Check final conditioning
    dt = DiagnosticsToolbox(m)
    dt.report_numerical_issues()
```

## Common Issues

### Issue: Variables Still Badly Scaled After calculate_scaling_factors

```python
# Problem: Some variables don't have values yet
# Solution: Set values before calculating scaling

# Initialize or set reasonable guesses
m.fs.heater.outlet.temperature[0].value = 350  # K
m.fs.heater.outlet.pressure[0].value = 101325  # Pa

# Then calculate scaling
iscale.calculate_scaling_factors(m)
```

### Issue: Scaling Factors Not Applied

```python
# Check if scaling factor is set
sf = iscale.get_scaling_factor(m.fs.heater.inlet.pressure)
if sf is None:
    print("No scaling factor set")
    iscale.set_scaling_factor(m.fs.heater.inlet.pressure, 1e-5)
```

### Issue: Conflicting Scales

```python
# Problem: Different parts of model use different scales for same quantity
# Solution: Use consistent scaling

# Create scaling function
def set_pressure_scaling(var, nominal_pressure=1e5):
    iscale.set_scaling_factor(var, 1/nominal_pressure)

# Apply consistently
set_pressure_scaling(m.fs.unit1.pressure)
set_pressure_scaling(m.fs.unit2.pressure)
set_pressure_scaling(m.fs.unit3.pressure)
```

## Quick Reference

### Basic Scaling

```python
from idaes.core.util import scaling as iscale

# Set variable scaling
iscale.set_scaling_factor(variable, 1e-5)

# Calculate constraint scaling
iscale.calculate_scaling_factors(m)

# Check scaling
badly_scaled = list(iscale.badly_scaled_var_generator(m))
```

### Check Quality

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)
dt.report_numerical_issues()
```

### Typical Values

```python
# Pressure (Pa): scale = 1e-5 to 1e-6
# Temperature (K): scale = 1e-2 to 1e-3
# Flow (mol/s): scale = 1e-2 to 1e-3
# Enthalpy (J/mol): scale = 1e-4 to 1e-6
# Heat duty (W): scale = 1e-4 to 1e-6
```
