# Flowsheet Construction

Building and connecting process flowsheets in IDAES.

## Overview

A flowsheet in IDAES represents a complete process model with unit operations connected by material and energy streams. This guide covers how to construct, connect, and manage flowsheets.

## Basic Flowsheet Creation

### Simple Flowsheet

```python
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock

# Create Pyomo model
m = ConcreteModel()

# Add IDAES flowsheet
m.fs = FlowsheetBlock(dynamic=False)
```

### Dynamic Flowsheet

```python
# With explicit time points
time_points = [0, 1, 2, 5, 10, 20, 30]
m.fs = FlowsheetBlock(dynamic=True, time_set=time_points)

# With continuous time using Pyomo DAE
from pyomo.dae import ContinuousSet
m.time = ContinuousSet(bounds=(0, 100))
m.fs = FlowsheetBlock(dynamic=True, time=m.time)
```

## Adding Unit Models

### Single Unit

```python
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater

# Add property package
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Add unit model
m.fs.heater = Heater(property_package=m.fs.properties)
```

### Multiple Units

```python
from idaes.models.unit_models import (
    Mixer, Heater, Flash, Compressor, Pump
)

# Add several units
m.fs.mixer = Mixer(
    property_package=m.fs.properties,
    num_inlets=2
)

m.fs.heat_exchanger = Heater(
    property_package=m.fs.properties,
    has_pressure_change=True
)

m.fs.flash = Flash(
    property_package=m.fs.properties,
    has_heat_transfer=True,
    has_pressure_change=True
)

m.fs.compressor = Compressor(
    property_package=m.fs.properties
)
```

## Connecting Units with Arcs

### Basic Connection

Arcs connect the outlet of one unit to the inlet of another.

```python
from pyomo.network import Arc

# Connect heater outlet to flash inlet
m.fs.s01 = Arc(
    source=m.fs.heater.outlet,
    destination=m.fs.flash.inlet
)

# Connect flash vapor outlet to compressor inlet
m.fs.s02 = Arc(
    source=m.fs.flash.vap_outlet,
    destination=m.fs.compressor.inlet
)
```

### Expanding Arcs

After defining all arcs, you must expand them to create the equality constraints:

```python
from pyomo.network import expand_arcs

# This creates the actual connection constraints
expand_arcs(m)
```

**Important:** Always call `expand_arcs()` before initializing or solving.

### Multiple Connections

```python
# Mixer with two inlets
m.fs.s01 = Arc(source=m.fs.unit1.outlet, destination=m.fs.mixer.inlet_1)
m.fs.s02 = Arc(source=m.fs.unit2.outlet, destination=m.fs.mixer.inlet_2)
m.fs.s03 = Arc(source=m.fs.mixer.outlet, destination=m.fs.unit3.inlet)

expand_arcs(m)
```

## Specifying Flowsheet Inputs

### Inlet Stream Specifications

Fix the required state variables at flowsheet inlets:

```python
# For FHP state variables (flow, enthalpy, pressure)
m.fs.heater.inlet.flow_mol[0].fix(100)  # mol/s
m.fs.heater.inlet.enth_mol[0].fix(5000)  # J/mol
m.fs.heater.inlet.pressure[0].fix(101325)  # Pa

# For FTP state variables (flow, temperature, pressure)
m.fs.heater.inlet.flow_mol[0].fix(100)  # mol/s
m.fs.heater.inlet.temperature[0].fix(300)  # K
m.fs.heater.inlet.pressure[0].fix(101325)  # Pa

# For multicomponent systems
m.fs.mixer.inlet_1.flow_mol[0].fix(50)
m.fs.mixer.inlet_1.temperature[0].fix(298.15)
m.fs.mixer.inlet_1.pressure[0].fix(101325)
m.fs.mixer.inlet_1.mole_frac_comp[0, "H2O"].fix(0.9)
m.fs.mixer.inlet_1.mole_frac_comp[0, "CO2"].fix(0.1)
```

### Unit Operation Specifications

Fix additional variables to fully specify the model:

```python
# Heater: specify heat duty
m.fs.heater.heat_duty[0].fix(10000)  # W

# Pump: specify pressure increase
m.fs.pump.deltaP[0].fix(500000)  # Pa

# Compressor: specify pressure ratio or outlet pressure
m.fs.compressor.ratioP[0].fix(3.0)  # dimensionless
# OR
m.fs.compressor.outlet.pressure[0].fix(303975)  # Pa

# Flash: specify temperature and pressure
m.fs.flash.heat_duty[0].fix(0)  # Adiabatic
m.fs.flash.deltaP[0].fix(0)  # No pressure drop
```

## Hierarchical Flowsheets

For complex processes, create sub-flowsheets:

```python
# Main flowsheet
m.fs = FlowsheetBlock(dynamic=False)

# Sub-flowsheet for separation section
m.fs.separation = FlowsheetBlock(dynamic=False)
m.fs.separation.props = m.fs.properties  # Share property package

# Add units to sub-flowsheet
m.fs.separation.flash1 = Flash(property_package=m.fs.properties)
m.fs.separation.flash2 = Flash(property_package=m.fs.properties)

# Connect within sub-flowsheet
m.fs.separation.s01 = Arc(
    source=m.fs.separation.flash1.liq_outlet,
    destination=m.fs.separation.flash2.inlet
)

# Connect to main flowsheet
m.fs.s_to_sep = Arc(
    source=m.fs.heater.outlet,
    destination=m.fs.separation.flash1.inlet
)

expand_arcs(m)
```

## Flowsheet Patterns

### Linear Chain

```python
# Unit1 → Unit2 → Unit3 → Unit4
m.fs.s01 = Arc(source=m.fs.unit1.outlet, destination=m.fs.unit2.inlet)
m.fs.s02 = Arc(source=m.fs.unit2.outlet, destination=m.fs.unit3.inlet)
m.fs.s03 = Arc(source=m.fs.unit3.outlet, destination=m.fs.unit4.inlet)
expand_arcs(m)
```

### Split and Merge

```python
from idaes.models.unit_models import Separator

# Split stream
m.fs.splitter = Separator(
    property_package=m.fs.properties,
    num_outlets=2,
    split_basis=SplittingType.totalFlow
)

m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.splitter.inlet)
m.fs.s02 = Arc(source=m.fs.splitter.outlet_1, destination=m.fs.unit1.inlet)
m.fs.s03 = Arc(source=m.fs.splitter.outlet_2, destination=m.fs.unit2.inlet)

# Specify split fraction
m.fs.splitter.split_fraction[0, "outlet_1"].fix(0.6)

# Merge streams
m.fs.mixer = Mixer(
    property_package=m.fs.properties,
    num_inlets=2
)

m.fs.s04 = Arc(source=m.fs.unit1.outlet, destination=m.fs.mixer.inlet_1)
m.fs.s05 = Arc(source=m.fs.unit2.outlet, destination=m.fs.mixer.inlet_2)

expand_arcs(m)
```

### Recycle Loop

```python
# Forward path: feed → reactor → separator
m.fs.s01 = Arc(source=m.fs.feed_mixer.outlet, destination=m.fs.reactor.inlet)
m.fs.s02 = Arc(source=m.fs.reactor.outlet, destination=m.fs.separator.inlet)

# Product stream
m.fs.s03 = Arc(source=m.fs.separator.outlet_1, destination=m.fs.product.inlet)

# Recycle stream (not an Arc for initialization purposes)
# Connect manually after initialization or use tear stream approach
```

**Note:** Recycle loops require special handling during initialization. See `initialization.md` for details.

## Checking Flowsheet Structure

### Degrees of Freedom

```python
from idaes.core.util.model_statistics import degrees_of_freedom

# Total flowsheet DOF
dof = degrees_of_freedom(m)
print(f"Flowsheet DOF: {dof}")

# Should be 0 for a fully specified flowsheet
assert dof == 0, f"Flowsheet is not fully specified (DOF = {dof})"

# Check individual units
for unit_name in ['heater', 'flash', 'compressor']:
    unit = getattr(m.fs, unit_name)
    unit_dof = degrees_of_freedom(unit)
    print(f"{unit_name} DOF: {unit_dof}")
```

### Model Statistics

```python
from idaes.core.util.model_statistics import (
    number_variables,
    number_activated_constraints,
    number_unused_variables
)

print(f"Variables: {number_variables(m)}")
print(f"Constraints: {number_activated_constraints(m)}")
print(f"Unused variables: {number_unused_variables(m)}")
```

### Stream Report

```python
# Display all streams in flowsheet
from idaes.core.util import stream_table

print(stream_table.generate_table(m))
```

## Visualization

### Display Connections

```python
# Show all arcs
m.fs.display_arcs()

# List all units
for component in m.fs.component_objects(descend_into=False):
    print(component.name)
```

### Export to Diagram Tools

```python
# Generate visualization data
from idaes.core.util.model_serializer import to_json
import json

# Export model structure
model_data = to_json(m, return_dict=True)
with open('flowsheet_structure.json', 'w') as f:
    json.dump(model_data, f, indent=2)
```

## Best Practices

### Naming Conventions

Use descriptive names:
```python
# Good
m.fs.boiler_feed_pump = Pump(...)
m.fs.main_steam_line = Arc(...)

# Avoid
m.fs.pump1 = Pump(...)
m.fs.s01 = Arc(...)  # OK for simple flowsheets, but consider descriptive names
```

### Property Package Management

Share property packages across units:
```python
# Create once
m.fs.water_props = iapws95.Iapws95ParameterBlock()

# Use in all water/steam units
m.fs.boiler = Heater(property_package=m.fs.water_props)
m.fs.turbine = Turbine(property_package=m.fs.water_props)
m.fs.condenser = Heater(property_package=m.fs.water_props)
```

### Modular Construction

Build flowsheets in functions for reusability:
```python
def build_heat_recovery_section(m):
    """Add heat recovery units to flowsheet."""
    m.fs.hrsg = HeatExchanger(property_package=m.fs.properties)
    m.fs.economizer = HeatExchanger(property_package=m.fs.properties)
    # ... connect units ...
    return m

def build_power_cycle(m):
    """Add steam turbine power cycle."""
    m.fs.hp_turbine = Turbine(property_package=m.fs.properties)
    m.fs.lp_turbine = Turbine(property_package=m.fs.properties)
    # ... connect units ...
    return m

# Use functions
m = build_flowsheet()
m = build_heat_recovery_section(m)
m = build_power_cycle(m)
```

### Documentation

Document specifications and assumptions:
```python
# Specify feed conditions
# Basis: 100 kmol/hr mixed feed at 25 C, 1 atm
# Composition: 80% N2, 20% O2 (molar basis)
m.fs.feed.flow_mol[0].fix(100/3600)  # Convert kmol/hr to mol/s
m.fs.feed.temperature[0].fix(298.15)  # K
m.fs.feed.pressure[0].fix(101325)  # Pa
m.fs.feed.mole_frac_comp[0, "N2"].fix(0.8)
m.fs.feed.mole_frac_comp[0, "O2"].fix(0.2)
```

### Systematic Build Process

1. **Create flowsheet block**
2. **Add property packages**
3. **Add unit models** (test individually if possible)
4. **Connect with arcs**
5. **Expand arcs**
6. **Specify inputs** (check DOF after each specification)
7. **Verify DOF = 0**
8. **Initialize**
9. **Solve**

## Common Flowsheet Examples

### Simple Process

```python
# Feed → Heater → Flash → Products
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater, Flash
from pyomo.network import Arc, expand_arcs

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.props = iapws95.Iapws95ParameterBlock()

m.fs.heater = Heater(property_package=m.fs.props)
m.fs.flash = Flash(property_package=m.fs.props)

m.fs.s01 = Arc(source=m.fs.heater.outlet, destination=m.fs.flash.inlet)
expand_arcs(m)
```

### Heat Exchanger Network

```python
# Multiple heat exchangers in series and parallel
m.fs.hx1 = HeatExchanger(property_package=m.fs.props)
m.fs.hx2 = HeatExchanger(property_package=m.fs.props)
m.fs.hx3 = HeatExchanger(property_package=m.fs.props)

# Hot side: hx1 → hx2 → hx3
m.fs.hot_s1 = Arc(source=m.fs.hx1.hot_outlet, destination=m.fs.hx2.hot_inlet)
m.fs.hot_s2 = Arc(source=m.fs.hx2.hot_outlet, destination=m.fs.hx3.hot_inlet)

# Cold side: hx3 → hx2 → hx1 (counter-current)
m.fs.cold_s1 = Arc(source=m.fs.hx3.cold_outlet, destination=m.fs.hx2.cold_inlet)
m.fs.cold_s2 = Arc(source=m.fs.hx2.cold_outlet, destination=m.fs.hx1.cold_inlet)

expand_arcs(m)
```

## Troubleshooting

### Arc Connection Errors

```python
# Check port names
print(m.fs.heater.inlet.vars)
print(m.fs.heater.outlet.vars)

# Verify units have same property package
assert m.fs.unit1.config.property_package is m.fs.unit2.config.property_package
```

### Missing expand_arcs

```python
# Symptom: DOF much higher than expected
# Solution: Call expand_arcs(m) before checking DOF

from pyomo.network import expand_arcs
expand_arcs(m)
```

### Over-specification

```python
# Find fixed variables
fixed_vars = []
for v in m.component_data_objects(Var, descend_into=True):
    if v.fixed:
        fixed_vars.append((v.name, v.value))

print(f"Fixed variables: {len(fixed_vars)}")
for name, value in fixed_vars:
    print(f"  {name} = {value}")
```

## Quick Reference

### Create and Connect

```python
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from pyomo.network import Arc, expand_arcs

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Add units (with property package)
m.fs.unit1 = UnitModel(property_package=m.fs.props)
m.fs.unit2 = UnitModel(property_package=m.fs.props)

# Connect
m.fs.stream = Arc(source=m.fs.unit1.outlet, destination=m.fs.unit2.inlet)
expand_arcs(m)
```

### Check Status

```python
from idaes.core.util.model_statistics import degrees_of_freedom

# DOF should be 0
print(f"DOF: {degrees_of_freedom(m)}")

# Stream table
from idaes.core.util import stream_table
print(stream_table.generate_table(m))
```
