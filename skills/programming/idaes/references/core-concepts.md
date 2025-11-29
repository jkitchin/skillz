# IDAES Core Concepts

Understanding IDAES architecture and fundamental modeling concepts.

## Overview

IDAES is built on top of Pyomo, a Python-based optimization modeling language. Understanding the relationship between Pyomo and IDAES is essential for effective modeling.

### Pyomo Foundation

- **ConcreteModel**: Root container for all model components
- **Variables (Var)**: Decision variables and state variables
- **Constraints (Constraint)**: Equations and inequalities
- **Expressions**: Reusable algebraic expressions
- **Objectives**: Functions to minimize or maximize

### IDAES Extensions

IDAES adds process systems engineering specific components on top of Pyomo:
- FlowsheetBlock: Container for process models
- Unit Models: Equipment representations
- Property Packages: Thermodynamic calculations
- Control Volumes: Material and energy balance frameworks
- Ports and Arcs: Stream connectivity

## Fundamental Components

### 1. FlowsheetBlock

The main container for process models in IDAES.

```python
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock

# Create model
m = ConcreteModel()

# Add flowsheet
m.fs = FlowsheetBlock(dynamic=False)
```

**Key features:**
- Contains all unit models and property packages
- Manages time domain for dynamic simulations
- Provides configuration options
- Can be nested for hierarchical models

**Dynamic vs. Steady-State:**
```python
# Steady-state flowsheet
m.fs = FlowsheetBlock(dynamic=False)

# Dynamic flowsheet with time domain
m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 10, 20, 30])
```

### 2. Property Packages

Define thermodynamic and transport properties for chemical species.

```python
from idaes.models.properties import iapws95

# Create property parameter block
m.fs.properties = iapws95.Iapws95ParameterBlock()
```

**Components of a property package:**
- **Parameter Block**: Container for property parameters and methods
- **State Block**: Calculates properties at specific conditions
- **Component definitions**: Species in the system
- **Phase definitions**: Vapor, liquid, solid phases
- **Property calculations**: Enthalpy, entropy, density, etc.

**State Variables:**
Common sets of state variables:
- (flow, temperature, pressure) - FTP
- (flow, enthalpy, pressure) - FHP (used by IAPWS)
- (flow, mole fraction, temperature, pressure) - FcTP

```python
# Accessing state variables
m.fs.heater.inlet.flow_mol[0].fix(100)  # mol/s
m.fs.heater.inlet.pressure[0].fix(101325)  # Pa
m.fs.heater.inlet.enth_mol[0].fix(5000)  # J/mol
```

### 3. Unit Models

Represent process equipment with material and energy balances.

```python
from idaes.models.unit_models import Heater

# Create unit model
m.fs.heater = Heater(property_package=m.fs.properties)
```

**Unit model structure:**
- **Inlet/Outlet Ports**: Connection points for streams
- **Control Volume**: Material and energy balance equations
- **Configuration Options**: Customize behavior
- **Performance Equations**: Equipment-specific constraints
- **Variables**: Operating conditions and specifications

**Common configuration options:**
```python
from idaes.models.unit_models import Heater
from idaes.core import MaterialBalanceType, EnergyBalanceType

m.fs.heater = Heater(
    property_package=m.fs.properties,
    has_pressure_change=True,
    material_balance_type=MaterialBalanceType.componentTotal,
    energy_balance_type=EnergyBalanceType.enthalpyTotal
)
```

### 4. Control Volumes

The foundation for material and energy balances in unit models.

**Control Volume Types:**
- **ControlVolume0D**: Well-mixed, no spatial variation
- **ControlVolume1D**: Distributed along one dimension
- **ControlVolume0DBlock**: Multiple 0D control volumes

**Balance Types:**

Material Balance:
- `componentTotal`: Track each component separately
- `componentPhase`: Track each component in each phase
- `total`: Track total material only
- `elementTotal`: Track by elemental composition
- `useDefault`: Use package default

Energy Balance:
- `enthalpyTotal`: Total enthalpy balance
- `enthalpyPhase`: Separate balance per phase
- `energyTotal`: Total energy (includes kinetic/potential)
- `useDefault`: Use package default

### 5. Ports and Arcs

Connect unit models to represent process streams.

**Ports:**
```python
# Ports are automatically created by unit models
# Access port variables:
m.fs.heater.inlet.flow_mol
m.fs.heater.inlet.temperature
m.fs.heater.inlet.pressure
```

**Arcs:**
```python
from pyomo.network import Arc

# Create connection
m.fs.stream1 = Arc(
    source=m.fs.unit1.outlet,
    destination=m.fs.unit2.inlet
)

# Expand all arcs (required before solving)
from pyomo.network import expand_arcs
expand_arcs(m)
```

### 6. Time Domain

For dynamic simulations, IDAES uses a time domain.

```python
from pyomo.environ import Set

# Discrete time points
time_set = [0, 1, 2, 3, 4, 5]
m.fs = FlowsheetBlock(dynamic=True, time_set=time_set)

# Continuous time
from pyomo.dae import ContinuousSet
m.time = ContinuousSet(bounds=(0, 100))
m.fs = FlowsheetBlock(dynamic=True, time=m.time)

# Accessing time-indexed variables
m.fs.heater.inlet.flow_mol[0]  # at time = 0
m.fs.heater.inlet.flow_mol[5]  # at time = 5
```

## Model Hierarchy

IDAES models are organized hierarchically:

```
ConcreteModel
├── FlowsheetBlock (m.fs)
│   ├── Property Parameter Block (m.fs.props)
│   ├── Unit Model 1 (m.fs.unit1)
│   │   ├── Inlet Port (m.fs.unit1.inlet)
│   │   ├── Outlet Port (m.fs.unit1.outlet)
│   │   ├── Control Volume (m.fs.unit1.control_volume)
│   │   └── State Blocks
│   ├── Unit Model 2 (m.fs.unit2)
│   └── Arcs (m.fs.stream1, etc.)
```

## Degrees of Freedom

Critical concept for model specification.

**Definition:** Number of variables minus number of equations
- DOF = 0: Exactly specified (square system)
- DOF > 0: Under-specified (need more constraints)
- DOF < 0: Over-specified (redundant/conflicting constraints)

**Checking degrees of freedom:**
```python
from idaes.core.util.model_statistics import degrees_of_freedom

# Overall model
dof = degrees_of_freedom(m)
print(f"Degrees of freedom: {dof}")

# Specific unit
dof_heater = degrees_of_freedom(m.fs.heater)
print(f"Heater DOF: {dof_heater}")
```

**Typical specifications for a heater:**
```python
# 3 inlet conditions + heat duty = 4 specifications
m.fs.heater.inlet.flow_mol.fix(100)
m.fs.heater.inlet.pressure.fix(101325)
m.fs.heater.inlet.enth_mol.fix(5000)
m.fs.heater.heat_duty.fix(10000)

# This gives DOF = 0 for the heater
```

## Configuration Framework

IDAES uses a configuration system for customizing models.

```python
from idaes.core import useDefault

# View available configuration options
print(m.fs.heater.config)

# Set during construction
m.fs.heater = Heater(
    property_package=m.fs.properties,
    has_pressure_change=True,
    has_phase_equilibrium=False
)

# Some configurations can be changed after construction
m.fs.heater.config.has_pressure_change = True
```

**Common configuration options:**
- `dynamic`: Enable time-dependent behavior
- `has_holdup`: Include material holdup terms
- `has_pressure_change`: Allow pressure drop
- `has_phase_equilibrium`: Enable phase equilibrium calculations
- `has_heat_transfer`: Include heat transfer
- `has_work_transfer`: Include shaft work

## State Block Indexing

State blocks calculate properties at specific conditions.

```python
# Steady-state: indexed by time (single point)
m.fs.heater.control_volume.properties_in[0]  # Inlet at t=0
m.fs.heater.control_volume.properties_out[0]  # Outlet at t=0

# Dynamic: indexed by time
m.fs.heater.control_volume.properties_in[t]  # for each t in time set

# Access properties
T_in = m.fs.heater.control_volume.properties_in[0].temperature
h_out = m.fs.heater.control_volume.properties_out[0].enth_mol
```

## Naming Conventions

IDAES follows consistent naming:

**Variables:**
- `flow_mol`: Molar flow rate
- `flow_mass`: Mass flow rate
- `flow_vol`: Volumetric flow rate
- `mole_frac_comp`: Mole fraction by component
- `pressure`: Pressure
- `temperature`: Temperature
- `enth_mol`: Molar enthalpy
- `enth_mass`: Mass enthalpy

**Components:**
- Lowercase with underscores: `heat_exchanger`, `feed_stream`
- Descriptive names: `m.fs.boiler_feed_pump` not `m.fs.pump1`

## Units of Measurement

IDAES uses SI units by default:
- Length: meters (m)
- Time: seconds (s)
- Temperature: Kelvin (K)
- Pressure: Pascal (Pa)
- Energy: Joules (J)
- Power: Watts (W)
- Molar quantities: moles (mol)
- Mass: kilograms (kg)

**Unit conversion:**
```python
# Pressure conversion
P_atm = 101325  # Pa (1 atm)
P_bar = 1e5     # Pa (1 bar)
P_psi = 6894.76  # Pa (1 psi)

# Temperature conversion
T_C = 25  # Celsius
T_K = T_C + 273.15  # Kelvin

# Energy conversion
kW_to_W = 1000
MW_to_W = 1e6
```

## Best Practices

### Model Organization
- Use clear, descriptive names for units and streams
- Group related units in sub-flowsheets for complex models
- Document assumptions and specifications in comments

### Configuration Selection
- Choose appropriate balance types for your application
- Enable only necessary features (avoid unnecessary complexity)
- Use default configurations when possible

### Variable Specification
- Fix exactly the required number of variables (DOF = 0)
- Specify inlet conditions before outlet conditions
- Use reasonable initial guesses for unfixed variables

### Property Package Selection
- Choose property packages appropriate for your system conditions
- IAPWS95 for water/steam
- Ideal gas for low-pressure gas mixtures
- Modular property framework for custom systems

### Debugging
- Always check degrees of freedom before solving
- Use `.display()` to inspect variable values
- Check configuration with `.config`
- Verify connections with `.display_arcs()`

## Quick Reference

### Creating a Basic Model
```python
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.props = iapws95.Iapws95ParameterBlock()
m.fs.heater = Heater(property_package=m.fs.props)
```

### Fixing Variables
```python
m.fs.heater.inlet.flow_mol.fix(100)
m.fs.heater.inlet.pressure.fix(101325)
m.fs.heater.inlet.enth_mol.fix(5000)
m.fs.heater.heat_duty.fix(10000)
```

### Checking DOF
```python
from idaes.core.util.model_statistics import degrees_of_freedom
print(f"DOF: {degrees_of_freedom(m)}")
```

### Solving
```python
from idaes.core.solvers import get_solver
m.fs.heater.initialize()
solver = get_solver()
results = solver.solve(m)
```

### Displaying Results
```python
m.fs.heater.report()
m.fs.heater.outlet.display()
```
