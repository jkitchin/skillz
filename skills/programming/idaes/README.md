# IDAES Skill

Comprehensive guidance for using IDAES (Institute for the Design of Advanced Energy Systems) for process systems engineering, flowsheet modeling, and optimization.

## Overview

This skill provides expert assistance with IDAES for:
- Process flowsheet development
- Chemical process modeling
- Energy systems simulation
- Process optimization
- Property modeling
- Model diagnostics and debugging

## Quick Start

### Installation

```bash
# Create conda environment
conda create -n idaes python=3.11
conda activate idaes

# Install IDAES
pip install idaes-pse

# Get solver binaries (IPOPT, etc.)
idaes get-extensions

# Verify installation
idaes --version
```

### First Example

See `examples/simple_heater.py` for a complete working example.

```python
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater
from idaes.core.solvers import get_solver

# Create model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# Add property package
m.fs.properties = iapws95.Iapws95ParameterBlock()

# Add unit model
m.fs.heater = Heater(property_package=m.fs.properties)

# Specify inputs
m.fs.heater.inlet.flow_mol[0].fix(100)
m.fs.heater.inlet.pressure[0].fix(101325)
m.fs.heater.inlet.enth_mol[0].fix(5000)
m.fs.heater.heat_duty[0].fix(10000)

# Initialize and solve
m.fs.heater.initialize()
solver = get_solver()
results = solver.solve(m)

# View results
m.fs.heater.report()
```

## Documentation Structure

### Reference Files

- **`core-concepts.md`** - IDAES architecture, flowsheet blocks, property packages, unit models
- **`flowsheets.md`** - Building and connecting process flowsheets
- **`initialization.md`** - Strategies for initializing complex models
- **`diagnostics.md`** - Debugging models, finding structural and numerical issues
- **`scaling.md`** - Applying proper scaling for numerical robustness
- **`property-packages.md`** - Working with thermodynamic properties
- **`unit-models.md`** - Using and configuring unit operation models
- **`generic-models.md`** - Generic model library (mixers, heaters, separators, etc.)
- **`power-generation.md`** - Power plant specific models
- **`gas-solid-models.md`** - Gas-solid contactor models
- **`solving.md`** - Solving flowsheets and handling convergence issues
- **`optimization.md`** - Process optimization and design studies
- **`parameter-estimation.md`** - Fitting models to data
- **`costing.md`** - Economic analysis and process costing
- **`dynamic-modeling.md`** - Time-dependent simulations
- **`custom-models.md`** - Developing custom unit and property models

### Examples

- **`simple_heater.py`** - Basic single-unit model
- **`flash_separation.py`** - Multi-unit flowsheet with separation
- More examples available at https://idaes.github.io/examples-pse/latest/

## Common Workflows

### Building a Flowsheet

1. Create model and flowsheet block
2. Add property packages
3. Add unit models
4. Connect units with Arcs
5. Expand arcs
6. Specify inputs (fix variables)
7. Check degrees of freedom (should be 0)
8. Initialize units sequentially
9. Solve
10. Analyze results

### Debugging a Model

1. Run structural diagnostics
2. Check degrees of freedom
3. Verify specifications
4. Check scaling
5. Initialize carefully
6. Use verbose solver output
7. Check constraint residuals
8. Apply fixes based on diagnostics

See `references/diagnostics.md` for detailed debugging procedures.

### Optimization

1. Build and solve base case
2. Define objective function
3. Set decision variable bounds
4. Configure optimization solver
5. Solve optimization problem
6. Perform sensitivity analysis

See `references/optimization.md` for optimization workflows.

## Key Concepts

### Degrees of Freedom

Your model must be "square" - same number of variables as equations:

```python
from idaes.core.util.model_statistics import degrees_of_freedom
dof = degrees_of_freedom(m)
# Should be 0 for fully specified model
```

### Initialization

Models must be initialized before solving:

```python
# Sequential initialization
m.fs.unit1.initialize()
propagate_state(arc=m.fs.stream1)
m.fs.unit2.initialize()

# Or automatic
from idaes.core.util.initializer import BlockTriangularizationInitializer
initializer = BlockTriangularizationInitializer()
initializer.initialize(m.fs)
```

### Scaling

Apply scaling for numerical robustness:

```python
from idaes.core.util import scaling as iscale

iscale.set_scaling_factor(m.fs.heater.inlet.pressure, 1e-5)
iscale.set_scaling_factor(m.fs.heater.inlet.flow_mol, 1e-2)
iscale.calculate_scaling_factors(m)
```

### Diagnostics

Use diagnostics to find and fix issues:

```python
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

dt = DiagnosticsToolbox(m)
dt.report_structural_issues()
dt.report_numerical_issues()
```

## Common Property Packages

- **`iapws95`** - High-accuracy water/steam properties
- **`ideal gas`** - Ideal gas mixtures
- **`cubic EOS`** - Cubic equations of state (SRK, PR)
- **`modular`** - Flexible framework for custom properties

## Common Unit Models

### Generic Models
- Mixer, Splitter, Separator
- Heater, HeatExchanger
- Pump, Compressor, Turbine
- Flash, Distillation
- CSTR, PFR, EquilibriumReactor

### Power Generation
- Boiler, FiresideHeater
- SteamTurbine
- HeatRecoverySteamGenerator (HRSG)
- FeedWaterHeater

### Gas-Solid
- FixedBed
- FluidizedBed
- MovingBed

## External Resources

- **Official Documentation:** https://idaes-pse.readthedocs.io/
- **Examples:** https://idaes.github.io/examples-pse/latest/
- **GitHub Repository:** https://github.com/IDAES/idaes-pse
- **Support:** idaes-support@idaes.org
- **Tutorials:** https://idaes-pse.readthedocs.io/en/stable/tutorials/

## Getting Help

When working with IDAES, common questions can be answered by:

1. **Check degrees of freedom** - Most issues stem from incorrect specifications
2. **Run diagnostics** - DiagnosticsToolbox finds structural and numerical problems
3. **Review initialization** - Poor initialization is a common cause of solver failures
4. **Apply scaling** - Numerical conditioning is critical for complex models
5. **Check solver output** - Use `tee=True` to see what the solver is doing

## Tips for Success

1. **Start simple** - Build and test units individually before connecting
2. **Check DOF frequently** - After every specification change
3. **Initialize sequentially** - Follow process flow direction
4. **Scale early** - Apply scaling before initialization
5. **Use diagnostics** - Run diagnostics when anything goes wrong
6. **Save working states** - Use `to_json()` to save successful configurations
7. **Read solver output** - The solver tells you what's wrong
8. **Test incrementally** - Make one change at a time

## License

IDAES is open source software distributed under BSD-3 license.
