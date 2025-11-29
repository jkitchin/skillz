# pycalphad Quick Reference

## Installation

```bash
pip install pycalphad
# or
conda install -c conda-forge pycalphad
```

## Essential Imports

```python
from pycalphad import Database, equilibrium, calculate, binplot, variables as v
import numpy as np
import matplotlib.pyplot as plt
```

## Core Functions

### Load Database
```python
db = Database('database.tdb')
```

### Binary Phase Diagram
```python
binplot(db, comps=['AL', 'ZN', 'VA'],
        phases=['LIQUID', 'FCC_A1', 'HCP_A3'],
        conditions={'P': 101325, 'T': (300, 1000, 10), 'X(ZN)': (0, 1, 0.01)})
```

### Equilibrium Calculation
```python
result = equilibrium(db,
                    comps=['FE', 'C', 'VA'],
                    phases=['FCC_A1', 'BCC_A2'],
                    conditions={v.T: 1200, v.P: 101325, v.X('C'): 0.01})
```

### Property Calculation (No Minimization)
```python
result = calculate(db,
                  comps=['AL', 'CU', 'VA'],
                  phases=['FCC_A1'],
                  T=800, P=101325,
                  N=1, X={'AL': 0.5, 'CU': 0.5})
```

## State Variables

```python
v.T         # Temperature (K)
v.P         # Pressure (Pa)
v.X('C')    # Mole fraction of component C
v.MU('FE')  # Chemical potential of Fe
v.N         # Moles
```

## Result Variables (xarray Dataset)

```python
result.Phase     # Phase names
result.NP        # Phase fractions (moles)
result.X         # Overall compositions
result.Y         # Site fractions (sublattice)
result.GM        # Gibbs energy (J/mol-atom)
result.HM        # Enthalpy (J/mol-atom)
result.SM        # Entropy (J/mol-atom-K)
result.CPM       # Heat capacity (J/mol-atom-K)
result.MU        # Chemical potentials (J/mol)
result.ACR_*     # Activity (component *)
```

## Common Patterns

### T-X Property Map
```python
temps = np.linspace(500, 1500, 50)
x_vals = np.linspace(0, 1, 50)
T_grid, X_grid = np.meshgrid(temps, x_vals)

result = equilibrium(db, comps, phases,
                    {v.T: T_grid.flatten(),
                     v.P: 101325,
                     v.X('AL'): X_grid.flatten()})

# Reshape and plot
property_map = result.GM.values.reshape(T_grid.shape)
plt.contourf(X_grid, T_grid, property_map)
```

### Activity Calculation
```python
result = equilibrium(db, comps, phases, conditions)
activity_fe = result['ACR_FE'].values.squeeze()
mu_fe = result['MU_FE'].values.squeeze()
```

### Driving Force
```python
# Metastable state (parent phase only)
parent_eq = equilibrium(db, comps, ['FCC_A1'], conditions)
gm_parent = parent_eq.GM.values.squeeze()

# Stable state (all phases)
stable_eq = equilibrium(db, comps, all_phases, conditions)
gm_stable = stable_eq.GM.values.squeeze()

# Driving force = energy reduction
df = gm_parent - gm_stable
```

## Units

- Temperature: Kelvin (K)
- Pressure: Pascal (Pa), 101325 Pa = 1 atm
- Energy: J/mol or J/mol-atom
- Composition: Mole fraction (0-1)

## Essential Tips

1. **Always include 'VA' (vacancy)** for phases with vacant sites
2. **Check phase names**: `list(db.phases.keys())` (case-sensitive!)
3. **Specify n-1 compositions**: For n components, last composition is automatic
4. **Use pdens for accuracy**: Higher pdens = better phase boundaries (default 2000)
5. **Validate results**: Check `result.NP.sum() ≈ 1` and `result.NP >= 0`

## Common Errors

**"Component not found"**
- Check `db.elements` for exact spelling/case
- Remember to include 'VA' if needed

**"Phase not found"**
- Check `list(db.phases.keys())` for exact names
- Phase names are case-sensitive

**Composition doesn't sum to 1**
- Specify only n-1 independent fractions
- pycalphad calculates the dependent fraction

**Negative phase fractions**
- Calculation didn't converge
- Try adjusting pdens or initial conditions

**Slow calculations**
- Reduce pdens for faster (less accurate) results
- Reduce grid density for property maps
- Use fewer phases if some are unlikely

## xarray Dataset Operations

```python
# Select specific conditions
subset = result.sel(T=1200, method='nearest')

# Slice range
temps_slice = result.sel(T=slice(1000, 1500))

# Boolean filtering
liquid = result.where(result.Phase == 'LIQUID', drop=True)

# Convert to pandas
df = result.to_dataframe()

# Statistical operations
mean_gm = result.GM.mean()
std_cp = result.CPM.std()
```

## Debugging

```python
# Check what's in database
print("Elements:", db.elements)
print("Phases:", list(db.phases.keys()))

# Inspect phase definition
phase = db.phases['FCC_A1']
print("Constituents:", phase.constituents)

# Check for NaN in results
print("Failed points:", result.GM.isnull().sum().values)

# Visualize all phase energies
calc_all = calculate(db, comps, phases,
                     {v.T: 1000, v.P: 101325, v.X('AL'): (0, 1, 0.01)})
for phase in phases:
    phase_data = calc_all.where(calc_all.Phase == phase, drop=True)
    plt.plot(phase_data.X.sel(component='AL'), phase_data.GM, label=phase)
plt.legend()
plt.show()
```

## Resources

- **Docs**: https://pycalphad.org/docs/latest/
- **Examples**: https://pycalphad.org/docs/latest/examples/
- **GitHub**: https://github.com/pycalphad/pycalphad
- **Google Group**: pycalphad@googlegroups.com

## Citation

```
Otis, R. & Liu, Z.-K., (2017). pycalphad: CALPHAD-based
Computational Thermodynamics in Python. Journal of Open Research
Software. 5(1), p.1. DOI: http://doi.org/10.5334/jors.140
```
