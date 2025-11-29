# CALPHAD Integration

## Overview

CALPHAD (CALculation of PHAse Diagrams) combines thermodynamic modeling with computational databases. Integration with DFT provides ab initio phase diagrams.

## pycalphad Basics

```python
from pycalphad import Database, equilibrium, variables as v

# Load thermodynamic database
db = Database('database.tdb')

# Define system
comps = ['CU', 'AU', 'VA']  # VA = vacancy
phases = ['FCC_A1', 'LIQUID']

# Calculate equilibrium
eq = equilibrium(db, comps, phases,
                 {v.X('AU'): (0, 1, 0.01),
                  v.T: 1000,
                  v.P: 101325})

# Plot phase diagram
from pycalphad.plot import eqplot
eqplot(eq)
```

## DFT → CALPHAD Workflow

### 1. Calculate Formation Energies

```python
# For Cu-Au system
compositions = [
    ('Cu', 1.0, 0.0),
    ('Cu3Au', 0.75, 0.25),
    ('CuAu', 0.5, 0.5),
    ('CuAu3', 0.25, 0.75),
    ('Au', 0.0, 1.0)
]

H_form = {}
for name, x_Cu, x_Au in compositions:
    structure = create_structure(name)
    E = calculate_dft_energy(structure)
    H_f = E - x_Cu*E_Cu - x_Au*E_Au
    H_form[name] = H_f
```

### 2. Fit Interaction Parameters

Redlich-Kister polynomial:
```
G_excess = x_Cu * x_Au * Σᵢ Lᵢ (x_Cu - x_Au)ⁱ
```

Fit L₀, L₁, L₂, ... to DFT formation energies.

### 3. Create TDB File

```
ELEMENT CU FCC_A1 63.546 !
ELEMENT AU FCC_A1 196.967 !

PHASE FCC_A1 %  2 1   1 !
CONSTITUENT FCC_A1 :CU,AU : VA :  !

PARAMETER G(FCC_A1,CU;0) 298.15 +GHSERCU; 6000 N !
PARAMETER G(FCC_A1,AU;0) 298.15 +GHSERAU; 6000 N !
PARAMETER L(FCC_A1,CU,AU;0) 298.15 -12000; 6000 N !
```

## Phonon Contributions

Include vibrational free energy:

```
G(T) = H_DFT + F_vib(T)
```

Calculate F_vib via phonopy at each composition.

## Applications

- Multi-component phase diagrams
- Solidification paths
- Precipitation predictions
- High-throughput screening

## Databases

- **SGTE:** Pure elements
- **NIMS:** Multi-component alloys
- **NIST:** Various systems

## Tools

- **pymatgen:** TDB generation from DFT
- **ESPEI:** Automatic parameterization
- **Thermo-Calc:** Commercial software

## References

- Lukas, Fries & Sundman, "Computational Thermodynamics," Cambridge (2007)
- Otis & Liu, "pycalphad," J. Open Res. Softw. 5, 1 (2017)

## See Also

- `cluster_expansion.md`: Configurational energetics
- `formation_energy.md`: Thermodynamic stability
