# Thermal Transport Properties

## Overview

Calculate thermal conductivity from lattice vibrations (phonons) and electronic contributions.

## Lattice Thermal Conductivity

### Phonon Boltzmann Transport Equation

```
κ_lattice = (1/3) Σ_qλ C_qλ v_qλ² τ_qλ
```

Where:
- C: Mode heat capacity
- v: Group velocity
- τ: Relaxation time (scattering)

### Calculation with phono3py

```python
from phono3py import Phono3py

# Setup
phonon3 = Phono3py(atoms,
                   supercell_matrix=[[2,0,0],[0,2,0],[0,0,2]])

# Generate displacements
phonon3.generate_displacements()

# Calculate forces for 2nd and 3rd order
# (many force calculations required!)

# Produce force constants
phonon3.produce_fc2()
phonon3.produce_fc3()

# Calculate lattice thermal conductivity
phonon3.run_thermal_conductivity(
    temperatures=[300, 400, 500],
    mesh=[19, 19, 19]
)

# Get results
kappa = phonon3.thermal_conductivity.kappa
```

### Scattering Mechanisms

**Anharmonic (phonon-phonon):**
- Three-phonon: q₁ + q₂ ↔ q₃
- Requires 3rd order force constants

**Isotope:**
- Mass disorder scattering
- Can be included in phono3py

**Defect:**
- Point defects, dislocations
- Reduce conductivity

## Electronic Contribution

Wiedemann-Franz law:
```
κ_electron = L σ T
```

Where:
- L = π²k_B²/3e² (Lorenz number)
- σ: Electrical conductivity

Typically dominant in metals.

## Typical Values (W/m/K at 300K)

- Diamond: ~2000 (very high)
- Cu: ~400 (metals)
- Si: ~150 (semiconductors)
- Al₂O₃: ~30 (insulators)
- Polymers: ~0.1-1

## Simplified Estimates

### Debye-Callaway Model

For quick estimates:
```python
from ase.thermochemistry import crystal_thermal_conductivity_estimate

# Requires:
# - Debye temperature
# - Sound velocity
# - Volume per atom

kappa = crystal_thermal_conductivity_estimate(
    theta_D=θ_D,
    v_sound=v_m,
    volume=V_atom
)
```

### Slack Model

For high-T limit:
```
κ = A M V^(1/3) θ_D³ / (γ² n^(2/3) T)
```

Parameters:
- M: Average mass
- V: Volume
- θ_D: Debye temperature
- γ: Grüneisen parameter
- n: Number of atoms in primitive cell

## Applications

- Thermoelectric materials (want low κ)
- Thermal management (electronics)
- Heat exchangers
- Thermal barrier coatings

## Tools

- **phono3py:** Lattice thermal conductivity
- **BoltzTraP:** Electronic transport
- **almaBTE:** Phonon Boltzmann solver
- **ShengBTE:** Third-order force constants

## References

- Togo, Chaput & Tanaka, "Distributions of phonon lifetimes," Phys. Rev. B 91, 094306 (2015)
- Slack, "Nonmetallic crystals with high thermal conductivity," J. Phys. Chem. Solids 34, 321 (1973)

## See Also

- `phonon_properties.md`: Phonon calculations with phonopy
- `electronic_structure.md`: Electronic transport
