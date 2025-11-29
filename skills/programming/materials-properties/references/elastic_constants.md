# Elastic Constants and Mechanical Properties

## Overview

Elastic constants characterize a material's response to applied stress/strain. Essential for:
- Mechanical stability
- Sound velocities
- Debye temperature
- Thermal expansion
- Phase transitions

## Theoretical Background

### Hooke's Law (Generalized)

```
σᵢⱼ = Σₖₗ Cᵢⱼₖₗ εₖₗ
```

Where:
- σ: Stress tensor (force/area)
- ε: Strain tensor (dimensionless)
- C: Elastic constant tensor (4th rank)

### Voigt Notation

Reduce to 6×6 matrix using symmetry:

```
C = | C₁₁ C₁₂ C₁₃ C₁₄ C₁₅ C₁₆ |
    | C₁₂ C₂₂ C₂₃ C₂₄ C₂₅ C₂₆ |
    | C₁₃ C₂₃ C₃₃ C₃₄ C₃₅ C₃₆ |
    | C₁₄ C₂₄ C₃₄ C₄₄ C₄₅ C₄₆ |
    | C₁₅ C₂₅ C₃₅ C₄₅ C₅₅ C₅₆ |
    | C₁₆ C₂₆ C₃₆ C₄₆ C₅₆ C₆₆ |
```

Indices: 1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy

### Crystal Symmetry

**Cubic (FCC, BCC, diamond):**
- 3 independent constants: C₁₁, C₁₂, C₄₄

**Hexagonal:**
- 5 independent: C₁₁, C₁₂, C₁₃, C₃₃, C₄₄

**Tetragonal:**
- 6 or 7 independent (depends on subgroup)

**Orthorhombic:**
- 9 independent

**Triclinic:**
- 21 independent (full tensor)

## Calculation Methods

### 1. Energy-Strain Method

Apply strain, calculate energy change:

```
E(ε) = E₀ + V₀ Σᵢⱼ σᵢⱼ εᵢⱼ + ½ V₀ Σᵢⱼₖₗ Cᵢⱼₖₗ εᵢⱼ εₖₗ
```

For small strain:
```
Cᵢⱼₖₗ = (1/V₀) ∂²E/∂εᵢⱼ∂εₖₗ
```

**Implementation:**
```python
from ase.constraints import ExpCellFilter

# Optimize to zero stress
atoms.calc = calculator
ecf = ExpCellFilter(atoms)
opt = BFGS(ecf)
opt.run(fmax=0.01)

# Apply strain
δ = 0.01  # Strain magnitude
E0 = atoms.get_potential_energy()
V0 = atoms.get_volume()

# Example: C₁₁ from εxx
atoms_strain = atoms.copy()
cell = atoms_strain.get_cell()
cell[0] *= (1 + δ)
atoms_strain.set_cell(cell, scale_atoms=True)
atoms_strain.calc = calculator

stress = atoms_strain.get_stress(voigt=True)
C11_approx = -stress[0] / δ
```

### 2. Stress-Strain Method

More accurate, use stress response:

```
σᵢ = Σⱼ Cᵢⱼ εⱼ
```

Fit stress vs strain to extract Cᵢⱼ.

### 3. DFPT (Density Functional Perturbation Theory)

Direct calculation without finite differences:
- More accurate
- No strain magnitude dependence
- Available in VASP (IBRION=6), Quantum ESPRESSO

### 4. Using elastic Package

Automates strain patterns:

```python
from elastic import get_elastic_tensor

atoms.calc = calculator
Cij, Bij = get_elastic_tensor(atoms, delta=0.01)
```

Returns both elastic (C) and compliance (B = C⁻¹) tensors.

## Derived Properties

### Bulk Modulus (B)

**Definition:**
```
B = -V (∂P/∂V)_T
```

**From Elastic Constants:**

Cubic:
```
B_Voigt = (C₁₁ + 2C₁₂) / 3
B_Reuss = 1 / (S₁₁ + 2S₁₂)
B_Hill = (B_Voigt + B_Reuss) / 2
```

### Shear Modulus (G)

**Voigt Average (upper bound):**
```
G_V = (C₁₁ - C₁₂ + 3C₄₄) / 5  (cubic)
```

**Reuss Average (lower bound):**
```
G_R = 5 / (4S₁₁ - 4S₁₂ + 3S₄₄)  (cubic)
```

**Hill Average:**
```
G = (G_V + G_R) / 2
```

### Young's Modulus (E)

```
E = 9BG / (3B + G)
```

### Poisson's Ratio (ν)

```
ν = (3B - 2G) / (2(3B + G))
```

Range: -1 to 0.5 (typically 0.2-0.4 for metals)

### Ductility Indicator

**Pugh's Ratio:**
```
B/G > 1.75  →  Ductile
B/G < 1.75  →  Brittle
```

### Anisotropy

**Zener Ratio (cubic):**
```
A = 2C₄₄ / (C₁₁ - C₁₂)
```

A = 1: Isotropic
A ≠ 1: Anisotropic

**Universal Anisotropy:**
```
A^U = 5 G_V/G_R + B_V/B_R - 6 ≥ 0
```

## Mechanical Stability

### Born Stability Criteria

**Cubic:**
```
C₁₁ - C₁₂ > 0
C₁₁ + 2C₁₂ > 0
C₄₄ > 0
```

**Hexagonal:**
```
C₁₁ > |C₁₂|
2C₁₃² < C₃₃(C₁₁ + C₁₂)
C₄₄ > 0
C₆₆ > 0
```

Violation → mechanically unstable (soft mode, phase transition)

## Sound Velocities

### From Elastic Constants

**Longitudinal:**
```
v_l = √(C₁₁ / ρ)
```

**Transverse:**
```
v_t = √(C₄₄ / ρ)
```

ρ: Mass density

**Average Sound Velocity:**
```
v_m = [(1/v_l³ + 2/v_t³) / 3]^(-1/3)
```

### Debye Temperature

```
θ_D = (h/k_B) v_m (3N / 4πV)^(1/3)
```

Important for:
- Thermal properties
- Superconductivity (T_c correlations)
- Phonon-limited transport

## Typical Values

### Bulk Modulus (GPa):
- Aluminum: 76
- Copper: 140
- Iron: 170
- Diamond: 440
- MgO: 160

### Young's Modulus (GPa):
- Aluminum: 70
- Copper: 130
- Iron: 210
- Diamond: 1050
- Steel: 200

### Poisson's Ratio:
- Most metals: 0.25-0.35
- Cork: ~0 (highly compressible)
- Rubber: ~0.5 (incompressible)

## Convergence and Accuracy

### Strain Magnitude

Test convergence:
```python
strains = [0.005, 0.01, 0.015, 0.02]
for δ in strains:
    C = calculate_elastic(delta=δ)
```

**Guidelines:**
- Too small: numerical noise
- Too large: anharmonic effects
- Sweet spot: 0.005-0.015 for most materials

### DFT Parameters

**k-points:**
- Denser than structural optimization
- Forces more sensitive than energies

**Energy Cutoff:**
- Converged for forces
- Stress tensor very sensitive

**Numerical Precision:**
- Use PREC=Accurate in VASP
- Tight SCF convergence

## Temperature Dependence

Elastic constants decrease with temperature:

```
C(T) = C₀ - α T
```

Calculate via:
1. Quasiharmonic approximation
2. Ab initio molecular dynamics
3. Thermodynamic integration

## Pressure Dependence

**Birch-Murnaghan EOS:**

```
B'(P) = ∂B/∂P
```

Typical B' ~ 4-6 for most materials.

## Practical Example

```python
from ase.build import bulk
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
import numpy as np

# Optimize structure
atoms = bulk('Cu', 'fcc', a=3.6)
atoms.calc = calculator

ecf = ExpCellFilter(atoms)
opt = BFGS(ecf)
opt.run(fmax=0.001)

# Apply strains for cubic system
δ = 0.01

# Isotropic strain for C₁₁ + 2C₁₂
atoms_iso = atoms.copy()
atoms_iso.set_cell(atoms.cell * (1 + δ), scale_atoms=True)
atoms_iso.calc = calculator
stress_iso = atoms_iso.get_stress(voigt=True)
C_sum = -stress_iso[0] / δ / 160.2177  # to GPa

# Volume-conserving for C₁₁ - C₁₂
atoms_vc = atoms.copy()
cell_vc = atoms.cell.copy()
cell_vc[0] *= (1 + δ)
cell_vc[1] *= (1 - δ)
atoms_vc.set_cell(cell_vc, scale_atoms=True)
atoms_vc.calc = calculator
stress_vc = atoms_vc.get_stress(voigt=True)
C_diff = -(stress_vc[0] - stress_vc[1]) / (2*δ) / 160.2177

# Shear for C₄₄
atoms_shear = atoms.copy()
cell_shear = atoms.cell.copy()
cell_shear[0, 1] += cell_shear[0, 0] * δ
atoms_shear.set_cell(cell_shear, scale_atoms=True)
atoms_shear.calc = calculator
stress_shear = atoms_shear.get_stress(voigt=True)
C44 = -stress_shear[5] / (2*δ) / 160.2177

# Extract C11, C12
C11 = (C_sum + 2*C_diff) / 3
C12 = (C_sum - C_diff) / 3

print(f"C₁₁ = {C11:.1f} GPa")
print(f"C₁₂ = {C12:.1f} GPa")
print(f"C₄₄ = {C44:.1f} GPa")

# Derived properties
B = (C11 + 2*C12) / 3
G = (C11 - C12 + 3*C44) / 5
E = 9*B*G / (3*B + G)
nu = (3*B - 2*G) / (2*(3*B + G))

print(f"\nBulk modulus: {B:.1f} GPa")
print(f"Shear modulus: {G:.1f} GPa")
print(f"Young's modulus: {E:.1f} GPa")
print(f"Poisson's ratio: {nu:.3f}")
```

## References

1. **Theory:**
   - Nye, "Physical Properties of Crystals," Oxford (1985)
   - Wallace, "Thermodynamics of Crystals," Wiley (1972)

2. **Computational:**
   - Nielsen & Martin, "First-principles calculation of stress," Phys. Rev. Lett. 50, 697 (1983)
   - Le Page & Saxe, "Symmetry-general least-squares extraction of elastic data," Phys. Rev. B 65, 104104 (2002)

3. **Databases:**
   - Elastic constants database: http://www.materialsproject.org
   - Landolt-Börnstein tables

## See Also

- `equation_of_state.md`: Bulk modulus from E-V curves
- `phonon_properties.md`: Connection to phonon frequencies
- `examples/elastic_constants.py`: Complete working examples
