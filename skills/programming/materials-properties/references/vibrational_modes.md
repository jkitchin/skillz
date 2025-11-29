# Vibrational Modes and Frequencies

## Overview

Vibrational analysis determines normal modes and frequencies of atomic oscillations. Applications include:
- Thermodynamic properties (entropy, free energy, heat capacity)
- Infrared and Raman spectroscopy
- Transition state verification
- Zero-point energy corrections
- Isotope effects

## Theoretical Background

### Harmonic Approximation

Expand potential energy around minimum:

```
V(R) = V(R₀) + Σᵢ ∂V/∂Rᵢ|₀ ΔRᵢ + ½ Σᵢⱼ ∂²V/∂Rᵢ∂Rⱼ|₀ ΔRᵢΔRⱼ + ...
```

At minimum: ∂V/∂Rᵢ = 0, keep only quadratic term.

### Dynamical Matrix

The Hessian (force constant matrix):

```
Dᵢⱼ = ∂²V / ∂Rᵢ∂Rⱼ
```

Mass-weighted dynamical matrix:

```
D̃ᵢⱼ = Dᵢⱼ / √(mᵢ mⱼ)
```

### Normal Modes

Solve eigenvalue problem:

```
D̃ · e = ω² e
```

Where:
- `e`: Eigenvector (normal mode)
- `ω²`: Eigenvalue (squared frequency)

**Frequencies:**
```
ν = ω / (2π) = (1/2π) √(eigenvalue)
```

Units: commonly cm⁻¹ for spectroscopy

**Conversion:**
```
1 eV = 8065.5 cm⁻¹
```

### Zero-Point Energy

Quantum mechanical ground state energy:

```
ZPE = ½ Σᵢ ℏωᵢ
```

Can be significant correction (~0.1-0.5 eV for molecules).

## Computational Methods

### 1. Finite Differences

Numerically calculate Hessian:

**Central Differences:**
```
∂²V/∂Rᵢ∂Rⱼ ≈ [V(Rᵢ+δ, Rⱼ+δ) - V(Rᵢ+δ, Rⱼ-δ)
                - V(Rᵢ-δ, Rⱼ+δ) + V(Rᵢ-δ, Rⱼ-δ)] / (4δ²)
```

**Displacement:**
- Typically δ = 0.01 Å
- Smaller for tight potentials
- Larger for soft modes

**Cost:**
- 6N energy/force calculations for N atoms
- Can be parallelized

### 2. Using ASE

```python
from ase.vibrations import Vibrations

# System must be at minimum (tight convergence!)
atoms.calc = calculator

vib = Vibrations(atoms, indices=None, delta=0.01)
vib.run()
vib.summary()
vib.write_mode(mode_index)
```

**Parameters:**
- `indices`: Specify atoms to displace (None = all)
- `delta`: Displacement magnitude (Å)
- `nfree`: 2 (central diff, default) or 4 (more accurate)

### 3. Partial Mode Analysis

For large systems, only compute modes for subset:

```python
# Only adsorbate vibrations
adsorbate_indices = list(range(len(slab), len(slab_ads)))
vib = Vibrations(slab_ads, indices=adsorbate_indices)
vib.run()
```

**Note:** Coupling to substrate neglected if substrate fixed.

## Analysis

### 1. Mode Classification

**Molecular System:**
- 3N total degrees of freedom
- 3 translations (ω = 0)
- 3 rotations (ω = 0 for non-linear)
- 3N-6 vibrations (3N-5 for linear)

**Adsorbed Molecule:**
- Frustrated translations → low-frequency modes
- Frustrated rotations → librations
- Internal vibrations redshifted

**Imaginary Frequencies:**
- Negative eigenvalues → imaginary ω
- Indicates system not at minimum
- One imaginary mode → transition state
- Multiple imaginary → saddle point or poor geometry

### 2. Spectroscopic Properties

**IR Intensity:**
Requires Born effective charges:

```
I ∝ |∂μ/∂Q|²
```

Not directly from ASE Vibrations (need Berry phase calculations).

**Raman Intensity:**
Requires polarizability derivatives:

```
I ∝ |∂α/∂Q|²
```

More complex, typically needs specialized codes.

### 3. Thermodynamic Properties

**Partition Function:**
```
q_vib = Πᵢ 1 / (1 - exp(-ℏωᵢ/kT))
```

**Vibrational Contributions:**
- Free energy: F_vib = kT Σᵢ ln(1 - exp(-ℏωᵢ/kT))
- Entropy: S_vib = Σᵢ [ℏωᵢ/T × 1/(exp(ℏωᵢ/kT)-1) - k ln(1-exp(-ℏωᵢ/kT))]
- Heat capacity: C_v = Σᵢ k(ℏωᵢ/kT)² exp(ℏωᵢ/kT) / (exp(ℏωᵢ/kT)-1)²

```python
from ase.thermochemistry import HarmonicThermo

# Get frequencies in eV
vib_energies = vib.get_energies()

# Calculate thermodynamic properties
thermo = HarmonicThermo(vib_energies, potentialenergy=E)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325)
```

## Typical Frequency Ranges

### Molecular Vibrations (cm⁻¹):
- C-H stretch: 2850-3000
- O-H stretch: 3200-3600
- C=O stretch: 1700-1750
- C=C stretch: 1620-1680
- C-C stretch: 800-1200

### Metal-Adsorbate (cm⁻¹):
- M-CO: 400-500
- M-O: 400-600
- M-H: 1500-2100

### Surface Phonons:
- Rayleigh wave: 100-300 cm⁻¹
- Optical modes: 100-500 cm⁻¹

## Convergence and Accuracy

### 1. Geometry Optimization

**Critical:**
- Vibrations very sensitive to geometry
- `fmax < 0.01 eV/Å` minimum
- `fmax < 0.001 eV/Å` for accurate thermochemistry

### 2. Displacement Size

Test δ dependence:
```python
deltas = [0.005, 0.01, 0.015, 0.02]
for d in deltas:
    vib = Vibrations(atoms, delta=d)
    vib.run()
```

Default δ = 0.01 Å usually good.

### 3. DFT Accuracy

**Systematic Errors:**
- GGA: typically overestimates by 50-100 cm⁻¹
- Hybrid functionals more accurate
- Scaling factors often applied

**Scaling:**
```
ν_scaled = 0.96 × ν_calculated  (for PBE)
```

### 4. Numerical Noise

Clean forces required:
- Tight SCF convergence
- Dense k-point grids
- High energy cutoffs

## Advanced Topics

### 1. Anharmonic Effects

Beyond harmonic approximation:
- Vibrational Self-Consistent Field (VSCF)
- Vibrational Perturbation Theory (VPT2)
- Path integral molecular dynamics

Significant at high temperature or for soft modes.

### 2. Temperature-Dependent Frequencies

Implicit anharmonicity:
```
ω(T) = ω(0) + Δω_thermal(T)
```

Calculate via quasiharmonic approximation or AIMD.

### 3. Isotope Effects

Mass effects on frequencies:

```
ω_isotope / ω = √(m / m_isotope)
```

Example: D₂O vs H₂O
```
ν(O-D) / ν(O-H) = √(m_H / m_D) ≈ 0.71
```

### 4. Mode Projection

Project out translations and rotations for molecules:
```python
# ASE automatically handles this for finite systems
# For surfaces: no translational modes to remove
```

## Practical Tips

1. **Tight Optimization:**
   - Cannot overemphasize this!
   - Spurious imaginary modes from poor geometry
   - Use tighter fmax than typical

2. **Check for Imaginary Modes:**
   - Zero for minimum
   - One for transition state
   - Multiple → reoptimize

3. **Large Systems:**
   - Compute only relevant modes (adsorbate)
   - Use model systems for full analysis
   - Consider resonance Raman

4. **Thermochemistry:**
   - Include ZPE corrections
   - Important for reaction energies
   - Temperature-dependent properties

5. **Visualization:**
   - Use ASE GUI or other tools
   - Helps identify mode character
   - Check symmetry

## Example Workflow

```python
from ase.build import molecule
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo

# 1. Optimize geometry (tight!)
mol = molecule('CO')
mol.calc = calculator

opt = BFGS(mol, trajectory='opt.traj')
opt.run(fmax=0.001)  # Very tight!

E = mol.get_potential_energy()

# 2. Calculate vibrations
vib = Vibrations(mol, delta=0.01)
vib.run()
vib.summary()

# Get frequencies
vib_energies = vib.get_energies()

# Check for imaginary modes
for i, energy in enumerate(vib_energies):
    if energy.imag != 0:
        print(f"Warning: Mode {i} is imaginary!")

# 3. Thermochemistry
thermo = HarmonicThermo(vib_energies, potentialenergy=E)

# At standard conditions
T = 298.15  # K
p = 101325  # Pa

F = thermo.get_helmholtz_energy(T)
G = thermo.get_gibbs_energy(T, p)
S = thermo.get_entropy(T, p)

print(f"Helmholtz free energy: {F:.3f} eV")
print(f"Gibbs free energy: {G:.3f} eV")
print(f"Entropy: {S:.6f} eV/K")

# 4. Visualize mode
vib.write_mode(0)  # Write mode to trajectory
```

## Common Pitfalls

1. **Loose Geometry:**
   - Most common error
   - Leads to imaginary modes
   - Must have fmax < 0.01 eV/Å

2. **Numerical Noise:**
   - From loose DFT settings
   - Increase k-points, cutoff
   - Tighter SCF convergence

3. **Wrong System:**
   - Molecule in periodic box has translations
   - Use cell = None or large box
   - Check for rigid translations

4. **Large Displacements:**
   - δ too large violates harmonic approximation
   - Especially for soft modes
   - Test convergence

5. **Partial Hessian:**
   - Computing only subset of atoms
   - Valid only if coupling negligible
   - Check by computing full Hessian

## Transition State Verification

**Requirements:**
- Exactly ONE imaginary frequency
- Mode corresponds to reaction coordinate
- Endpoints connected by IRC

```python
vib.summary()
# Should see one negative frequency

# Visualize the mode
vib.write_mode(-1)  # Imaginary mode usually first
```

## References

1. **Fundamental Theory:**
   - Wilson, Decius & Cross, "Molecular Vibrations," McGraw-Hill (1955)
   - Califano, "Vibrational States," Wiley (1976)

2. **Computational Methods:**
   - Allouche & Ferro, "Ab initio geometry optimization and vibrational calculations," J. Phys. Chem. 100, 1820 (1996)

3. **Surface Vibrations:**
   - Ibach, "Physics of Surfaces and Interfaces," Springer (2006)
   - Einstein, "Planck's Theory of Radiation and the Theory of the Specific Heat," Ann. Phys. 22, 180 (1907)

4. **Thermochemistry:**
   - Cramer, "Essentials of Computational Chemistry," 2nd ed. Wiley (2004)
   - McQuarrie, "Statistical Mechanics," University Science Books (2000)

5. **ASE Implementation:**
   - ASE documentation: https://wiki.fysik.dtu.dk/ase/ase/vibrations/vibrations.html

## See Also

- `phonon_properties.md`: Bulk phonon calculations with phonopy
- `neb_barriers.md`: Transition state location
- `thermal_transport.md`: Thermal conductivity from phonons
- `examples/phonon_calculation.py`: Phonopy integration examples
