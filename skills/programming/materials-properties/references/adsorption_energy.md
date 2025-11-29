# Adsorption Energy Calculations

## Overview

Adsorption energy quantifies the strength of binding between an adsorbate (atom, molecule, or cluster) and a surface. Critical for:
- Heterogeneous catalysis
- Gas sensing
- Surface functionalization
- Corrosion and oxidation

## Theoretical Background

### Definition

The adsorption energy is defined as:

```
E_ads = E_slab+ads - E_slab - E_molecule
```

Where:
- `E_slab+ads`: Total energy of slab with adsorbate
- `E_slab`: Energy of clean slab
- `E_molecule`: Energy of isolated molecule/atom

**Sign Convention:**
- `E_ads < 0`: Exothermic, favorable adsorption
- `E_ads > 0`: Endothermic, unfavorable adsorption

### Differential vs Integral Adsorption Energy

**Differential (Local):**
```
E_ads(θ) = ∂E/∂N |_{θ}
```

Energy to add one adsorbate at coverage θ.

**Integral (Average):**
```
E_ads,avg(θ) = [E(θ) - E(0)] / N
```

Average energy per adsorbate at coverage θ.

## Adsorption Sites

### Common Sites on FCC(111):
1. **Top/Ontop:** Above a surface atom
2. **Bridge:** Between two surface atoms
3. **FCC hollow:** Above third-layer atom
4. **HCP hollow:** No atom directly below

### Typical Binding Strength Order:
```
Hollow (fcc/hcp) > Bridge > Top
```

But depends on adsorbate electronic structure!

## Computational Procedure

### 1. Preparation

**Clean Slab:**
```python
from ase.build import fcc111
from ase.constraints import FixAtoms

slab = fcc111('Pt', size=(3, 3, 4), vacuum=10.0)

# Fix bottom layers
z_pos = slab.get_positions()[:, 2]
mask = [z < z_pos.min() + 4.0 for z in z_pos]
slab.set_constraint(FixAtoms(mask=mask))

# Relax
slab.calc = calculator
opt.run(fmax=0.02)
E_slab = slab.get_potential_energy()
```

**Gas Phase Reference:**
```python
from ase.build import molecule

mol = molecule('CO')
mol.center(vacuum=10.0)
mol.calc = calculator
opt_mol.run(fmax=0.01)
E_mol = mol.get_potential_energy()
```

### 2. Adsorbate Placement

**Using add_adsorbate:**
```python
from ase.build import add_adsorbate

slab_ads = slab.copy()
add_adsorbate(slab_ads, mol, height=2.0, position='fcc')
```

**Manual Placement:**
```python
from ase import Atom

# Add atom at specific position
position = [x, y, z]
slab_ads.append(Atom('O', position=position))
```

### 3. Optimization Strategy

**Atomic Adsorbates:**
- Typically converge quickly
- Use BFGS or LBFGS

**Molecular Adsorbates:**
- May need careful initialization
- Consider different orientations
- Possible dissociation during optimization

**Convergence Criteria:**
- `fmax = 0.02-0.05 eV/Å` usually sufficient
- Tighter for vibrational analysis (`fmax < 0.01 eV/Å`)

### 4. Coverage Effects

**Low Coverage:**
- Use large supercells (> 3×3)
- Minimize adsorbate-adsorbate interactions

**High Coverage:**
- Smaller cells acceptable
- Model ordered adsorbate structures

**Coverage Definition:**
```
θ = N_ads / N_surface_atoms
```

Monolayer (ML): θ = 1

## Molecular vs Dissociative Adsorption

### Energy Comparison

**Molecular:**
```
E_ads(O2) = E_slab+O2 - E_slab - E_O2(gas)
```

**Dissociative:**
```
E_ads(2O) = E_slab+2O - E_slab - 2 × E_O(atom)
```

Or referenced to molecule:
```
E_diss = E_slab+2O - E_slab - E_O2(gas)
```

**Dissociation Barrier:**
Use NEB to find transition state for O₂ → 2O

## Convergence Considerations

### 1. Slab Size

Test lateral supercell size:
```python
sizes = [(2,2), (3,3), (4,4)]
for nx, ny in sizes:
    slab = fcc111('Cu', size=(nx, ny, 4))
    E_ads = calculate_adsorption()
```

Converged when |E_ads(n) - E_ads(n-1)| < 0.05 eV

### 2. Slab Thickness

Generally same as clean surface calculations (7-9 layers)

### 3. k-point Sampling

- Reduce k-points with supercell size
- For 3×3 surface: 4×4×1 often sufficient
- Always test convergence

### 4. Initial Height

Scan initial adsorbate heights:
```python
heights = np.linspace(1.5, 3.0, 7)
for h in heights:
    add_adsorbate(slab, mol, height=h)
    # optimize and compare
```

Avoid local minima!

## Analysis

### 1. Geometric Analysis

**Bond Lengths:**
```python
d = slab_ads.get_distance(surf_atom, ads_atom)
```

**Adsorption Height:**
```python
h = z_ads - max(z_surface_atoms)
```

**Molecular Deformation:**
- Compare geometry to gas phase
- Bond elongation indicates activation

### 2. Electronic Structure

**Charge Transfer:**
- Bader charge analysis
- Löwdin population analysis

**Density of States:**
- Project onto adsorbate orbitals
- Identify bonding/antibonding features

**Work Function Change:**
```
Δφ = φ(clean) - φ(adsorbate)
```

### 3. Vibrational Modes

Calculate vibrational frequencies:
```python
from ase.vibrations import Vibrations

vib = Vibrations(slab_ads, indices=[ads_indices])
vib.run()
vib.summary()
```

Redshifted modes indicate weakened bonds.

## Typical Values

### Atomic Adsorption (eV):
- O on Pt(111): -1.0 to -2.0
- H on Pd(111): -0.5 to -1.0
- N on Fe(100): -1.5 to -2.5

### Molecular Adsorption (eV):
- CO on Cu(111): -0.5 to -0.8
- NO on Pt(111): -1.5 to -2.0
- H2 on surfaces: typically < -0.3 (physisorption)

## Advanced Topics

### 1. Solvent Effects

Include implicit solvent models:
- VASPsol (VASP)
- iSMS (Quantum ESPRESSO)
- Hybrid implicit/explicit water models

### 2. Electric Field Effects

Apply external field to model electrode:
```python
from ase.constraints import ExternalForce
```

### 3. van der Waals Corrections

Essential for physisorption:
- DFT-D3 (Grimme)
- vdW-DF (non-local correlation)
- MBD (many-body dispersion)

### 4. Coverage-Dependent Energetics

**Lateral Interactions:**
```
E_ads(θ) = E_ads(0) + α θ + β θ²
```

Fit to multiple coverage calculations.

**Mean-Field Models:**
- Langmuir isotherm
- Frumkin isotherm (includes interactions)

## Practical Tips

1. **Multiple Orientations:**
   - Test different molecular orientations
   - Consider parallel vs perpendicular
   - May have multiple local minima

2. **Site Comparison:**
   - Systematically test all high-symmetry sites
   - Use similar optimization settings
   - Preferred site may be coverage-dependent

3. **Gas Phase Reference:**
   - Use large box (>10 Å vacuum)
   - Check for artificial constraints
   - Spin polarization for radicals

4. **Convergence:**
   - Tighter than bulk (fmax < 0.02 eV/Å)
   - Check forces on adsorbate AND surface atoms

5. **Initial Guess:**
   - Start with reasonable geometry
   - Use chemical intuition
   - Consider coordination preferences

## Common Pitfalls

1. **Small Supercells:**
   - Adsorbate-adsorbate interactions
   - Overestimated binding

2. **Wrong Gas Phase Reference:**
   - Use same DFT settings as surface
   - Check spin state (O₂ is triplet!)
   - Same XC functional

3. **Incomplete Optimization:**
   - Check all atoms relaxed properly
   - May need tighter fmax for molecules

4. **Neglecting Dispersion:**
   - Important for large molecules
   - Physisorption dominated by vdW

5. **Ignoring Coverage:**
   - Real systems often have finite coverage
   - Lateral interactions can be significant

## Example Workflow

```python
from ase.build import fcc111, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.constraints import FixAtoms

# 1. Clean slab
slab = fcc111('Cu', size=(3, 3, 4), vacuum=10.0)
z_pos = slab.get_positions()[:, 2]
mask = [z < z_pos.min() + 3.0 for z in z_pos]
slab.set_constraint(FixAtoms(mask=mask))
slab.calc = calculator

opt = BFGS(slab)
opt.run(fmax=0.02)
E_slab = slab.get_potential_energy()

# 2. Gas phase
mol = molecule('CO')
mol.center(vacuum=10.0)
mol.calc = calculator

opt_mol = BFGS(mol)
opt_mol.run(fmax=0.01)
E_mol = mol.get_potential_energy()

# 3. Adsorbed system
sites = ['ontop', 'bridge', 'fcc', 'hcp']
for site in sites:
    slab_ads = slab.copy()
    add_adsorbate(slab_ads, mol, height=2.0, position=site)
    slab_ads.set_constraint(FixAtoms(mask=mask))
    slab_ads.calc = calculator

    opt = BFGS(slab_ads)
    opt.run(fmax=0.02)

    E_total = slab_ads.get_potential_energy()
    E_ads = E_total - E_slab - E_mol

    print(f"{site}: E_ads = {E_ads:.3f} eV")
```

## References

1. **Fundamental Theory:**
   - Hammer & Nørskov, "Theoretical surface science and catalysis," Adv. Catal. 45, 71 (2000)
   - Nørskov et al., "Origin of the overpotential for oxygen reduction at a fuel-cell cathode," J. Phys. Chem. B 108, 17886 (2004)

2. **Computational Methods:**
   - Feibelman et al., "The CO/Pt(111) puzzle," J. Phys. Chem. B 105, 4018 (2001)
   - Wellendorff et al., "Density functionals for surface science," Phys. Rev. B 85, 235149 (2012)

3. **Coverage Effects:**
   - Reuter & Scheffler, "Composition, structure, and stability of RuO₂(110) as a function of oxygen pressure," Phys. Rev. B 65, 035406 (2001)

4. **Dispersion Corrections:**
   - Tkatchenko & Scheffler, "Accurate molecular van der Waals interactions from ground-state electron density," Phys. Rev. Lett. 102, 073005 (2009)

## See Also

- `surface_energy.md`: Surface stability
- `neb_barriers.md`: Reaction pathways on surfaces
- `examples/adsorption_energy.py`: Working examples with site comparison
