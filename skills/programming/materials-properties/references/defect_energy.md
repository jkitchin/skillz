# Defect Energy Calculations

## Overview

Comprehensive guide for calculating formation energies of various defects:
- Point defects (vacancies, interstitials, substitutionals)
- Line defects (dislocations, Peierls barriers)
- Planar defects (stacking faults, grain boundaries)
- Solute-dislocation interactions

## 1. Vacancy Formation Energy

### Definition

```
E_f^v = E(N-1) - (N-1)/N × E(N) - μ
```

Where μ is the chemical potential (energy per atom in reservoir).

### Computational Procedure

```python
# Perfect supercell
supercell = make_supercell(bulk, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
supercell.calc = calculator
opt = BFGS(supercell)
opt.run(fmax=0.02)
E_perfect = supercell.get_potential_energy()
N = len(supercell)

# Create vacancy
vacancy_cell = supercell.copy()
del vacancy_cell[N//2]  # Remove central atom
vacancy_cell.calc = calculator
opt = BFGS(vacancy_cell)
opt.run(fmax=0.02)
E_vacancy = vacancy_cell.get_potential_energy()

# Formation energy
E_f = E_vacancy - ((N-1)/N) * E_perfect
```

### Convergence

- **Supercell size:** Typically 4×4×4 minimum
- **Charged defects:** Require special treatment (see below)
- **Relaxation volume:** ΔV = V_defect - V_perfect

**Typical Values (eV):**
- Cu: ~1.3
- Al: ~0.7
- Fe: ~2.0

## 2. Interstitial Defects

### Formation Energy

```
E_f^i = E(N+1) - (N+1)/N × E(N) + μ
```

### Interstitial Sites

**FCC:**
- Octahedral: (½, ½, 0) - usually most stable
- Tetrahedral: (¼, ¼, ¼) - higher energy

**BCC:**
- Octahedral: (½, 0, 0)
- Tetrahedral: (½, ¼, 0)

### Calculation

```python
# Add atom at octahedral site
from ase import Atom

atoms_int = supercell.copy()
cell = atoms_int.get_cell()
int_site = (cell[0] + cell[1]) / 2  # Octahedral

atoms_int.append(Atom('Cu', position=int_site))
atoms_int.calc = calculator

opt = BFGS(atoms_int)
opt.run(fmax=0.02)
E_int = atoms_int.get_potential_energy()

E_f = E_int - ((N+1)/N) * E_perfect
```

**Typical Values (eV):**
- Self-interstitials: 3-5 eV (high)
- Often metastable, important for radiation damage

## 3. Substitutional Defects

### Formation Energy

```
E_f = E(host+impurity) - E(host) - μ_impurity + μ_host
```

### Chemical Potentials

Reference to bulk reservoirs:
```python
# Host reference
host_bulk = bulk('Cu', 'fcc', a=3.6)
host_bulk.calc = calculator
mu_host = host_bulk.get_potential_energy() / len(host_bulk)

# Impurity reference
imp_bulk = bulk('Au', 'fcc', a=4.08)
imp_bulk.calc = calculator
mu_imp = imp_bulk.get_potential_energy() / len(imp_bulk)
```

### Calculation

```python
# Substitute central atom
subst_cell = supercell.copy()
subst_cell[N//2].symbol = 'Au'

subst_cell.calc = calculator
opt = BFGS(subst_cell)
opt.run(fmax=0.02)
E_subst = subst_cell.get_potential_energy()

E_f = E_subst - E_perfect - mu_imp + mu_host
```

### Solution Energy

Directly measures mixing tendency:
- E_f < 0: Favorable substitution
- E_f > 0: Segregation/precipitation expected

## 4. Stacking Fault Energy

### Intrinsic Stacking Fault (ISF)

FCC stacking: ...ABCABC...
ISF: ...ABCABABC... (remove layer C)

```
γ_sf = (E_fault - E_perfect) / A
```

### Calculation

```python
from ase.build import fcc111

# Perfect slab
slab = fcc111('Cu', size=(4, 4, 12), vacuum=0)
slab.calc = calculator
opt = BFGS(slab)
opt.run(fmax=0.02)
E_perfect = slab.get_potential_energy()

# Create fault by shifting upper half
positions = slab.get_positions()
z_mid = (positions[:, 2].max() + positions[:, 2].min()) / 2

# Shift by burgers vector component
b_sf = slab.cell[0] / 3  # For FCC: a/3[112]

for i, pos in enumerate(positions):
    if pos[2] > z_mid:
        slab.positions[i] += b_sf

slab.calc = calculator
opt = BFGS(slab)
opt.run(fmax=0.02)
E_fault = slab.get_potential_energy()

# Surface area
A = np.linalg.norm(np.cross(slab.cell[0], slab.cell[1]))
gamma_sf = (E_fault - E_perfect) / A

# Convert to mJ/m²
gamma_sf_SI = gamma_sf * 16.0217
```

**Typical Values (mJ/m²):**
- Cu: ~40
- Al: ~120
- Ag: ~20
- Stainless steel: ~50-100

### Extrinsic Stacking Fault (ESF)

FCC ESF: ...ABCABCBCABC... (extra layer)

Calculate similarly with different shift.

## 5. Peierls Barrier

Energy barrier for dislocation glide between equivalent lattice positions.

### Definition

```
E_Peierls = max(E(x)) - min(E(x))
```

Along glide direction x.

### Simplified Calculation

```python
# Create slab with dislocation
# (In practice: use atomman or matscipy for proper setup)

# Reference: perfect slab
slab_ref = create_slab()
E_ref = calculate_energy(slab_ref)

# Displaced configurations
n_images = 11
positions_fraction = np.linspace(0, 1, n_images)

energies = []
for frac in positions_fraction:
    slab_disp = apply_dislocation_displacement(slab_ref, frac)
    E = calculate_energy(slab_disp)
    energies.append(E - E_ref)

# Peierls barrier
E_P = max(energies) - min(energies)
```

**Key Considerations:**
- Requires large cells (>1000 atoms)
- Proper boundary conditions (flexible/fixed)
- Anisotropic elasticity

**Typical Values (meV/b):** (per unit Burgers vector length)
- Al: 10-50 meV
- Cu: 50-200 meV
- Fe (screw): ~1000 meV (high Peierls stress)

### Peierls Stress

```
τ_P = 2E_P / (b² √3)
```

Stress required to move dislocation.

## 6. Solute-Dislocation Interaction

### Interaction Energy

```
E_int = [E(disloc+solute) - E(disloc)] - [E(bulk+solute) - E(bulk)]
```

Measures binding of solute to dislocation.

### Physical Origin

**Size Misfit:**
- Oversized solute: attracted to tension regions
- Undersized solute: attracted to compression regions

**Modulus Misfit:**
- Soft solute: reduces strain energy
- Hard solute: increases local stiffness

**Electronic Effects:**
- Charge transfer
- d-band interactions
- Specific to material pair

### Calculation Procedure

**1. Bulk References:**
```python
# Pure bulk
bulk_pure = optimize_bulk('Al')
E_bulk = bulk_pure.get_potential_energy()

# Bulk with substitutional solute
bulk_solute = bulk_pure.copy()
bulk_solute[N//2].symbol = 'Mg'
E_bulk_solute = optimize(bulk_solute).get_potential_energy()

E_solution = E_bulk_solute - E_bulk
```

**2. Dislocation Model:**
```python
# Create dislocation (use specialized tools!)
from matscipy.dislocation import CubicCrystalDislocation

# Edge dislocation
disloc = CubicCrystalDislocation(
    'Al', alat, C11, C12, C44,
    crack_surface=[1,1,1],
    crack_direction=[1,-1,0],
    dislocation_type='edge'
)

# Or approximate with strained slab
slab_disloc = create_dislocation_slab()
E_disloc = optimize(slab_disloc).get_potential_energy()
```

**3. Solute Near Dislocation:**
```python
# Add solute at core
slab_disloc_solute = slab_disloc.copy()

# Find core position (max Von Mises stress)
core_idx = identify_dislocation_core(slab_disloc)
slab_disloc_solute[core_idx].symbol = 'Mg'

E_disloc_solute = optimize(slab_disloc_solute).get_potential_energy()
```

**4. Interaction Energy:**
```python
E_int = (E_disloc_solute - E_disloc) - E_solution
```

**Interpretation:**
- E_int < 0: Solute attracted to dislocation (solute atmosphere)
- E_int > 0: Solute repelled from dislocation
- |E_int| ~ 0.1-1 eV typical

### Binding Energy Profile

Calculate E_int vs distance from core:
```python
distances = [0, 5, 10, 15, 20]  # Angstroms
E_binding = []

for d in distances:
    site_idx = find_atom_at_distance(core_idx, d)
    E_int = calculate_interaction(site_idx)
    E_binding.append(E_int)

# Fit to continuum model
# E_bind(r) ∝ 1/r for edge dislocation
```

## Advanced Topics

### 1. Charged Defects

**Electrostatic Corrections:**
- Makov-Payne correction
- Freysoldt et al. (FNV) scheme
- Potential alignment

```
E_f^charged = E_raw + E_corr^Coulomb + E_corr^align
```

Implemented in pymatgen, pylada-defects.

### 2. Defect Complexes

Multiple defects interact:
- Vacancy clusters
- Impurity-vacancy pairs
- Divacancies

**Binding Energy:**
```
E_bind = E(complex) - Σᵢ E(isolated_i)
```

### 3. Migration Barriers

Use NEB between defect configurations:
```python
# Vacancy migration
initial = vacancy_at_site_A
final = vacancy_at_site_B

# Run NEB
barrier = calculate_neb(initial, final)
```

### 4. Thermodynamics

**Concentration:**
```
c = exp(-E_f / kT)
```

**Temperature Dependence:**
- Vibrational entropy contributions
- Include phonons at defect site

## Practical Tools

### Specialized Packages

**atomman (NIST):**
- Dislocation generation
- Elastic solutions
- Defect analysis

**matscipy:**
- Crack/dislocation tools
- Neighbor lists
- Elastic constants

**PyLammps:**
- Large-scale MD
- Defect dynamics
- Classical potentials

### Recommended Workflow

1. **Small system (DFT):**
   - Accurate energies
   - Electronic structure
   - ~100-500 atoms

2. **Large system (empirical potentials):**
   - Long-range strain fields
   - Realistic boundary conditions
   - ~10,000-100,000 atoms

3. **Hybrid QM/MM:**
   - DFT at core
   - Classical far field
   - Best of both

## References

1. **Point Defects:**
   - Freysoldt et al., "First-principles calculations for point defects in solids," Rev. Mod. Phys. 86, 253 (2014)

2. **Dislocations:**
   - Hirth & Lothe, "Theory of Dislocations," 2nd ed. Wiley (1982)
   - Olmsted et al., "Atomistic simulations of dislocation mobility in Al, Ni and Al/Mg alloys," Model. Simul. Mater. Sci. Eng. 13, 371 (2005)

3. **Stacking Faults:**
   - Vitek, "Intrinsic stacking faults in body-centred cubic crystals," Phil. Mag. 18, 773 (1968)

4. **Solute-Dislocation:**
   - Clouet et al., "Dislocation interaction with C in α-Fe," Phys. Rev. B 78, 064109 (2008)
   - Leyson, Curtin & Hector, "Solute strengthening from first principles," Phys. Rev. B 87, 174107 (2013)

5. **Computational Methods:**
   - Dederichs et al., "Lattice theory of point defects," J. Nucl. Mater. 69-70, 176 (1978)

## See Also

- `formation_energy.md`: General formation energy concepts
- `neb_barriers.md`: Migration barriers via NEB
- `examples/defect_formation.py`: Complete working examples for all defect types
