# NEB Reaction Barrier Calculations

## Overview

The Nudged Elastic Band (NEB) method finds minimum energy paths (MEPs) and transition states for chemical reactions and diffusion processes. Essential for:
- Catalytic reaction mechanisms
- Diffusion barriers
- Phase transitions
- Structural transformations

## Theoretical Background

### Minimum Energy Path (MEP)

The MEP is the path connecting reactants to products that has minimum energy perpendicular to the path at each point.

**Transition State:** Saddle point on MEP (maximum along reaction coordinate, minimum in all perpendicular directions)

### NEB Method

**Elastic Band:**
- Chain of images connecting initial and final states
- Images connected by spring forces
- Springs keep images evenly distributed

**Nudging:**
- Spring forces act only parallel to path
- True forces act only perpendicular to path
- Prevents "corner cutting"

**Energy:**
```
E_NEB = Σᵢ [V(Rᵢ) + k/2 |Rᵢ₊₁ - Rᵢ|²]
```

**Forces:**
```
F_i = F^perp_i + F^spring_i

F^perp = F^true - (F^true · τ̂) τ̂
F^spring = k [(|Rᵢ₊₁ - Rᵢ| - |Rᵢ - Rᵢ₋₁|)] τ̂
```

Where τ̂ is the local tangent to the path.

### Climbing Image NEB (CI-NEB)

**Improved Method:**
- Highest energy image "climbs" uphill
- Gives accurate transition state
- No spring force on climbing image

**Forces on Climbing Image:**
```
F_climb = F^true - 2(F^true · τ̂) τ̂
```

Maximizes energy along path, minimizes perpendicular.

## Computational Procedure

### 1. Endpoint Optimization

**Critical First Step:**
```python
# Initial state
initial.calc = calculator
opt_init = BFGS(initial)
opt_init.run(fmax=0.02)
E_init = initial.get_potential_energy()

# Final state
final.calc = calculator
opt_final = BFGS(final)
opt_final.run(fmax=0.02)
E_final = final.get_potential_energy()
```

**Requirements:**
- Tight convergence (fmax < 0.02 eV/Å)
- Same constraints on initial and final
- Verify correct final state (not local minimum)

### 2. Image Generation

**Linear Interpolation:**
```python
from ase.neb import NEB

n_images = 7  # Include endpoints
images = [initial]

for i in range(n_images - 2):
    image = initial.copy()
    image.calc = calculator
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
```

**Image Number Guidelines:**
- Simple reactions: 5-9 images
- Complex paths: 11-15 images
- Test convergence with respect to n_images

### 3. NEB Optimization

**Standard NEB:**
```python
from ase.optimize import BFGS

neb = NEB(images, climb=False)
neb.interpolate()

opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=0.05)
```

**Climbing Image NEB:**
```python
neb = NEB(images, climb=True)
neb.interpolate()

opt = BFGS(neb, trajectory='neb_climb.traj')
opt.run(fmax=0.05)
```

**Optimizer Choice:**
- BFGS: Generally reliable
- FIRE: Fast for rough paths
- MDMin: Good for difficult cases

### 4. Barrier Analysis

**Extract Barrier Heights:**
```python
from ase.neb import NEBTools

nebtools = NEBTools(images)

# Get barrier
Ef_forward, dE = nebtools.get_barrier()
Ef_reverse = Ef_forward - dE

print(f"Forward barrier: {Ef_forward:.3f} eV")
print(f"Reverse barrier: {Ef_reverse:.3f} eV")
print(f"Reaction energy: {dE:.3f} eV")
```

**Plot Energy Profile:**
```python
fig = nebtools.plot_band()
fig.savefig('neb_barrier.png')
```

## Convergence Considerations

### 1. Number of Images

Test convergence:
```python
n_images_list = [5, 7, 9, 11]
barriers = []

for n in n_images_list:
    neb = run_neb(n_images=n)
    Ef, _ = get_barrier(neb)
    barriers.append(Ef)
```

Converged when |Ef(n) - Ef(n-2)| < 0.05 eV

### 2. Spring Constant

Default k = 0.1 eV/Å² usually good.

**Adjust if:**
- Images bunch up: increase k
- Images too evenly spaced: decrease k

**Variable Spring Constants:**
```python
k = [0.1] * (n_images - 1)
k[peak_index-1] = 0.2
k[peak_index] = 0.2
neb = NEB(images, k=k, climb=True)
```

### 3. Parallel vs Sequential

**Sequential (default):**
- Optimize images one after another
- Good for initial convergence

**Parallel (faster):**
- Optimize all images simultaneously
- Requires parallel ASE/calculator
- Much faster for expensive calculators

### 4. Force Convergence

**Rough Optimization:**
- `fmax = 0.1 eV/Å`: Quick initial path
- Locate approximate transition state

**Final Optimization:**
- `fmax = 0.05 eV/Å`: Standard accuracy
- `fmax = 0.01 eV/Å`: Vibrational analysis

## Advanced Techniques

### 1. AutoNEB

Automatically adds images in regions with high curvature:
```python
from ase.neb import AutoNEB

autoneb = AutoNEB(calculator,
                  'initial.traj',
                  'final.traj',
                  n_max=15,
                  fmax=0.05)
autoneb.run()
```

### 2. String Method

Alternative to NEB:
- Reparameterizes path to equal arc length
- Can be more robust for certain systems
- Requires different implementation

### 3. Dimer Method

Find transition states without knowing final state:
```python
from ase.dimer import DimerControl, MinModeAtoms
```

Useful when final state is unknown.

### 4. Free Energy Barriers

Include entropic contributions:
- Metadynamics
- Umbrella sampling
- Thermodynamic integration

## Analysis Tools

### 1. Transition State Verification

**Vibrational Analysis:**
```python
from ase.vibrations import Vibrations

ts_image = images[peak_index]
vib = Vibrations(ts_image)
vib.run()
vib.summary()

# Should have exactly ONE imaginary frequency
```

**Intrinsic Reaction Coordinate (IRC):**
- Follow gradient from TS to verify endpoints
- Not directly in ASE (use special tools)

### 2. Reaction Coordinate

**Path Length:**
```python
path_length = nebtools.get_path_length()
```

**Collective Variables:**
- Bond distances
- Angles
- Coordination numbers

### 3. Energy Decomposition

Track contributions along path:
- Adsorbate-surface interaction
- Adsorbate internal energy
- Surface relaxation energy

## Typical Barriers

### Surface Diffusion:
- Metals on metals: 0.1-0.5 eV
- Adatoms on terraces: 0.2-1.0 eV
- Vacancy migration: 0.5-1.5 eV

### Chemical Reactions:
- H₂ dissociation: 0-1.5 eV (catalyst-dependent)
- CO oxidation: 0.5-1.5 eV
- N₂ dissociation: 1-3 eV (very challenging)

### Diffusion in Bulk:
- Vacancy migration: 0.5-2.0 eV
- Interstitial diffusion: 0.3-1.5 eV
- Ionic conductors: 0.3-1.0 eV

## Practical Tips

1. **Good Initial/Final States:**
   - Optimize tightly (fmax < 0.02 eV/Å)
   - Verify they're correct minima
   - Use same constraints

2. **Initial Path Guess:**
   - Linear interpolation usually sufficient
   - For complex reactions: hand-craft better guess
   - IDPP interpolation available in ASE

3. **Start with Rough NEB:**
   - Use fmax = 0.1 eV/Å initially
   - Locate approximate barrier
   - Refine with tighter convergence

4. **Monitor Convergence:**
   - Plot energy vs iteration
   - Check path evolution
   - Watch for image bunching

5. **Climbing Image:**
   - Use after initial convergence
   - Gives accurate TS geometry
   - Essential for vibrational analysis

## Common Pitfalls

1. **Poor Endpoints:**
   - Loose optimization → wrong barrier
   - Different constraints → artificial forces
   - Local minima → wrong reaction

2. **Too Few Images:**
   - Miss transition state
   - Inaccurate barrier
   - Always test convergence

3. **Wrong Final State:**
   - Multiple possible products
   - May need several NEB calculations
   - Check energetics carefully

4. **Artificial Constraints:**
   - Inconsistent constraints between endpoints
   - Over-constraining reaction coordinate
   - Let system find natural path

5. **Ignoring Multiple Paths:**
   - May be several mechanisms
   - Calculate all reasonable paths
   - Lowest barrier dominates kinetics

## Example Workflow

```python
from ase.build import fcc111, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.neb import NEB, NEBTools
from ase.constraints import FixAtoms

# 1. Create initial state
slab = fcc111('Cu', size=(4, 4, 4), vacuum=10)
mask = [atom.z < 6.0 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

initial = slab.copy()
# Add adsorbate at site A
add_adsorbate(initial, 'O', height=1.5, position=(0, 0))
initial.calc = calculator
opt = BFGS(initial)
opt.run(fmax=0.02)

# 2. Create final state
final = slab.copy()
# Add adsorbate at site B
add_adsorbate(final, 'O', height=1.5, position=(2.5, 0))
final.calc = calculator
opt = BFGS(final)
opt.run(fmax=0.02)

# 3. Run NEB
images = [initial]
for i in range(5):
    image = initial.copy()
    image.calc = calculator
    images.append(image)
images.append(final)

neb = NEB(images, climb=True)
neb.interpolate()

opt_neb = BFGS(neb, trajectory='neb.traj')
opt_neb.run(fmax=0.05)

# 4. Analyze
nebtools = NEBTools(images)
Ef, dE = nebtools.get_barrier()

print(f"Activation energy: {Ef:.3f} eV")
print(f"Reaction energy: {dE:.3f} eV")

fig = nebtools.plot_band()
fig.savefig('barrier.png')
```

## References

1. **Original NEB:**
   - Jónsson, Mills & Jacobsen, "Nudged elastic band method for finding minimum energy paths of transitions," in Classical and Quantum Dynamics in Condensed Phase Simulations (1998)

2. **Climbing Image:**
   - Henkelman, Uberuaga & Jónsson, "A climbing image nudged elastic band method for finding saddle points and minimum energy paths," J. Chem. Phys. 113, 9901 (2000)

3. **Improved NEB:**
   - Henkelman & Jónsson, "Improved tangent estimate in the nudged elastic band method," J. Chem. Phys. 113, 9978 (2000)

4. **AutoNEB:**
   - Kolsbjerg, Groves & Hammer, "An automated nudged elastic band method," J. Chem. Phys. 145, 094107 (2016)

5. **Reviews:**
   - Sheppard, Terrell & Henkelman, "Optimization methods for finding minimum energy paths," J. Chem. Phys. 128, 134106 (2008)

## See Also

- `adsorption_energy.md`: Initial/final states for surface reactions
- `vibrational_modes.md`: Transition state verification
- `examples/neb_barrier.py`: Working examples with CI-NEB
