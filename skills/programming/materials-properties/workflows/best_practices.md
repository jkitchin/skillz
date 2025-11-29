# Best Practices for Materials Properties Calculations

## Project Organization

### Directory Structure

```
project/
├── 00-bulk/
│   ├── optimization/
│   ├── convergence/
│   └── results/
├── 01-surface/
│   ├── 111-facet/
│   ├── 100-facet/
│   └── analysis/
├── 02-adsorption/
│   ├── CO/
│   ├── O2/
│   └── comparison/
├── 03-barriers/
│   ├── neb-calculations/
│   └── analysis/
├── scripts/
│   ├── setup.py
│   ├── analyze.py
│   └── plot.py
└── README.md
```

### File Naming

```
# Good
cu_fcc_relaxed.traj
co_adsorption_fcc_site.traj
neb_co_diffusion_001.traj

# Bad
final.traj
test2.traj
new_result_updated_v3.traj
```

## Code Quality

### Reproducible Scripts

```python
#!/usr/bin/env python3
"""
CO adsorption on Cu(111) surface.

Author: Your Name
Date: 2025-01-09
DFT Code: GPAW 23.9.1
Purpose: Calculate adsorption energies for CO at different sites
"""

from ase.build import fcc111, molecule, add_adsorbate
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT  # Replace with production calculator
import json

# Parameters (easy to modify)
METAL = 'Cu'
LATTICE_CONSTANT = 3.6  # Å
SLAB_SIZE = (3, 3, 4)
VACUUM = 10.0  # Å
FMAX = 0.02  # eV/Å
SITES = ['ontop', 'bridge', 'fcc', 'hcp']

# Calculator settings
def get_calculator():
    """Return configured calculator."""
    return EMT()  # Replace with production settings

# Main calculation
def main():
    """Calculate CO adsorption on Cu(111)."""

    print("="*60)
    print("CO Adsorption Study")
    print("="*60)

    # Clean slab
    slab = fcc111(METAL, size=SLAB_SIZE, vacuum=VACUUM, a=LATTICE_CONSTANT)

    z_pos = slab.get_positions()[:, 2]
    mask = [z < z_pos.min() + 3.0 for z in z_pos]
    slab.set_constraint(FixAtoms(mask=mask))

    slab.calc = get_calculator()
    opt = BFGS(slab, trajectory='slab_clean.traj')
    opt.run(fmax=FMAX)

    E_slab = slab.get_potential_energy()

    # CO molecule
    co = molecule('CO')
    co.center(vacuum=10.0)
    co.calc = get_calculator()
    opt_co = BFGS(co, trajectory='co_molecule.traj')
    opt_co.run(fmax=FMAX/2)

    E_co = co.get_potential_energy()

    # Adsorption
    results = {}

    for site in SITES:
        print(f"\nCalculating {site} site...")

        slab_ads = slab.copy()
        add_adsorbate(slab_ads, co, height=2.0, position=site)
        slab_ads.set_constraint(FixAtoms(mask=mask))

        slab_ads.calc = get_calculator()
        opt = BFGS(slab_ads, trajectory=f'co_{site}.traj')
        opt.run(fmax=FMAX)

        E_total = slab_ads.get_potential_energy()
        E_ads = E_total - E_slab - E_co

        results[site] = {
            'E_ads': E_ads,
            'E_total': E_total,
            'n_steps': opt.get_number_of_steps()
        }

        print(f"  E_ads = {E_ads:.4f} eV")

    # Save results
    with open('adsorption_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print("\nResults saved to adsorption_results.json")

if __name__ == "__main__":
    main()
```

### Configuration Files

```yaml
# config.yaml
system:
  element: Cu
  lattice_constant: 3.6
  crystal: fcc

slab:
  facet: 111
  size: [3, 3, 4]
  vacuum: 10.0
  fixed_layers: 2

calculator:
  name: GPAW
  mode: PW
  encut: 500
  kpts: [4, 4, 1]
  xc: PBE

convergence:
  fmax: 0.02
  energy: 1.0e-6
```

## Data Management

### Save All Data

```python
import pickle

# Save calculation details
metadata = {
    'calculator': str(atoms.calc),
    'parameters': atoms.calc.parameters,
    'energy': atoms.get_potential_energy(),
    'forces': atoms.get_forces().tolist(),
    'stress': atoms.get_stress().tolist(),
    'date': datetime.now().isoformat(),
    'hostname': socket.gethostname()
}

with open('metadata.pkl', 'wb') as f:
    pickle.dump(metadata, f)
```

### Trajectories

```python
from ase.io import write, read

# Save optimization trajectory
opt = BFGS(atoms, trajectory='optimization.traj')
opt.run(fmax=0.02)

# Analyze later
traj = read('optimization.traj', ':')
energies = [a.get_potential_energy() for a in traj]
```

## Error Handling

### Robust Calculations

```python
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    filename='calculation.log')

def safe_optimize(atoms, fmax=0.02, max_steps=200):
    """Optimize with error handling."""
    try:
        opt = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')
        opt.run(fmax=fmax, steps=max_steps)

        if opt.get_number_of_steps() >= max_steps:
            logging.warning("Optimization did not converge!")
            return False

        logging.info(f"Converged in {opt.get_number_of_steps()} steps")
        return True

    except Exception as e:
        logging.error(f"Optimization failed: {e}")
        return False
```

### Checkpointing

```python
import os

def calculate_with_restart(atoms, checkpoint='checkpoint.traj'):
    """Resume from checkpoint if it exists."""

    if os.path.exists(checkpoint):
        print(f"Resuming from {checkpoint}")
        atoms = read(checkpoint)

    opt = BFGS(atoms, trajectory=checkpoint)

    try:
        opt.run(fmax=0.02)
    except KeyboardInterrupt:
        print("\nInterrupted! Restart to resume.")
        write(checkpoint, atoms)
        raise

    return atoms
```

## Performance Optimization

### Parallel Calculations

```python
from multiprocessing import Pool

def calculate_structure(params):
    """Calculate single structure (for parallel execution)."""
    structure, site = params

    atoms = setup_calculation(structure, site)
    E = atoms.get_potential_energy()

    return site, E

# Parallel execution
sites = ['ontop', 'bridge', 'fcc', 'hcp']
params = [(base_structure, site) for site in sites]

with Pool(processes=4) as pool:
    results = pool.map(calculate_structure, params)
```

### Resource Management

```python
# Set calculator threads
calc.set(txt='output.txt',
         parallel={'domain': 1,
                   'band': 1,
                   'kpt': 4})  # 4 k-point parallelization
```

## Validation and Testing

### Unit Tests

```python
import unittest

class TestBulkProperties(unittest.TestCase):
    """Test bulk calculations."""

    def setUp(self):
        """Set up test structure."""
        self.atoms = bulk('Cu', 'fcc', a=3.6)
        self.atoms.calc = EMT()

    def test_lattice_constant(self):
        """Test optimized lattice constant."""
        opt = BFGS(self.atoms)
        opt.run(fmax=0.01)

        a = self.atoms.cell.cellpar()[0]
        self.assertAlmostEqual(a, 3.6, delta=0.1)

    def test_cohesive_energy(self):
        """Test cohesive energy is negative."""
        E = self.atoms.get_potential_energy() / len(self.atoms)
        self.assertLess(E, 0.0)

if __name__ == '__main__':
    unittest.main()
```

### Continuous Integration

```yaml
# .github/workflows/tests.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
    - uses: actions/checkout@v2

    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: 3.9

    - name: Install dependencies
      run: |
        pip install ase pytest

    - name: Run tests
      run: |
        pytest tests/
```

## Documentation

### Docstrings

```python
def calculate_surface_energy(slab, bulk_energy, fix_layers=2):
    """
    Calculate surface energy from slab model.

    Parameters
    ----------
    slab : ase.Atoms
        Slab structure with vacuum
    bulk_energy : float
        Energy per atom in bulk (eV)
    fix_layers : int
        Number of bottom layers to fix

    Returns
    -------
    gamma : float
        Surface energy in eV/Å²

    Examples
    --------
    >>> slab = fcc111('Cu', size=(3,3,7), vacuum=10)
    >>> gamma = calculate_surface_energy(slab, -3.5)
    >>> print(f"Surface energy: {gamma:.4f} eV/Å²")

    References
    ----------
    .. [1] Fiorentini & Methfessel, J. Phys.: Condens. Matter 8, 6525 (1996)
    """
    # Implementation
    pass
```

### README Files

````markdown
# CO Adsorption on Cu(111)

## Overview
This project calculates CO adsorption energies on Cu(111) surface.

## Requirements
- ASE >= 3.22
- GPAW >= 23.9
- Python >= 3.8

## Installation
```bash
pip install -r requirements.txt
```

## Usage
```bash
python calculate_adsorption.py
```

## Results
- `adsorption_results.json`: Adsorption energies
- `*.traj`: Optimization trajectories
- `plots/`: Visualization

## Parameters
See `config.yaml` for calculation settings.

## Citation
If you use this data, please cite:
- Author et al., Journal, Volume, Page (Year)
````

## Common Mistakes to Avoid

1. **Inconsistent Settings:** Use same DFT parameters for all related calculations
2. **Poor Geometry:** Always verify starting structures
3. **Missing Convergence Tests:** Never skip convergence studies
4. **No Version Control:** Use git for code and important data
5. **Inadequate Documentation:** Document assumptions and decisions
6. **Ignoring Errors:** Check for warnings in output files
7. **Not Saving Raw Data:** Keep trajectories and full outputs
8. **Hard-coded Parameters:** Use configuration files
9. **No Backups:** Regularly backup calculation data
10. **Publishing Without Validation:** Compare with experiments/benchmarks

## Publication Checklist

- [ ] All calculations converged properly
- [ ] Convergence tests documented
- [ ] Settings clearly specified
- [ ] Raw data archived
- [ ] Code available (GitHub/Zenodo)
- [ ] Supplementary information prepared
- [ ] Figures have proper labels and units
- [ ] Compared with experimental data when available
- [ ] Statistical uncertainties reported
- [ ] Computational cost documented

## Resources

- **ASE Documentation:** https://wiki.fysik.dtu.dk/ase/
- **GPAW Tutorials:** https://wiki.fysik.dtu.dk/gpaw/tutorials/tutorials.html
- **Materials Project:** https://materialsproject.org
- **Binder for reproducible notebooks:** https://mybinder.org

## See Also

- `convergence_testing.md`: Systematic convergence procedures
- `calculator_comparison.md`: Comparing different DFT codes
