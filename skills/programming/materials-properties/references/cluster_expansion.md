# Cluster Expansion

## Overview

Cluster expansion (CE) provides an efficient parameterization of configurational energetics in alloys and compounds. Maps DFT data to a rapid predictive model.

## Theory

### CE Hamiltonian

```
E(σ) = J₀ + Σᵢ Jᵢ ⟨σᵢ⟩ + Σᵢⱼ Jᵢⱼ ⟨σᵢσⱼ⟩ + Σᵢⱼₖ Jᵢⱼₖ ⟨σᵢσⱼσₖ⟩ + ...
```

Where:
- σᵢ: Occupation variable (+1 or -1 for binary)
- Jᵢ: Effective cluster interactions (ECIs)
- ⟨...⟩: Correlation functions

### Cluster Types

**Point:** Single site (composition)
**Pair:** Two sites (nearest-neighbor, next-nearest, etc.)
**Triplet:** Three sites (angles)
**Quadruplet:** Four sites

## Workflow with icet

```python
from icet import ClusterSpace, StructureContainer, ClusterExpansion
from icet.tools import enumerate_structures

# 1. Define cluster space
from ase.build import bulk

prim = bulk('Cu', 'fcc', a=3.6)
cutoffs = [6.0, 5.0, 4.0]  # pair, triplet, quadruplet

cs = ClusterSpace(prim, cutoffs, chemical_symbols=[['Cu', 'Au']])

# 2. Generate training structures
structures = enumerate_structures(
    prim,
    sizes=range(1, 9),
    chemical_symbols=[['Cu', 'Au']]
)

# 3. Calculate DFT energies
sc = StructureContainer(cs)

for structure in structures:
    structure.calc = calculator
    E = structure.get_potential_energy()

    # Formation energy
    E_form = calculate_formation_energy(structure)

    sc.add_structure(structure, properties={'energy': E_form})

# 4. Train
from trainstation import Optimizer

opt = Optimizer(sc.get_fit_data())
opt.train()

# 5. Create CE
ce = ClusterExpansion(cs, opt.parameters)

print(f"RMSE: {opt.rmse:.6f} eV/atom")
print(f"ECIs: {opt.parameters}")

# 6. Predict
new_structure = create_structure()
E_pred = ce.predict(new_structure)
```

## Monte Carlo Simulations

```python
from mchammer.ensembles import CanonicalEnsemble
from mchammer.calculators import ClusterExpansionCalculator

# Setup
supercell = make_supercell(prim, 4*np.eye(3))
calc = ClusterExpansionCalculator(supercell, ce)

# Run MC
mc = CanonicalEnsemble(supercell, calc, temperature=600)
mc.run(10000)

# Analyze
print(f"Final energy: {calc.calculate_total(supercell):.4f} eV")
```

## Ground State Search

```python
# Enumerate and evaluate all structures
results = []
for structure in enumerate_structures(...):
    E = ce.predict(structure)
    results.append((structure, E))

# Sort by energy
results.sort(key=lambda x: x[1])

# Ground states (convex hull)
ground_states = extract_hull(results)
```

## Applications

- Phase diagram prediction
- Ordering transitions
- Short-range order
- Precipitation in alloys

## Validation

**Cross-validation:**
- Hold out test set
- Predict vs DFT
- Check extrapolation

**Leave-one-out CV:**
```python
from sklearn.model_selection import LeaveOneOut

cv_scores = []
for train_idx, test_idx in LeaveOneOut().split(structures):
    # Train on train_idx, test on test_idx
    ce_cv = train_ce(structures[train_idx])
    error = predict(structures[test_idx])
    cv_scores.append(error)
```

## References

- Sanchez, Ducastelle & Gratias, "Generalized cluster description," Physica A 128, 334 (1984)
- Van de Walle et al., "Automating first-principles phase diagram calculations," J. Phase Equilibria 23, 348 (2002)
- icet documentation: https://icet.materialsmodeling.org

## See Also

- `examples/cluster_expansion.py`: Complete working examples
- `calphad_integration.md`: Thermodynamic modeling
