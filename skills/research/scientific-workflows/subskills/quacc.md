# quacc Subskill

## Overview

quacc (Quantum Accelerators) is a **high-level platform for computational materials science and quantum chemistry workflows**. Choose quacc for:
- High-throughput materials screening
- Pre-built computational chemistry recipes
- Multi-backend workflow execution (Parsl, Dask, Prefect)
- ASE-based materials calculations

## Quick Start

```python
from ase.build import bulk
from quacc.recipes.emt.core import relax_job

# Simple structure relaxation
atoms = bulk("Cu")
result = relax_job(atoms)

print(result["atoms"])  # Relaxed structure
print(result["results"]["energy"])  # Final energy
```

## Pre-Built Recipes

### Relaxation

```python
from quacc import job
from ase.build import molecule

@job
def my_relax(atoms):
    from ase.optimize import BFGS
    atoms.calc = get_calculator()
    opt = BFGS(atoms)
    opt.run(fmax=0.05)
    return {"atoms": atoms, "energy": atoms.get_potential_energy()}
```

### High-Throughput Workflows

```python
from quacc import flow, job
from quacc.recipes.emt.core import relax_job, static_job

@flow
def screening_workflow(atoms_list):
    results = []
    for atoms in atoms_list:
        relaxed = relax_job(atoms)
        energy = static_job(relaxed["atoms"])
        results.append(energy)
    return results
```

## Backend Selection

```python
from quacc import QuaccConfig

# Use Parsl for HPC
QuaccConfig.workflow_engine = "parsl"

# Use Dask for local parallel
QuaccConfig.workflow_engine = "dask"

# Use Prefect for monitoring
QuaccConfig.workflow_engine = "prefect"
```

## When to Use quacc

✓ Materials science workflows
✓ Quantum chemistry calculations
✓ High-throughput screening
✓ Pre-built computational recipes
✓ ASE integration needed

## Limitations

❌ General-purpose workflows (use Prefect/Parsl)
❌ Simple caching (use joblib)
❌ Non-materials applications

## References

- Docs: https://quantum-accelerators.github.io/quacc/
- GitHub: https://github.com/Quantum-Accelerators/quacc
- Examples: `../examples/materials_workflow.py`
