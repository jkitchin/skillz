# Covalent Subskill

## Overview

Covalent orchestrates **quantum computing and cloud-native workflows** with infrastructure abstraction. Use for:
- Cloud-agnostic deployment (AWS, Azure, GCP)
- Quantum computing workflows
- ML/HPC across heterogeneous resources

## Quick Start

```python
import covalent as ct

@ct.electron
def task(x):
    return x * 2

@ct.lattice
def workflow(x):
    return task(x)

# Dispatch and retrieve
dispatch_id = ct.dispatch(workflow)(5)
result = ct.get_result(dispatch_id, wait=True)
print(result.result)  # 10
```

## Executors (Infrastructure)

```python
import covalent as ct

# Local
@ct.electron(executor="local")
def local_task(x):
    return x

# AWS
@ct.electron(executor="awsbatch")
def aws_task(x):
    return x

# SLURM
@ct.electron(executor="slurm")
def hpc_task(x):
    return x
```

## When to Use

✓ Multi-cloud deployment
✓ Quantum workflows
✓ Infrastructure independence
✓ Serverless HPC

## References

- GitHub: https://github.com/AgnostiqHQ/covalent
- Docs: https://docs.covalent.xyz/
