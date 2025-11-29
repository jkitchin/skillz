# Parsl Subskill

## Overview

Parsl (Parallel Scripting Library) enables **HPC scientific workflows** with implicit dataflow. Choose Parsl for:
- Running on HPC clusters (SLURM, PBS, SGE)
- Implicit task parallelism with data dependencies
- "Write once, run anywhere" across computing resources
- Interactive parallel computing in Jupyter

## Quick Start

```python
import parsl
from parsl.app.app import python_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor

# Configure for local execution
config = Config(
    executors=[HighThroughputExecutor(max_workers=4)]
)
parsl.load(config)

@python_app
def compute(x):
    return x ** 2

# Execute (returns AppFuture)
future = compute(10)
result = future.result()  # Get result (blocking)
```

## HPC Configuration

### SLURM Cluster

```python
from parsl.providers import SlurmProvider
from parsl.executors import HighThroughputExecutor
from parsl.config import Config

config = Config(
    executors=[
        HighThroughputExecutor(
            label="slurm_exec",
            provider=SlurmProvider(
                partition='compute',
                nodes_per_block=2,
                cores_per_node=28,
                max_blocks=10,
                walltime='02:00:00',
                scheduler_options='#SBATCH --account=myproject'
            )
        )
    ]
)
parsl.load(config)
```

### Multi-Site Execution

```python
config = Config(
    executors=[
        HighThroughputExecutor(
            label="local",
            max_workers=4
        ),
        HighThroughputExecutor(
            label="cluster",
            provider=SlurmProvider(...)
        )
    ]
)

@python_app(executors=['cluster'])
def heavy_compute(x):
    # Runs on cluster
    pass

@python_app(executors=['local'])
def light_process(x):
    # Runs locally
    pass
```

## Dataflow Programming

```python
@python_app
def stage1(x):
    return x * 2

@python_app
def stage2(x):
    return x + 1

@python_app
def stage3(x, y):
    return x + y

# Dataflow graph created automatically
a = stage1(10)  # Returns future
b = stage2(5)   # Returns future
c = stage3(a, b)  # Waits for a and b

result = c.result()  # Get final result
```

## Bash Apps

```python
from parsl.app.app import bash_app

@bash_app
def run_simulation(inputs=[], outputs=[], stdout='sim.out'):
    return f"./simulation.exe {inputs[0]} > {outputs[0]}"

# Execute
future = run_simulation(
    inputs=['params.txt'],
    outputs=['results.dat']
)
future.result()  # Wait for completion
```

## Common Patterns

### Parameter Sweep

```python
@python_app
def simulate(param):
    import time
    time.sleep(1)  # Expensive simulation
    return param ** 2

# Launch all tasks (non-blocking)
futures = [simulate(i) for i in range(100)]

# Collect results
results = [f.result() for f in futures]
```

### File-Based Workflow

```python
from parsl import File

@bash_app
def preprocess(inputs=[], outputs=[]):
    return f"python preprocess.py {inputs[0]} {outputs[0]}"

@bash_app
def analyze(inputs=[], outputs=[]):
    return f"python analyze.py {inputs[0]} {outputs[0]}"

# Chain file-based tasks
cleaned = preprocess(
    inputs=[File('raw_data.csv')],
    outputs=[File('cleaned_data.csv')]
)

results = analyze(
    inputs=[cleaned.outputs[0]],  # Use output from preprocess
    outputs=[File('analysis.txt')]
)
```

## Monitoring

```python
from parsl.monitoring import MonitoringHub

config = Config(
    executors=[...],
    monitoring=MonitoringHub(
        hub_address="localhost",
        hub_port=55055,
        logging_endpoint='sqlite:///monitoring.db'
    )
)

# View at http://localhost:8080
```

## Best Practices

1. **Start Local:** Test with local executor before HPC
2. **Checkpoint:** Use Parsl checkpointing for long workflows
3. **Resource Estimation:** Specify task resource needs
4. **Error Handling:** Catch exceptions in app functions

## When to Use Parsl

✓ HPC cluster workflows
✓ Implicit dataflow parallelism
✓ Python-native task definition
✓ Cross-platform execution
✓ Jupyter notebook workflows

## Limitations

❌ Simple caching (use joblib)
❌ Web-based monitoring (use Prefect)
❌ Materials-specific (use quacc)
❌ Just local parallelism (use joblib)

## References

- Official: https://parsl-project.org/
- Examples: `../examples/hpc_workflow.py`
