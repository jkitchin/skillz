# FireWorks Subskill

## Overview

FireWorks is a **production workflow engine** for complex, long-running scientific workflows. Use for:
- Large-scale materials science production systems
- Complex failure recovery needed
- Thousands of interdependent jobs
- Materials Project-style workflows

## Architecture

- **LaunchPad:** Central MongoDB database managing workflows
- **FireWorkers:** Distributed agents executing jobs
- **Firetasks:** Atomic units of work

## Quick Start

```python
from fireworks import Firework, Workflow, LaunchPad, ScriptTask

# Define task
task = ScriptTask.from_str('echo "Hello"')

# Create Firework
fw = Firework(task, name="hello_fw")

# Create workflow
wf = Workflow([fw], name="test_workflow")

# Submit to LaunchPad
launchpad = LaunchPad()
launchpad.add_wf(wf)

# Execute
from fireworks.core.rocket_launcher import rapidfire
rapidfire(launchpad)
```

## Integration with atomate2

```python
from atomate2.vasp.flows.core import RelaxBandStructureFlow
from jobflow.managers.fireworks import flow_to_workflow

# Create atomate2 flow
flow = RelaxBandStructureFlow(structure)

# Convert to FireWorks
wf = flow_to_workflow(flow)

# Submit
launchpad.add_wf(wf)
```

## When to Use

✓ Production materials workflows
✓ Complex dependencies and retries
✓ Existing Materials Project infrastructure
✓ Thousands of jobs

## Limitations

❌ Simple scripts (use joblib)
❌ Learning/prototyping (too complex)
❌ Modern Python workflows (use Prefect)

## References

- Docs: https://materialsproject.github.io/fireworks/
- atomate2: https://github.com/materialsproject/atomate2
