# Prefect Subskill

## Overview

Prefect is a modern workflow orchestration platform for building **production-grade data pipelines** in pure Python. Choose Prefect when you need:
- Complex task dependencies (DAGs)
- Automatic retry logic and error handling
- Real-time monitoring and observability
- Dynamic, data-driven workflows

**Philosophy:** Python-first, no DAG syntax, full language capabilities.

## Core Concepts

### Flows and Tasks

```python
from prefect import flow, task

@task(retries=3, retry_delay_seconds=60)
def extract_data(source):
    """Tasks are reusable units of work."""
    return data

@task
def transform_data(data):
    return transformed_data

@flow(name="ETL Pipeline")
def etl_workflow(source):
    """Flows orchestrate tasks."""
    data = extract_data(source)
    result = transform_data(data)
    return result

# Execute
etl_workflow("database.db")
```

### Key Features

**1. Automatic Retries:**
```python
@task(retries=3, retry_delay_seconds=[10, 60, 300])
def unreliable_api_call():
    # Automatically retries with exponential backoff
    pass
```

**2. Caching:**
```python
from datetime import timedelta

@task(cache_key_fn=task_input_hash, cache_expiration=timedelta(hours=1))
def expensive_computation(params):
    # Results cached for 1 hour
    return result
```

**3. Dynamic Workflows:**
```python
@flow
def dynamic_pipeline(n_jobs):
    # Create tasks dynamically based on runtime data
    results = []
    for i in range(n_jobs):
        result = process_item.submit(i)  # Non-blocking
        results.append(result)

    return [r.result() for r in results]
```

**4. Parallel Execution:**
```python
from prefect import flow, task

@task
def process(item):
    return item ** 2

@flow
def parallel_workflow(items):
    # Submit tasks to run in parallel
    futures = [process.submit(item) for item in items]

    # Wait for all to complete
    results = [f.result() for f in futures]
    return results
```

## Common Patterns

### Pattern 1: ETL Pipeline

```python
from prefect import flow, task

@task
def extract(source):
    # Load data
    return data

@task
def transform(data, rules):
    # Apply transformations
    return transformed

@task
def load(data, destination):
    # Save to database/file
    pass

@flow
def etl(source, dest, rules):
    data = extract(source)
    transformed = transform(data, rules)
    load(transformed, dest)
```

### Pattern 2: Map-Reduce

```python
from prefect import flow, task

@task
def map_function(item):
    return processed_item

@task
def reduce_function(results):
    return combined_result

@flow
def map_reduce(items):
    # Map phase (parallel)
    mapped = [map_function.submit(item) for item in items]
    mapped_results = [f.result() for f in mapped]

    # Reduce phase
    final = reduce_function(mapped_results)
    return final
```

### Pattern 3: Conditional Workflows

```python
@flow
def analysis_pipeline(data, method):
    cleaned = clean_data(data)

    if method == "ml":
        result = ml_analysis(cleaned)
    elif method == "stats":
        result = statistical_analysis(cleaned)
    else:
        raise ValueError(f"Unknown method: {method}")

    return result
```

## Deployment

### Local Development

```python
# Just run the flow
if __name__ == "__main__":
    result = my_flow()
```

### Production Deployment

```python
# Create deployment
from prefect.deployments import Deployment
from prefect.server.schemas.schedules import CronSchedule

deployment = Deployment.build_from_flow(
    flow=etl_workflow,
    name="daily-etl",
    schedule=CronSchedule(cron="0 2 * * *"),  # 2 AM daily
    work_pool_name="default"
)
deployment.apply()
```

### Scheduled Workflows

```python
from prefect import flow
from prefect.server.schemas.schedules import IntervalSchedule
from datetime import timedelta

@flow
def scheduled_task():
    # Your workflow
    pass

# Schedule to run every hour
deployment = Deployment.build_from_flow(
    flow=scheduled_task,
    schedule=IntervalSchedule(interval=timedelta(hours=1))
)
```

## Best Practices

1. **Keep Tasks Small:** One logical unit per task
2. **Use Type Hints:** Helps with debugging
3. **Log Liberally:** Use logger, not print
4. **Handle Failures:** Design for task failures
5. **Test Locally:** Run flows before deployment

## When to Use Prefect

✓ Complex multi-step workflows
✓ Need retry/failure handling
✓ Want monitoring dashboard
✓ Dynamic workflow generation
✓ Event-driven pipelines
✓ Scheduled data pipelines

## Limitations

❌ Simple scripts (use joblib)
❌ Just need caching (use joblib)
❌ HPC-specific features (use Parsl)
❌ Materials science focus (use quacc)

## Migration from joblib

```python
# Before (joblib)
from joblib import Parallel, delayed

def task(x):
    return x ** 2

results = Parallel(n_jobs=4)(delayed(task)(i) for i in range(10))

# After (Prefect)
from prefect import flow, task

@task
def compute(x):
    return x ** 2

@flow
def workflow():
    futures = [compute.submit(i) for i in range(10)]
    return [f.result() for f in futures]

workflow()
```

## References

- Official Docs: https://docs.prefect.io/
- Examples: `../examples/ml_pipeline.py`

## See Also

- `joblib.md` - For simpler needs
- `parsl.md` - For HPC workflows
- `../SKILL.md` - Tool selection guide
