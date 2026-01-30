#!/usr/bin/env python3
"""
Parameter Sweep Comparison: joblib vs Prefect vs Parsl

Demonstrates how to run the same parameter sweep with different tools.
Shows when to use each tool based on complexity.
"""

import numpy as np
import time


def simulate(param):
    """Simulate an expensive computation."""
    time.sleep(0.1)  # Simulate work
    return param**2


# ============================================================================
# Option 1: joblib (SIMPLEST - Start here!)
# ============================================================================


def run_with_joblib():
    """Use joblib for simple parallel execution."""
    from joblib import Parallel, delayed

    print("\n" + "=" * 60)
    print("Option 1: joblib (Simplest)")
    print("=" * 60)

    parameters = range(10)

    start = time.time()
    results = Parallel(n_jobs=4, verbose=5)(delayed(simulate)(p) for p in parameters)
    elapsed = time.time() - start

    print(f"\nResults: {results}")
    print(f"Time: {elapsed:.2f} seconds")
    print("\nWhen to use: Simple parallel tasks, single machine")


# ============================================================================
# Option 2: Prefect (For complex workflows with monitoring)
# ============================================================================


def run_with_prefect():
    """Use Prefect for workflows with dependencies and monitoring."""
    try:
        from prefect import flow, task

        @task
        def simulate_task(param):
            time.sleep(0.1)
            return param**2

        @flow
        def parameter_sweep():
            parameters = range(10)
            futures = [simulate_task.submit(p) for p in parameters]
            return [f.result() for f in futures]

        print("\n" + "=" * 60)
        print("Option 2: Prefect (Complex workflows)")
        print("=" * 60)

        start = time.time()
        results = parameter_sweep()
        elapsed = time.time() - start

        print(f"\nResults: {results}")
        print(f"Time: {elapsed:.2f} seconds")
        print("\nWhen to use: Complex DAGs, need monitoring, error handling")

    except ImportError:
        print("\nPrefect not installed. Install with: pip install prefect")


# ============================================================================
# Option 3: Parsl (For HPC clusters)
# ============================================================================


def run_with_parsl():
    """Use Parsl for HPC workflows."""
    try:
        import parsl
        from parsl.app.app import python_app
        from parsl.config import Config
        from parsl.executors import HighThroughputExecutor

        config = Config(executors=[HighThroughputExecutor(max_workers=4)])
        parsl.load(config)

        @python_app
        def simulate_app(param):
            import time

            time.sleep(0.1)
            return param**2

        print("\n" + "=" * 60)
        print("Option 3: Parsl (HPC workflows)")
        print("=" * 60)

        parameters = range(10)

        start = time.time()
        futures = [simulate_app(p) for p in parameters]
        results = [f.result() for f in futures]
        elapsed = time.time() - start

        print(f"\nResults: {results}")
        print(f"Time: {elapsed:.2f} seconds")
        print("\nWhen to use: HPC clusters (SLURM/PBS), large-scale parallelism")

        parsl.clear()

    except ImportError:
        print("\nParsl not installed. Install with: pip install parsl")


# ============================================================================
# Main
# ============================================================================


def main():
    """Run parameter sweep with different tools."""

    print("=" * 60)
    print("Parameter Sweep Comparison")
    print("=" * 60)
    print("\nRunning the same parameter sweep with different tools.")
    print("All compute x^2 for x in range(10)")

    # Always works (only requires joblib)
    run_with_joblib()

    # Try Prefect if installed
    run_with_prefect()

    # Try Parsl if installed
    run_with_parsl()

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("\nRecommendation for this task:")
    print("  → Use joblib (simplest, sufficient for this use case)")
    print("\nWhen to escalate:")
    print("  → Prefect: If you need complex dependencies or monitoring")
    print("  → Parsl: If running on HPC cluster with 1000s of tasks")


if __name__ == "__main__":
    main()
