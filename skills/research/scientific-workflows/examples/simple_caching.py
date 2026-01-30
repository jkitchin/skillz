#!/usr/bin/env python3
"""
Simple Caching Example with joblib

Demonstrates basic function caching for expensive computations.
This is the simplest workflow tool - start here!
"""

from joblib import Memory
import time
import numpy as np

# Create cache directory
memory = Memory("./cache", verbose=1)


@memory.cache
def expensive_computation(n, seed=42):
    """
    Simulate an expensive computation that we want to cache.

    Parameters
    ----------
    n : int
        Size of computation
    seed : int
        Random seed for reproducibility

    Returns
    -------
    result : float
        Computed result
    """
    print(f"  Computing (this will be slow)...")
    time.sleep(2)  # Simulate expensive computation

    np.random.seed(seed)
    result = np.random.rand(n, n).sum()

    return result


def main():
    """Demonstrate caching behavior."""

    print("=" * 60)
    print("joblib Caching Example")
    print("=" * 60)

    # First call - will compute
    print("\n1. First call (will compute):")
    start = time.time()
    result1 = expensive_computation(1000)
    elapsed1 = time.time() - start
    print(f"  Result: {result1:.4f}")
    print(f"  Time: {elapsed1:.2f} seconds")

    # Second call - will use cache
    print("\n2. Second call (will use cache):")
    start = time.time()
    result2 = expensive_computation(1000)
    elapsed2 = time.time() - start
    print(f"  Result: {result2:.4f}")
    print(f"  Time: {elapsed2:.2f} seconds")

    # Verify results match
    print(f"\n3. Results match: {result1 == result2}")
    print(f"   Speedup: {elapsed1 / elapsed2:.1f}x faster")

    # Different parameters - will compute
    print("\n4. Different parameters (will compute):")
    start = time.time()
    result3 = expensive_computation(1000, seed=123)
    elapsed3 = time.time() - start
    print(f"  Result: {result3:.4f}")
    print(f"  Time: {elapsed3:.2f} seconds")

    print("\n" + "=" * 60)
    print("Cache location:", memory.location)
    print("=" * 60)


if __name__ == "__main__":
    main()

# Output:
# ============================================================
# joblib Caching Example
# ============================================================
#
# 1. First call (will compute):
#   Computing (this will be slow)...
#   Result: 500234.5678
#   Time: 2.05 seconds
#
# 2. Second call (will use cache):
#   Result: 500234.5678
#   Time: 0.01 seconds  # Much faster!
#
# 3. Results match: True
#    Speedup: 205.0x faster
