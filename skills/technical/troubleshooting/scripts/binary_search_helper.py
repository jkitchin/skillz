#!/usr/bin/env python3
"""
Binary Search Troubleshooting Helper

Interactive tool to guide binary search troubleshooting process.
Helps narrow down problem space by repeatedly dividing in half.

Usage:
    python binary_search_helper.py

The tool will guide you through:
1. Define the search space (working ↔ broken boundaries)
2. Test the midpoint
3. Determine which half contains the problem
4. Repeat until problem isolated

Examples:
    - Git bisect: Find which commit introduced a bug
    - Config: Find which of many settings causes issue
    - Data: Find which subset of data causes failure
    - Time: Find when problem started occurring
"""

import sys


def binary_search_session():
    """Run interactive binary search session."""
    print("="*80)
    print("BINARY SEARCH TROUBLESHOOTING")
    print("="*80)
    print("\nThis tool helps you efficiently narrow down the cause of a problem")
    print("by repeatedly dividing the search space in half.\n")

    # Define search space
    print("STEP 1: Define your search space")
    print("-"*80)

    space_type = input("What are you searching through?\n"
                      "  1. Git commits\n"
                      "  2. Configuration settings\n"
                      "  3. Data records\n"
                      "  4. Time range\n"
                      "  5. Code sections\n"
                      "  6. Other\n"
                      "Choose (1-6): ").strip()

    type_names = {
        "1": "commits",
        "2": "settings",
        "3": "records",
        "4": "time periods",
        "5": "code lines",
        "6": "items"
    }
    item_name = type_names.get(space_type, "items")

    # Get boundaries
    print(f"\nDefine the boundaries of your search space:")
    good_boundary = input(f"Working (good) boundary: ").strip()
    bad_boundary = input(f"Broken (bad) boundary: ").strip()

    # Get size
    try:
        size = int(input(f"\nHow many {item_name} between {good_boundary} and {bad_boundary}? "))
    except ValueError:
        print("Invalid number. Using 100 as default.")
        size = 100

    print(f"\n✅ Search space defined: {good_boundary} (working) ↔ {bad_boundary} (broken)")
    print(f"   Total {item_name}: {size}")

    # Binary search loop
    left = 0
    right = size
    iteration = 1

    while right - left > 1:
        mid = (left + right) // 2
        percentage = int((mid / size) * 100)

        print("\n" + "="*80)
        print(f"ITERATION {iteration}")
        print("="*80)
        print(f"Current range: positions {left} to {right} ({right - left} {item_name})")
        print(f"Midpoint: position {mid} ({percentage}% through)")

        if space_type == "1":  # Git commits
            print(f"\n📋 Test commit at position {mid}")
            print(f"   (Approximately {percentage}% between {good_boundary} and {bad_boundary})")
        elif space_type == "2":  # Config settings
            print(f"\n📋 Test with first {mid} settings enabled, rest disabled")
        elif space_type == "3":  # Data records
            print(f"\n📋 Test with records {left} to {mid}")
        elif space_type == "4":  # Time range
            print(f"\n📋 Test at midpoint time (position {mid})")
        elif space_type == "5":  # Code
            print(f"\n📋 Test with lines {left} to {mid} active")
        else:
            print(f"\n📋 Test with items {left} to {mid}")

        # Get test result
        while True:
            result = input(f"\nDoes the problem occur at this midpoint? (yes/no/quit): ").strip().lower()
            if result in ['yes', 'no', 'quit', 'q']:
                break
            print("Please enter 'yes', 'no', or 'quit'")

        if result in ['quit', 'q']:
            print("\n👋 Binary search interrupted.")
            sys.exit(0)

        # Update boundaries
        if result == 'yes':
            print(f"   Problem OCCURS → Issue is in first half (positions {left} to {mid})")
            right = mid
        else:
            print(f"   Problem DOESN'T OCCUR → Issue is in second half (positions {mid} to {right})")
            left = mid

        remaining = right - left
        if remaining <= 1:
            print(f"\n🎯 FOUND IT! The issue is at position {right}")
            if space_type == "1":
                print(f"   The problematic commit is at position {right}")
                print(f"   (Use 'git bisect' to find the exact commit)")
            elif space_type == "2":
                print(f"   The problematic setting is number {right}")
            elif space_type == "3":
                print(f"   The problematic record is at position {right}")
            elif space_type == "4":
                print(f"   The problem started at time position {right}")
            elif space_type == "5":
                print(f"   The problematic code is at line {right}")
            else:
                print(f"   The problematic item is at position {right}")
            break
        else:
            print(f"   Narrowed down to {remaining} {item_name}")

        iteration += 1

    # Summary
    print("\n" + "="*80)
    print("BINARY SEARCH COMPLETE")
    print("="*80)
    print(f"Iterations: {iteration}")
    print(f"Initial size: {size} {item_name}")
    print(f"Final range: {right - left} {item_name}")
    print(f"\nEfficiency: Found issue in {iteration} tests instead of {size} tests!")
    print(f"Time saved: {int((1 - iteration/size) * 100)}%")


def main():
    """Main entry point."""
    if len(sys.argv) > 1 and sys.argv[1] in ['-h', '--help', 'help']:
        print(__doc__)
        sys.exit(0)

    try:
        binary_search_session()
    except KeyboardInterrupt:
        print("\n\n👋 Binary search interrupted.")
        sys.exit(0)


if __name__ == "__main__":
    main()
