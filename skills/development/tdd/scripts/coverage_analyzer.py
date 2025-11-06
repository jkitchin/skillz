#!/usr/bin/env python3
"""
Coverage Analyzer

Parse coverage reports and identify untested code to prioritize next tests.

Usage:
    python coverage_analyzer.py analyze coverage.xml
    python coverage_analyzer.py analyze .coverage --format coverage.py
    python coverage_analyzer.py suggest coverage.xml --min-coverage 80

Supported formats:
    - coverage.xml (Cobertura XML format)
    - .coverage (coverage.py SQLite format)
    - coverage.json (coverage.py JSON format)
"""

import sys
import argparse
import json
from pathlib import Path
from collections import defaultdict


def parse_cobertura_xml(filepath):
    """Parse Cobertura XML coverage report."""
    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        print("Error: xml.etree not available")
        sys.exit(1)

    tree = ET.parse(filepath)
    root = tree.getroot()

    coverage_data = {
        'files': {},
        'summary': {
            'line_rate': float(root.get('line-rate', 0)),
            'branch_rate': float(root.get('branch-rate', 0)),
            'lines_covered': 0,
            'lines_total': 0
        }
    }

    for package in root.findall('.//package'):
        for class_elem in package.findall('.//class'):
            filename = class_elem.get('filename')
            lines = class_elem.findall('.//line')

            covered_lines = []
            uncovered_lines = []

            for line in lines:
                line_num = int(line.get('number'))
                hits = int(line.get('hits', 0))

                if hits > 0:
                    covered_lines.append(line_num)
                else:
                    uncovered_lines.append(line_num)

            coverage_data['files'][filename] = {
                'covered': covered_lines,
                'uncovered': uncovered_lines,
                'coverage': len(covered_lines) / len(lines) if lines else 0
            }

            coverage_data['summary']['lines_covered'] += len(covered_lines)
            coverage_data['summary']['lines_total'] += len(lines)

    return coverage_data


def parse_coverage_json(filepath):
    """Parse coverage.py JSON format."""
    with open(filepath, 'r') as f:
        data = json.load(f)

    coverage_data = {
        'files': {},
        'summary': {
            'line_rate': data['totals']['percent_covered'] / 100,
            'lines_covered': data['totals']['covered_lines'],
            'lines_total': data['totals']['num_statements']
        }
    }

    for filename, file_data in data['files'].items():
        executed = file_data.get('executed_lines', [])
        missing = file_data.get('missing_lines', [])

        coverage_data['files'][filename] = {
            'covered': executed,
            'uncovered': missing,
            'coverage': file_data['summary']['percent_covered'] / 100
        }

    return coverage_data


def parse_coverage_file(filepath, format_type=None):
    """Auto-detect and parse coverage file."""
    path = Path(filepath)

    if format_type == 'xml' or path.suffix == '.xml':
        return parse_cobertura_xml(filepath)
    elif format_type == 'json' or path.suffix == '.json':
        return parse_coverage_json(filepath)
    else:
        # Try to auto-detect
        if path.suffix == '.xml':
            return parse_cobertura_xml(filepath)
        elif path.suffix == '.json':
            return parse_coverage_json(filepath)
        else:
            print(f"Error: Unknown coverage format for {filepath}")
            print("Specify format with --format xml or --format json")
            sys.exit(1)


def analyze_coverage(coverage_data, min_coverage=80):
    """Analyze coverage data and identify gaps."""
    print("\n=== Coverage Summary ===")
    summary = coverage_data['summary']

    overall_coverage = (summary['lines_covered'] / summary['lines_total'] * 100
                       if summary['lines_total'] > 0 else 0)

    print(f"Overall Coverage: {overall_coverage:.1f}%")
    print(f"Lines Covered: {summary['lines_covered']}/{summary['lines_total']}")

    if 'branch_rate' in summary:
        print(f"Branch Coverage: {summary['branch_rate'] * 100:.1f}%")

    print(f"\nTarget Coverage: {min_coverage}%")

    if overall_coverage >= min_coverage:
        print(f"✅ Coverage target met!")
    else:
        gap = min_coverage - overall_coverage
        print(f"⚠️  Coverage gap: {gap:.1f}%")

    # Files with low coverage
    print("\n=== Files Below Target ===")
    low_coverage_files = []

    for filename, data in coverage_data['files'].items():
        file_coverage = data['coverage'] * 100
        if file_coverage < min_coverage:
            low_coverage_files.append((filename, file_coverage, data))

    low_coverage_files.sort(key=lambda x: x[1])  # Sort by coverage

    if not low_coverage_files:
        print("All files meet coverage target! 🎉")
    else:
        for filename, coverage, data in low_coverage_files[:10]:  # Top 10
            uncovered_count = len(data['uncovered'])
            print(f"\n{filename}")
            print(f"  Coverage: {coverage:.1f}%")
            print(f"  Uncovered lines: {uncovered_count}")

            # Show first few uncovered line ranges
            uncovered = sorted(data['uncovered'])
            if uncovered:
                ranges = group_line_ranges(uncovered)
                print(f"  Lines: {', '.join(ranges[:5])}")
                if len(ranges) > 5:
                    print(f"  ... and {len(ranges) - 5} more ranges")

    return low_coverage_files


def group_line_ranges(lines):
    """Group consecutive lines into ranges."""
    if not lines:
        return []

    ranges = []
    start = lines[0]
    end = lines[0]

    for line in lines[1:]:
        if line == end + 1:
            end = line
        else:
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = line
            end = line

    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")

    return ranges


def suggest_next_tests(coverage_data, low_coverage_files):
    """Suggest what to test next."""
    print("\n=== Suggested Next Tests ===")

    if not low_coverage_files:
        print("Coverage is good! Consider:")
        print("  • Edge cases and error conditions")
        print("  • Integration tests")
        print("  • End-to-end scenarios")
        return

    # Prioritize by impact
    for i, (filename, coverage, data) in enumerate(low_coverage_files[:5], 1):
        print(f"\n{i}. {Path(filename).name}")
        print(f"   Priority: {'🔴 High' if coverage < 50 else '🟡 Medium'}")
        print(f"   Current: {coverage:.1f}% coverage")

        uncovered = sorted(data['uncovered'])
        if uncovered:
            ranges = group_line_ranges(uncovered)
            print(f"   Focus on lines: {', '.join(ranges[:3])}")

        # Generic suggestions based on file type
        if 'test' not in filename.lower():
            print("   Suggestions:")
            if coverage < 30:
                print("     • Start with happy path tests")
                print("     • Test main public functions")
            elif coverage < 60:
                print("     • Add error condition tests")
                print("     • Test edge cases")
            else:
                print("     • Test remaining branches")
                print("     • Add integration tests")


def main():
    parser = argparse.ArgumentParser(description='Analyze test coverage reports')

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze coverage report')
    analyze_parser.add_argument('file', help='Coverage report file')
    analyze_parser.add_argument('--format', choices=['xml', 'json'],
                               help='Coverage format (auto-detected if not specified)')
    analyze_parser.add_argument('--min-coverage', type=float, default=80,
                               help='Minimum coverage target (default: 80)')

    # Suggest command
    suggest_parser = subparsers.add_parser('suggest', help='Suggest next tests')
    suggest_parser.add_argument('file', help='Coverage report file')
    suggest_parser.add_argument('--format', choices=['xml', 'json'],
                               help='Coverage format (auto-detected if not specified)')
    suggest_parser.add_argument('--min-coverage', type=float, default=80,
                               help='Minimum coverage target (default: 80)')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Parse coverage file
    coverage_data = parse_coverage_file(args.file, args.format)

    # Analyze
    low_coverage_files = analyze_coverage(coverage_data, args.min_coverage)

    # Suggest if requested
    if args.command == 'suggest':
        suggest_next_tests(coverage_data, low_coverage_files)


if __name__ == '__main__':
    main()
