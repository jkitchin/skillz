#!/usr/bin/env python3
"""
Hypothesis Tracker for Troubleshooting

Track hypotheses, tests, and results during investigation to maintain
clear investigation state and prevent circular testing.

Usage:
    python hypothesis_tracker.py add "Hypothesis text" --priority high
    python hypothesis_tracker.py test HYPOTHESIS_ID "Test description"
    python hypothesis_tracker.py result HYPOTHESIS_ID pass|fail "Notes"
    python hypothesis_tracker.py list
    python hypothesis_tracker.py summary

Examples:
    python hypothesis_tracker.py add "Config file has typo in database URL" --priority high
    python hypothesis_tracker.py test H001 "Check database URL in config.yaml"
    python hypothesis_tracker.py result H001 fail "URL is correct, not the issue"
    python hypothesis_tracker.py list
"""

import json
import sys
import os
from datetime import datetime
from pathlib import Path

DEFAULT_FILE = "troubleshooting_session.json"


def load_session(filename=DEFAULT_FILE):
    """Load or create troubleshooting session."""
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    return {
        "problem": "",
        "started": datetime.now().isoformat(),
        "hypotheses": [],
        "next_id": 1
    }


def save_session(session, filename=DEFAULT_FILE):
    """Save troubleshooting session."""
    with open(filename, 'w') as f:
        json.dump(session, f, indent=2)


def add_hypothesis(session, text, priority="medium"):
    """Add new hypothesis."""
    hyp_id = f"H{session['next_id']:03d}"
    hypothesis = {
        "id": hyp_id,
        "text": text,
        "priority": priority.lower(),
        "status": "untested",
        "added": datetime.now().isoformat(),
        "tests": []
    }
    session["hypotheses"].append(hypothesis)
    session["next_id"] += 1
    return hyp_id


def find_hypothesis(session, hyp_id):
    """Find hypothesis by ID."""
    for hyp in session["hypotheses"]:
        if hyp["id"] == hyp_id:
            return hyp
    return None


def add_test(session, hyp_id, test_description):
    """Add test for hypothesis."""
    hyp = find_hypothesis(session, hyp_id)
    if not hyp:
        return False

    test = {
        "description": test_description,
        "timestamp": datetime.now().isoformat(),
        "result": "pending"
    }
    hyp["tests"].append(test)
    hyp["status"] = "testing"
    return True


def record_result(session, hyp_id, result, notes=""):
    """Record test result."""
    hyp = find_hypothesis(session, hyp_id)
    if not hyp or not hyp["tests"]:
        return False

    # Update last test
    hyp["tests"][-1]["result"] = result.lower()
    hyp["tests"][-1]["notes"] = notes

    # Update hypothesis status
    if result.lower() == "pass":
        hyp["status"] = "confirmed"
    elif result.lower() == "fail":
        hyp["status"] = "falsified"

    return True


def list_hypotheses(session):
    """List all hypotheses."""
    if not session["hypotheses"]:
        print("No hypotheses yet. Add one with: add \"hypothesis text\" --priority high")
        return

    print("\n" + "="*80)
    print("HYPOTHESES")
    print("="*80)

    for hyp in sorted(session["hypotheses"], key=lambda h: {"high": 0, "medium": 1, "low": 2}[h["priority"]]):
        status_icon = {
            "untested": "❓",
            "testing": "🔬",
            "confirmed": "✅",
            "falsified": "❌"
        }[hyp["status"]]

        priority_str = hyp["priority"].upper()

        print(f"\n{status_icon} [{hyp['id']}] {priority_str} - {hyp['status'].upper()}")
        print(f"   {hyp['text']}")

        if hyp["tests"]:
            print(f"   Tests: {len(hyp['tests'])}")
            for i, test in enumerate(hyp["tests"], 1):
                result_icon = {"pending": "⏳", "pass": "✅", "fail": "❌"}[test["result"]]
                print(f"     {i}. {result_icon} {test['description']}")
                if test.get("notes"):
                    print(f"        Notes: {test['notes']}")


def print_summary(session):
    """Print investigation summary."""
    if not session["hypotheses"]:
        print("No hypotheses tracked yet.")
        return

    total = len(session["hypotheses"])
    untested = sum(1 for h in session["hypotheses"] if h["status"] == "untested")
    testing = sum(1 for h in session["hypotheses"] if h["status"] == "testing")
    confirmed = sum(1 for h in session["hypotheses"] if h["status"] == "confirmed")
    falsified = sum(1 for h in session["hypotheses"] if h["status"] == "falsified")

    print("\n" + "="*80)
    print("INVESTIGATION SUMMARY")
    print("="*80)
    print(f"Total Hypotheses: {total}")
    print(f"  ❓ Untested: {untested}")
    print(f"  🔬 Testing: {testing}")
    print(f"  ✅ Confirmed: {confirmed}")
    print(f"  ❌ Falsified: {falsified}")

    if confirmed > 0:
        print("\n🎯 CONFIRMED HYPOTHESES:")
        for hyp in session["hypotheses"]:
            if hyp["status"] == "confirmed":
                print(f"  • {hyp['text']}")

    if untested > 0:
        print("\n📋 NEXT TO TEST:")
        high_priority = [h for h in session["hypotheses"] if h["status"] == "untested" and h["priority"] == "high"]
        if high_priority:
            for hyp in high_priority[:3]:
                print(f"  • [{hyp['id']}] {hyp['text']}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    command = sys.argv[1].lower()
    session = load_session()

    if command == "add":
        if len(sys.argv) < 3:
            print("Usage: add \"hypothesis text\" [--priority high|medium|low]")
            sys.exit(1)

        text = sys.argv[2]
        priority = "medium"

        if len(sys.argv) > 3 and sys.argv[3] == "--priority" and len(sys.argv) > 4:
            priority = sys.argv[4]

        hyp_id = add_hypothesis(session, text, priority)
        save_session(session)
        print(f"✅ Added hypothesis {hyp_id}: {text}")
        print(f"   Priority: {priority.upper()}")

    elif command == "test":
        if len(sys.argv) < 4:
            print("Usage: test HYPOTHESIS_ID \"test description\"")
            sys.exit(1)

        hyp_id = sys.argv[2]
        test_desc = sys.argv[3]

        if add_test(session, hyp_id, test_desc):
            save_session(session)
            print(f"🔬 Added test for {hyp_id}: {test_desc}")
        else:
            print(f"❌ Hypothesis {hyp_id} not found")
            sys.exit(1)

    elif command == "result":
        if len(sys.argv) < 4:
            print("Usage: result HYPOTHESIS_ID pass|fail [\"notes\"]")
            sys.exit(1)

        hyp_id = sys.argv[2]
        result = sys.argv[3]
        notes = sys.argv[4] if len(sys.argv) > 4 else ""

        if result.lower() not in ["pass", "fail"]:
            print("Result must be 'pass' or 'fail'")
            sys.exit(1)

        if record_result(session, hyp_id, result, notes):
            save_session(session)
            icon = "✅" if result.lower() == "pass" else "❌"
            print(f"{icon} Recorded result for {hyp_id}: {result.upper()}")
            if notes:
                print(f"   Notes: {notes}")
        else:
            print(f"❌ Hypothesis {hyp_id} not found or no test recorded")
            sys.exit(1)

    elif command == "list":
        list_hypotheses(session)

    elif command == "summary":
        print_summary(session)

    elif command == "clear":
        if os.path.exists(DEFAULT_FILE):
            os.remove(DEFAULT_FILE)
        print("✅ Cleared troubleshooting session")

    else:
        print(f"Unknown command: {command}")
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
