#!/usr/bin/env python3
"""
Readability Analyzer for Scientific Text

Analyzes text readability and provides suggestions.

Usage:
    python readability_analyzer.py text_file.txt
    python readability_analyzer.py --text "Your text here"
"""

import sys
import argparse
import re


def count_syllables(word):
    """Count syllables in a word (simple heuristic)."""
    word = word.lower()
    vowels = 'aeiouy'
    syllable_count = 0
    previous_was_vowel = False

    for char in word:
        is_vowel = char in vowels
        if is_vowel and not previous_was_vowel:
            syllable_count += 1
        previous_was_vowel = is_vowel

    # Adjust for silent 'e'
    if word.endswith('e'):
        syllable_count -= 1
    if syllable_count == 0:
        syllable_count = 1

    return syllable_count


def flesch_reading_ease(text):
    """Calculate Flesch Reading Ease score."""
    sentences = re.split(r'[.!?]+', text)
    sentences = [s.strip() for s in sentences if s.strip()]

    words = text.split()
    word_count = len(words)

    if word_count == 0 or len(sentences) == 0:
        return 0

    syllable_count = sum(count_syllables(word) for word in words)

    avg_sentence_length = word_count / len(sentences)
    avg_syllables_per_word = syllable_count / word_count

    # Flesch Reading Ease formula
    score = 206.835 - (1.015 * avg_sentence_length) - (84.6 * avg_syllables_per_word)

    return score


def analyze_text(text):
    """Analyze text and return metrics."""
    sentences = re.split(r'[.!?]+', text)
    sentences = [s.strip() for s in sentences if s.strip()]

    words = text.split()
    word_count = len(words)

    if word_count == 0:
        return None

    # Basic metrics
    sentence_count = len(sentences)
    avg_sentence_length = word_count / sentence_count if sentence_count > 0 else 0

    # Complex words (3+ syllables)
    complex_words = [w for w in words if count_syllables(w) >= 3]
    complex_word_percent = (len(complex_words) / word_count) * 100

    # Passive voice detection (simple heuristic)
    passive_indicators = ['was', 'were', 'been', 'being']
    passive_count = sum(1 for word in words if word.lower() in passive_indicators)
    passive_percent = (passive_count / sentence_count) * 100 if sentence_count > 0 else 0

    # Readability score
    flesch_score = flesch_reading_ease(text)

    # Problem phrases
    problem_phrases = [
        'very', 'quite', 'rather', 'relatively', 'fairly', 'somewhat',
        'it is important to note that', 'it should be emphasized that',
        'due to the fact that', 'in order to', 'a number of'
    ]
    found_problems = [p for p in problem_phrases if p in text.lower()]

    return {
        'word_count': word_count,
        'sentence_count': sentence_count,
        'avg_sentence_length': avg_sentence_length,
        'complex_word_percent': complex_word_percent,
        'passive_percent': passive_percent,
        'flesch_score': flesch_score,
        'problem_phrases': found_problems
    }


def interpret_flesch_score(score):
    """Interpret Flesch Reading Ease score."""
    if score >= 90:
        return "Very easy (5th grade)"
    elif score >= 80:
        return "Easy (6th grade)"
    elif score >= 70:
        return "Fairly easy (7th grade)"
    elif score >= 60:
        return "Standard (8th-9th grade)"
    elif score >= 50:
        return "Fairly difficult (10th-12th grade)"
    elif score >= 30:
        return "Difficult (college)"
    else:
        return "Very difficult (college graduate)"


def print_results(metrics):
    """Print analysis results."""
    print("\n" + "="*60)
    print("READABILITY ANALYSIS")
    print("="*60)

    print(f"\nBasic Metrics:")
    print(f"  Words: {metrics['word_count']}")
    print(f"  Sentences: {metrics['sentence_count']}")
    print(f"  Average sentence length: {metrics['avg_sentence_length']:.1f} words")

    print(f"\nComplexity:")
    print(f"  Complex words (3+ syllables): {metrics['complex_word_percent']:.1f}%")

    # Flesch Reading Ease
    flesch = metrics['flesch_score']
    interpretation = interpret_flesch_score(flesch)
    print(f"\n  Flesch Reading Ease: {flesch:.1f}")
    print(f"  Interpretation: {interpretation}")

    # For scientific writing, target is 40-60 (difficult but acceptable)
    if flesch < 30:
        print("  ⚠️  Very difficult - consider simplifying")
    elif 30 <= flesch < 40:
        print("  ⚠️  Difficult - appropriate for scientific writing")
    elif 40 <= flesch <= 60:
        print("  ✓ Good range for scientific writing")
    else:
        print("  ✓ Easier to read (may be too simple for scientific audience)")

    print(f"\nStyle:")
    print(f"  Passive voice indicators: {metrics['passive_percent']:.1f}% of sentences")

    if metrics['problem_phrases']:
        print(f"\n⚠️  Problem Phrases Found:")
        for phrase in metrics['problem_phrases'][:5]:  # Show first 5
            print(f"    - '{phrase}'")
        if len(metrics['problem_phrases']) > 5:
            print(f"    ... and {len(metrics['problem_phrases']) - 5} more")

    print("\nSuggestions:")
    if metrics['avg_sentence_length'] > 25:
        print("  - Consider shortening sentences (avg > 25 words)")
    if metrics['complex_word_percent'] > 25:
        print("  - High complex word usage - ensure necessary for precision")
    if metrics['passive_percent'] > 40:
        print("  - High passive voice usage - consider more active voice")
    if metrics['problem_phrases']:
        print("  - Remove vague qualifiers and filler phrases")

    print("")


def main():
    parser = argparse.ArgumentParser(description='Analyze text readability')
    parser.add_argument('file', nargs='?', help='Text file to analyze')
    parser.add_argument('--text', help='Text string to analyze')

    args = parser.parse_args()

    if args.text:
        text = args.text
    elif args.file:
        try:
            with open(args.file, 'r', encoding='utf-8') as f:
                text = f.read()
        except FileNotFoundError:
            print(f"Error: File not found: {args.file}")
            sys.exit(1)
    else:
        parser.print_help()
        sys.exit(1)

    metrics = analyze_text(text)
    if metrics is None:
        print("Error: No text to analyze")
        sys.exit(1)

    print_results(metrics)


if __name__ == '__main__':
    main()
