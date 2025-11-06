#!/usr/bin/env python3
"""
Word Counter for Scientific Manuscripts

Count words by section for manuscript length requirements.

Usage:
    python word_counter.py manuscript.docx
    python word_counter.py manuscript.md --exclude-refs
"""

import sys
import argparse
from pathlib import Path


def count_words_text(text):
    """Count words in text."""
    return len(text.split())


def count_words_docx(filepath, exclude_refs=False):
    """Count words in DOCX file."""
    try:
        from docx import Document
    except ImportError:
        print("Error: python-docx not installed. Install with: pip install python-docx")
        sys.exit(1)

    doc = Document(filepath)

    sections = {
        'Total': 0,
        'Abstract': 0,
        'Introduction': 0,
        'Methods': 0,
        'Results': 0,
        'Discussion': 0,
        'Conclusions': 0,
        'Other': 0
    }

    current_section = 'Other'
    in_references = False

    for para in doc.paragraphs:
        text = para.text.strip()

        # Skip empty paragraphs
        if not text:
            continue

        # Detect sections
        text_lower = text.lower()
        if any(x in text_lower for x in ['abstract', 'summary']):
            current_section = 'Abstract'
            continue
        elif 'introduction' in text_lower:
            current_section = 'Introduction'
            continue
        elif any(x in text_lower for x in ['method', 'materials', 'experimental']):
            current_section = 'Methods'
            continue
        elif 'results' in text_lower:
            current_section = 'Results'
            continue
        elif 'discussion' in text_lower:
            current_section = 'Discussion'
            continue
        elif any(x in text_lower for x in ['conclusion', 'summary']):
            current_section = 'Conclusions'
            continue
        elif any(x in text_lower for x in ['references', 'bibliography']):
            in_references = True
            continue

        # Count words
        if exclude_refs and in_references:
            continue

        word_count = count_words_text(text)
        sections[current_section] += word_count
        sections['Total'] += word_count

    return sections


def count_words_markdown(filepath, exclude_refs=False):
    """Count words in markdown file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    sections = {
        'Total': 0,
        'Abstract': 0,
        'Introduction': 0,
        'Methods': 0,
        'Results': 0,
        'Discussion': 0,
        'Conclusions': 0,
        'Other': 0
    }

    current_section = 'Other'
    in_references = False

    for line in content.split('\n'):
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip markdown headings markers for section detection
        if line.startswith('#'):
            heading = line.lstrip('#').strip().lower()

            if any(x in heading for x in ['abstract', 'summary']):
                current_section = 'Abstract'
                continue
            elif 'introduction' in heading:
                current_section = 'Introduction'
                continue
            elif any(x in heading for x in ['method', 'material', 'experimental']):
                current_section = 'Methods'
                continue
            elif 'results' in heading:
                current_section = 'Results'
                continue
            elif 'discussion' in heading:
                current_section = 'Discussion'
                continue
            elif any(x in heading for x in ['conclusion', 'summary']):
                current_section = 'Conclusions'
                continue
            elif any(x in heading for x in ['reference', 'bibliography']):
                in_references = True
                continue

        # Count words
        if exclude_refs and in_references:
            continue

        word_count = count_words_text(line)
        sections[current_section] += word_count
        sections['Total'] += word_count

    return sections


def main():
    parser = argparse.ArgumentParser(description='Count words in manuscript by section')
    parser.add_argument('file', help='Manuscript file (docx or md)')
    parser.add_argument('--exclude-refs', action='store_true',
                       help='Exclude references section from count')

    args = parser.parse_args()

    filepath = Path(args.file)
    if not filepath.exists():
        print(f"Error: File not found: {filepath}")
        sys.exit(1)

    # Determine file type and count
    if filepath.suffix == '.docx':
        sections = count_words_docx(filepath, args.exclude_refs)
    elif filepath.suffix in ['.md', '.markdown']:
        sections = count_words_markdown(filepath, args.exclude_refs)
    else:
        print(f"Error: Unsupported file type: {filepath.suffix}")
        print("Supported: .docx, .md, .markdown")
        sys.exit(1)

    # Print results
    print(f"\n{'='*50}")
    print(f"Word Count for: {filepath.name}")
    print(f"{'='*50}\n")

    print(f"{'Section':<15} {'Words':>10}")
    print(f"{'-'*25}")

    for section in ['Abstract', 'Introduction', 'Methods', 'Results',
                    'Discussion', 'Conclusions', 'Other']:
        if sections[section] > 0:
            print(f"{section:<15} {sections[section]:>10}")

    print(f"{'-'*25}")
    print(f"{'TOTAL':<15} {sections['Total']:>10}\n")


if __name__ == '__main__':
    main()
