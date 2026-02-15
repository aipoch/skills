#!/usr/bin/env python3
"""Abstract Trimmer - Compresses abstracts to word limits while preserving key information."""

import argparse
import json
import re
import sys


class AbstractTrimmer:
    """Trims abstracts to meet word limits."""
    
    def trim(self, abstract: str, target_words: int, strategy: str = "balanced") -> dict:
        """Trim abstract to target word count.
        
        Args:
            abstract: Input abstract text
            target_words: Target word count
            strategy: Trimming strategy (conservative/balanced/aggressive)
            
        Returns:
            Dictionary with trimmed abstract and statistics
        """
        original_words = len(abstract.split())
        
        # Remove redundant phrases based on strategy
        trimmed = abstract
        
        if strategy in ["balanced", "aggressive"]:
            redundant = [
                r'\bin this study\b',
                r'\bit is important to note that\b',
                r'\bwe found that\b',
                r'\bthere was\b',
                r'\bit was observed that\b',
                r'\bit should be noted that\b',
                r'\bin conclusion\b',
                r'\bto the best of our knowledge\b',
            ]
            
            for pattern in redundant:
                trimmed = re.sub(pattern, '', trimmed, flags=re.IGNORECASE)
        
        if strategy == "aggressive":
            # More aggressive trimming
            aggressive_patterns = [
                r'\bhowever\b',
                r'\btherefore\b',
                r'\bfurthermore\b',
                r'\bmoreover\b',
            ]
            for pattern in aggressive_patterns:
                trimmed = re.sub(pattern, '', trimmed, flags=re.IGNORECASE)
        
        # Clean up extra spaces
        trimmed = re.sub(r'\s+', ' ', trimmed).strip()
        
        words = trimmed.split()
        if len(words) > target_words:
            # Keep first and last sentences, trim middle
            sentences = trimmed.split('. ')
            if len(sentences) > 2:
                middle_start = target_words // 3
                middle_end = target_words - len(sentences[0].split()) - len(sentences[-1].split()) - 5
                if middle_end > middle_start:
                    middle_words = words[middle_start:middle_end]
                    trimmed = f"{sentences[0]}. {' '.join(middle_words)}... {sentences[-1]}"
                else:
                    trimmed = ' '.join(words[:target_words])
            else:
                trimmed = ' '.join(words[:target_words])
        
        final_words = len(trimmed.split())
        
        return {
            "trimmed_abstract": trimmed,
            "original_words": original_words,
            "final_words": final_words,
            "reduction_percent": round((1 - final_words/original_words) * 100, 1)
        }
    
    def count_words(self, text: str) -> int:
        """Count words in text."""
        return len(text.split())


def main():
    parser = argparse.ArgumentParser(
        description="Abstract Trimmer - Compress abstracts to meet word limits",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --input abstract.txt --target 250
  python main.py --text "Your abstract here..." --target 200 --strategy aggressive
  python main.py --input abstract.txt --check-only
        """
    )
    
    parser.add_argument(
        "--input", "-i",
        type=str,
        help="Input file containing abstract"
    )
    parser.add_argument(
        "--text", "-t",
        type=str,
        help="Abstract text (alternative to --input)"
    )
    parser.add_argument(
        "--target", "-T",
        type=int,
        default=250,
        help="Target word count (default: 250)"
    )
    parser.add_argument(
        "--strategy", "-s",
        choices=["conservative", "balanced", "aggressive"],
        default="balanced",
        help="Trimming strategy (default: balanced)"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output file path"
    )
    parser.add_argument(
        "--check-only", "-c",
        action="store_true",
        help="Only check word count without trimming"
    )
    parser.add_argument(
        "--format",
        choices=["json", "text"],
        default="json",
        help="Output format (default: json)"
    )
    
    args = parser.parse_args()
    
    # Get input text
    if args.input:
        try:
            with open(args.input, 'r', encoding='utf-8') as f:
                abstract = f.read()
        except FileNotFoundError:
            print(f"Error: File not found: {args.input}", file=sys.stderr)
            sys.exit(1)
    elif args.text:
        abstract = args.text
    else:
        # Read from stdin
        abstract = sys.stdin.read()
    
    trimmer = AbstractTrimmer()
    
    # Check only mode
    if args.check_only:
        word_count = trimmer.count_words(abstract)
        excess = word_count - args.target
        print(f"Current word count: {word_count}")
        print(f"Target word count: {args.target}")
        if excess > 0:
            print(f"Excess words: {excess} ({excess/args.target*100:.1f}% over limit)")
        else:
            print(f"Within limit: {-excess} words under target")
        return
    
    # Trim abstract
    result = trimmer.trim(abstract, args.target, args.strategy)
    
    # Output results
    if args.format == "json":
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = result["trimmed_abstract"]
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Trimmed abstract saved to: {args.output}")
        print(f"Reduced from {result['original_words']} to {result['final_words']} words "
              f"({result['reduction_percent']}% reduction)")
    else:
        print(output)
        if args.format == "json":
            print(f"\nReduced from {result['original_words']} to {result['final_words']} words "
                  f"({result['reduction_percent']}% reduction)")


if __name__ == "__main__":
    main()
