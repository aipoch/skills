#!/usr/bin/env python3
"""Comparison Table Generator - Markdown table creator for research comparisons."""

import argparse
import json
import sys
from typing import List


class ComparisonTableGen:
    """Generates comparison tables in Markdown format."""
    
    def generate(self, items: List[str], attributes: List[str]) -> dict:
        """Generate comparison table.
        
        Args:
            items: List of items to compare (e.g., ["Drug A", "Drug B"])
            attributes: List of comparison attributes (e.g., ["Mechanism", "Dose"])
            
        Returns:
            Dictionary with markdown table and metadata
        """
        if not items or not attributes:
            raise ValueError("Both items and attributes must be non-empty lists")
        
        # Header
        header = "| Feature | " + " | ".join(items) + " |"
        separator = "|" + "|".join(["---"] * (len(items) + 1)) + "|"
        
        # Rows
        rows = []
        for attr in attributes:
            row = f"| {attr} |" + " | " * len(items)
            rows.append(row)
        
        markdown = "\n".join([header, separator] + rows)
        
        return {
            "markdown_table": markdown,
            "items": items,
            "attributes": attributes
        }


def parse_list(text: str) -> List[str]:
    """Parse comma-separated list into list of strings."""
    return [item.strip() for item in text.split(",") if item.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Comparison Table Generator - Create Markdown comparison tables",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two drugs
  python main.py --items "Drug A,Drug B" --attributes "Mechanism,Dose,Side Effects"
  
  # Compare three treatments
  python main.py --items "Surgery,Chemo,Radiation" --attributes "Cost,Efficacy,Side Effects" --output table.json
        """
    )
    
    parser.add_argument(
        "--items", "-i",
        type=str,
        required=True,
        help='Items to compare (comma-separated, e.g., "Drug A,Drug B")'
    )
    
    parser.add_argument(
        "--attributes", "-a",
        type=str,
        required=True,
        help='Comparison attributes (comma-separated, e.g., "Mechanism,Dose")'
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output JSON file path (optional, prints to stdout if not specified)"
    )
    
    args = parser.parse_args()
    
    try:
        # Parse lists
        items = parse_list(args.items)
        attributes = parse_list(args.attributes)
        
        # Generate table
        gen = ComparisonTableGen()
        result = gen.generate(items, attributes)
        
        # Output
        output = json.dumps(result, indent=2)
        
        if args.output:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.write(output)
            print(f"Table saved to: {args.output}")
        else:
            print(output)
            
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
