#!/usr/bin/env python3
"""Cover Letter Drafter - Professional cover letter generator for academic contexts.

Generates tailored cover letters for journal submissions, job applications, and fellowships.
"""

import argparse
import json
import sys
from typing import Dict, List


class CoverLetterDrafter:
    """Generates cover letters for academic and medical contexts."""
    
    TEMPLATES = {
        "journal": """Dear Editor,

I am pleased to submit our manuscript titled "{title}" for consideration for publication in {journal}.

{key_points}

This work represents {significance}. We believe it aligns with the scope and readership of your journal.

Thank you for your consideration.

Sincerely,
{author}""",
        "job": """Dear Hiring Committee,

I am writing to express my strong interest in the {position} position at {institution}.

{key_points}

I am excited about the opportunity to contribute to your team.

Thank you for your consideration.

Sincerely,
{applicant}""",
        "fellowship": """Dear Fellowship Committee,

I am writing to apply for the {fellowship} fellowship at {institution}.

{key_points}

I am eager to contribute to your research program.

Thank you for your consideration.

Sincerely,
{applicant}"""
    }
    
    def draft(self, purpose: str, recipient: str, key_points: List[str], **kwargs) -> Dict:
        """Generate cover letter.
        
        Args:
            purpose: Type of cover letter (journal, job, fellowship)
            recipient: Target journal or institution
            key_points: Main points to highlight
            **kwargs: Additional template variables
            
        Returns:
            Dictionary with cover letter and metadata
        """
        template = self.TEMPLATES.get(purpose, self.TEMPLATES["job"])
        
        key_points_text = "\n\n".join([f"• {point}" for point in key_points])
        
        letter = template.format(
            recipient=recipient,
            key_points=key_points_text,
            **kwargs
        )
        
        return {
            "cover_letter": letter,
            "purpose": purpose,
            "word_count": len(letter.split())
        }


def main():
    parser = argparse.ArgumentParser(
        description="Cover Letter Drafter - Generate professional cover letters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Journal submission cover letter
  python main.py --purpose journal --recipient "Nature Medicine" \\
    --key-points "Novel findings,Clinical relevance" \\
    --title "Study X" --significance "major advance" --author "Dr. Smith"
  
  # Job application cover letter
  python main.py --purpose job --recipient "Harvard Medical School" \\
    --key-points "10 years experience,Published 20 papers" \\
    --position "Assistant Professor" --applicant "Dr. Jones"
  
  # Fellowship application
  python main.py --purpose fellowship --recipient "NIH" \\
    --key-points "Research excellence,Leadership skills" \\
    --fellowship "K99" --applicant "Dr. Lee"
        """
    )
    
    parser.add_argument(
        "--purpose",
        type=str,
        choices=["journal", "job", "fellowship"],
        default="job",
        help="Type of cover letter (journal, job, fellowship)"
    )
    
    parser.add_argument(
        "--recipient", "-r",
        type=str,
        required=True,
        help="Target journal, institution, or committee"
    )
    
    parser.add_argument(
        "--key-points", "-k",
        type=str,
        required=True,
        help="Comma-separated key points to highlight"
    )
    
    parser.add_argument(
        "--title",
        type=str,
        help="Manuscript title (for journal submissions)"
    )
    
    parser.add_argument(
        "--significance",
        type=str,
        help="Significance statement (for journal submissions)"
    )
    
    parser.add_argument(
        "--author", "--applicant", "-a",
        type=str,
        default="Applicant",
        help="Author or applicant name"
    )
    
    parser.add_argument(
        "--position",
        type=str,
        help="Position title (for job applications)"
    )
    
    parser.add_argument(
        "--fellowship",
        type=str,
        help="Fellowship name (for fellowship applications)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output JSON file path (optional)"
    )
    
    args = parser.parse_args()
    
    drafter = CoverLetterDrafter()
    
    # Parse key points
    key_points = [k.strip() for k in args.key_points.split(",")]
    
    # Build kwargs based on purpose
    kwargs = {"author": args.author}
    
    if args.purpose == "journal":
        if not args.title:
            print("Error: --title required for journal cover letters", file=sys.stderr)
            sys.exit(1)
        kwargs["title"] = args.title
        kwargs["significance"] = args.significance or "important findings"
        kwargs["journal"] = args.recipient
    elif args.purpose == "job":
        if not args.position:
            print("Error: --position required for job cover letters", file=sys.stderr)
            sys.exit(1)
        kwargs["position"] = args.position
        kwargs["institution"] = args.recipient
        kwargs["applicant"] = args.author
    elif args.purpose == "fellowship":
        if not args.fellowship:
            print("Error: --fellowship required for fellowship cover letters", file=sys.stderr)
            sys.exit(1)
        kwargs["fellowship"] = args.fellowship
        kwargs["institution"] = args.recipient
        kwargs["applicant"] = args.author
    
    # Generate cover letter
    result = drafter.draft(args.purpose, args.recipient, key_points, **kwargs)
    
    # Output
    output = json.dumps(result, indent=2)
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Cover letter saved to: {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
