#!/usr/bin/env python3
"""
Authorship CRediT Generator
Generates standardized author contribution statements following CRediT taxonomy

ID: 160
"""

import argparse
import json
import sys
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict


# CRediT 14 roles standard definition
CREDIT_ROLES = {
    "C1": {
        "en": "Conceptualization",
        "zh": "Conceptualization",
        "desc": "Ideas; formulation or evolution of overarching research goals and aims"
    },
    "C2": {
        "en": "Data curation",
        "zh": "Data curation",
        "desc": "Management activities to annotate, scrub and maintain research data"
    },
    "C3": {
        "en": "Formal analysis",
        "zh": "Formal analysis",
        "desc": "Application of statistical, mathematical, computational techniques"
    },
    "C4": {
        "en": "Funding acquisition",
        "zh": "Funding acquisition",
        "desc": "Acquisition of the financial support for the project"
    },
    "C5": {
        "en": "Investigation",
        "zh": "Investigation",
        "desc": "Conducting a research and investigation process"
    },
    "C6": {
        "en": "Methodology",
        "zh": "Methodology",
        "desc": "Development or design of methodology"
    },
    "C7": {
        "en": "Project administration",
        "zh": "Project administration",
        "desc": "Management and coordination responsibility for the research"
    },
    "C8": {
        "en": "Resources",
        "zh": "Resources",
        "desc": "Provision of study materials, reagents, materials, patients, etc."
    },
    "C9": {
        "en": "Software",
        "zh": "Software",
        "desc": "Programming, software development"
    },
    "C10": {
        "en": "Supervision",
        "zh": "Supervision",
        "desc": "Oversight and leadership responsibility"
    },
    "C11": {
        "en": "Validation",
        "zh": "Validation",
        "desc": "Verification and replication of results"
    },
    "C12": {
        "en": "Visualization",
        "zh": "Visualization",
        "desc": "Preparation of figures and data presentation"
    },
    "C13": {
        "en": "Writing – original draft",
        "zh": "Writing – original draft",
        "desc": "Preparation and creation of the published work"
    },
    "C14": {
        "en": "Writing – review & editing",
        "zh": "Writing – review & editing",
        "desc": "Critical review, commentary or revision"
    }
}


@dataclass
class Author:
    """Author information class"""
    name: str
    roles: List[str]
    affiliation: str = ""
    
    def validate_roles(self) -> List[str]:
        """Validate if role codes are valid"""
        invalid = [r for r in self.roles if r not in CREDIT_ROLES]
        return invalid
    
    def get_role_names(self, lang: str = "en") -> List[str]:
        """Get list of role names"""
        return [CREDIT_ROLES[r][lang] for r in self.roles if r in CREDIT_ROLES]


@dataclass
class Contribution:
    """Contribution statement class"""
    authors: List[Author]
    equal_contribution: List[str] = None
    corresponding: List[str] = None
    language: str = "en"
    
    def __post_init__(self):
        if self.equal_contribution is None:
            self.equal_contribution = []
        if self.corresponding is None:
            self.corresponding = []


class CRediTGenerator:
    """CRediT contribution statement generator"""
    
    def __init__(self, contribution: Contribution):
        self.contribution = contribution
        self._validate()
    
    def _validate(self):
        """Validate input data"""
        for author in self.contribution.authors:
            invalid = author.validate_roles()
            if invalid:
                raise ValueError(f"Author '{author.name}' has invalid role codes: {', '.join(invalid)}")
    
    def generate_text(self) -> str:
        """Generate text format contribution statement"""
        lines = []
        lang = self.contribution.language
        
        # Header
        lines.append("Author Contributions")
        lines.append("=" * 20)
        lines.append("")
        
        # Author contribution list
        for author in self.contribution.authors:
            role_names = author.get_role_names(lang)
            if author.affiliation:
                lines.append(f"{author.name} ({author.affiliation}): {', '.join(role_names)}")
            else:
                lines.append(f"{author.name}: {', '.join(role_names)}")
        
        lines.append("")
        
        # Equal contribution note
        if self.contribution.equal_contribution:
            names = ", ".join(self.contribution.equal_contribution)
            lines.append(f"*{names} contributed equally to this work")
            lines.append("")
        
        # Corresponding author note
        if self.contribution.corresponding:
            names = ", ".join(self.contribution.corresponding)
            lines.append(f"Corresponding author(s): {names}")
        
        return "\n".join(lines)
    
    def generate_bilingual(self) -> str:
        """Generate bilingual contribution statement"""
        lines = []
        lines.append("Author Contributions")
        lines.append("=" * 40)
        lines.append("")
        
        for author in self.contribution.authors:
            en_roles = author.get_role_names("en")
            
            if author.affiliation:
                lines.append(f"{author.name} ({author.affiliation}):")
            else:
                lines.append(f"{author.name}:")
            
            lines.append(f"  {', '.join(en_roles)}")
            lines.append("")
        
        # Equal contribution
        if self.contribution.equal_contribution:
            names = ", ".join(self.contribution.equal_contribution)
            lines.append(f"*Equal contribution: {names}")
            lines.append("")
        
        # Corresponding
        if self.contribution.corresponding:
            names = ", ".join(self.contribution.corresponding)
            lines.append(f"Corresponding: {names}")
        
        return "\n".join(lines)
    
    def generate_json(self) -> str:
        """Generate JSON format contribution statement"""
        data = {
            "authors": [
                {
                    "name": a.name,
                    "affiliation": a.affiliation,
                    "roles": [
                        {
                            "code": r,
                            "name_en": CREDIT_ROLES[r]["en"]
                        }
                        for r in a.roles
                    ]
                }
                for a in self.contribution.authors
            ],
            "equal_contribution": self.contribution.equal_contribution,
            "corresponding_authors": self.contribution.corresponding
        }
        return json.dumps(data, ensure_ascii=False, indent=2)
    
    def generate_xml(self) -> str:
        """Generate CRediT XML format contribution statement"""
        lines = []
        lines.append('<?xml version="1.0" encoding="UTF-8"?>')
        lines.append('<contrib-group>')
        
        for author in self.contribution.authors:
            lines.append('  <contrib contrib-type="author">')
            lines.append(f'    <name>{author.name}</name>')
            if author.affiliation:
                lines.append(f'    <aff>{author.affiliation}</aff>')
            
            for role_code in author.roles:
                role = CREDIT_ROLES[role_code]
                lines.append(f'    <role vocab="credit" vocab-identifier="http://credit.niso.org/">')
                lines.append(f'      {role["en"]}')
                lines.append(f'    </role>')
            
            # Corresponding author
            if author.name in self.contribution.corresponding:
                lines.append('    <email>Corresponding Author</email>')
            
            # Equal contribution
            if author.name in self.contribution.equal_contribution:
                lines.append('    <xref ref-type="equal"/>')
            
            lines.append('  </contrib>')
        
        lines.append('</contrib-group>')
        return "\n".join(lines)
    
    def generate(self, format_type: str = "text") -> str:
        """Generate output based on format type"""
        if format_type == "json":
            return self.generate_json()
        elif format_type == "xml":
            return self.generate_xml()
        elif format_type == "bilingual":
            return self.generate_bilingual()
        else:
            return self.generate_text()


def parse_short_format(text: str) -> List[Author]:
    """Parse short format author information"""
    # Format: Name1:Role1,Role2,...|Name2:Role3,Role4,...
    authors = []
    parts = text.split("|")
    
    for part in parts:
        if ":" not in part:
            continue
        name, roles_str = part.split(":", 1)
        roles = [r.strip() for r in roles_str.split(",")]
        authors.append(Author(name=name.strip(), roles=roles))
    
    return authors


def interactive_mode():
    """Interactive mode"""
    print("=" * 50)
    print("CRediT Author Contribution Generator")
    print("=" * 50)
    print()
    
    # Display roles
    print("CRediT 14 Standard Roles:")
    print("-" * 30)
    for code, info in CREDIT_ROLES.items():
        print(f"  {code}: {info['en']}")
    print()
    
    # Number of authors
    num_authors = int(input("Enter number of authors: "))
    
    authors = []
    for i in range(num_authors):
        print(f"\n--- Author {i+1} ---")
        name = input("Name: ")
        affiliation = input("Affiliation (optional): ")
        print("Role codes (comma-separated, e.g., C1,C5,C13):")
        roles_input = input("Roles: ")
        roles = [r.strip() for r in roles_input.split(",")]
        
        authors.append(Author(name=name, roles=roles, affiliation=affiliation))
    
    # Equal contribution
    print("\n--- Equal Contribution ---")
    equal_input = input("Equal contribution authors (comma-separated, leave blank if none): ")
    equal_contribution = [n.strip() for n in equal_input.split(",")] if equal_input else []
    
    # Corresponding authors
    print("\n--- Corresponding Authors ---")
    corres_input = input("Corresponding author names (comma-separated): ")
    corresponding = [n.strip() for n in corres_input.split(",")] if corres_input else []
    
    # Output format
    print("\n--- Output Format ---")
    print("1. Text")
    print("2. JSON")
    print("3. XML")
    format_choice = input("Select (1-3): ")
    
    format_map = {"1": "text", "2": "json", "3": "xml"}
    format_type = format_map.get(format_choice, "text")
    
    # Generate output
    contribution = Contribution(
        authors=authors,
        equal_contribution=equal_contribution,
        corresponding=corresponding,
        language="en"
    )
    
    generator = CRediTGenerator(contribution)
    output = generator.generate(format_type)
    
    print("\n" + "=" * 50)
    print("Generated Contribution Statement:")
    print("=" * 50)
    print(output)
    
    # Save option
    save = input("\nSave to file? (y/n): ")
    if save.lower() == 'y':
        filename = input("Filename: ")
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Saved to: {filename}")


def main():
    parser = argparse.ArgumentParser(
        description='CRediT Author Contribution Statement Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --interactive
  %(prog)s --authors "John:C1,C5,C13|Jane:C2,C6,C9,C14"
  %(prog)s --input team.json --format json
        """
    )
    
    parser.add_argument('--authors', type=str, 
                       help='Author roles shorthand (format: Name1:Role1,Role2,...|Name2:...)')
    parser.add_argument('--input', '-i', type=str,
                       help='Input JSON file path')
    parser.add_argument('--output', '-o', type=str,
                       help='Output file path (default: stdout)')
    parser.add_argument('--format', '-f', type=str, 
                       choices=['text', 'json', 'xml', 'bilingual'],
                       default='text',
                       help='Output format (default: text)')
    parser.add_argument('--language', '-l', type=str,
                       choices=['en'],
                       default='en',
                       help='Output language (default: en)')
    parser.add_argument('--interactive', action='store_true',
                       help='Interactive mode')
    parser.add_argument('--corresponding', type=str,
                       help='Corresponding author names (comma-separated)')
    parser.add_argument('--equal', type=str,
                       help='Equal contribution author names (comma-separated)')
    
    args = parser.parse_args()
    
    try:
        # Interactive mode
        if args.interactive:
            interactive_mode()
            return
        
        # Read from file
        if args.input:
            with open(args.input, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            authors = [
                Author(
                    name=a['name'],
                    roles=a['roles'],
                    affiliation=a.get('affiliation', '')
                )
                for a in data.get('authors', [])
            ]
            
            contribution = Contribution(
                authors=authors,
                equal_contribution=data.get('equal_contribution', []),
                corresponding=data.get('corresponding', []),
                language=data.get('language', 'en')
            )
        
        # From command line arguments
        elif args.authors:
            authors = parse_short_format(args.authors)
            equal_contribution = [n.strip() for n in args.equal.split(",")] if args.equal else []
            corresponding = [n.strip() for n in args.corresponding.split(",")] if args.corresponding else []
            
            contribution = Contribution(
                authors=authors,
                equal_contribution=equal_contribution,
                corresponding=corresponding,
                language=args.language
            )
        
        else:
            parser.print_help()
            return
        
        # Generate output
        generator = CRediTGenerator(contribution)
        output = generator.generate(args.format)
        
        # Output result
        if args.output:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.write(output)
            print(f"Saved to: {args.output}")
        else:
            print(output)
    
    except FileNotFoundError:
        print(f"Error: File not found '{args.input}'", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: JSON parsing failed - {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
