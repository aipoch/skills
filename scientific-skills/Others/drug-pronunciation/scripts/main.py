#!/usr/bin/env python3
"""Drug Pronunciation - Medical drug name pronunciation assistant with IPA phonetics.

Provides correct pronunciation guides for complex drug generic names with 
IPA phonetic transcriptions and syllable breakdowns.
"""

import argparse
import json
import sys
from typing import Dict, List


class DrugPronunciation:
    """Drug pronunciation guide with IPA and syllable breakdown."""
    
    # Simplified drug database
    DRUG_DATABASE = {
        "metformin": {
            "ipa": "mɛtˈfɔrmɪn",
            "syllables": ["met", "for", "min"],
            "emphasis": "FOR",
            "common_errors": ["met-FOR-min", "MET-for-min"]
        },
        "atorvastatin": {
            "ipa": "əˌtɔrvəˈstætɪn",
            "syllables": ["a", "tor", "va", "sta", "tin"],
            "emphasis": "STA",
            "common_errors": ["a-TOR-va-stat-in", "ator-VAS-ta-tin"]
        },
        "lisinopril": {
            "ipa": "laɪˈsɪnəprɪl",
            "syllables": ["ly", "sin", "o", "pril"],
            "emphasis": "SIN",
            "common_errors": ["LIS-in-o-pril", "ly-SIN-o-pril"]
        },
        "omeprazole": {
            "ipa": "oʊˈmɛprəzoʊl",
            "syllables": ["o", "mep", "ra", "zole"],
            "emphasis": "MEP",
            "common_errors": ["OH-me-pra-zole", "o-me-PRA-zole"]
        },
        "amoxicillin": {
            "ipa": "əˌmɒksɪˈsɪlɪn",
            "syllables": ["a", "mox", "i", "cil", "lin"],
            "emphasis": "CIL",
            "common_errors": ["a-MOX-i-cil-in", "amo-XI-cil-lin"]
        }
    }
    
    def get_pronunciation(self, drug_name: str, format_type: str = "detailed") -> Dict:
        """Get pronunciation guide for a drug.
        
        Args:
            drug_name: Generic or brand drug name
            format_type: Output format (ipa, simple, detailed)
            
        Returns:
            Dictionary with pronunciation information
        """
        drug_lower = drug_name.lower()
        
        if drug_lower not in self.DRUG_DATABASE:
            return {
                "drug_name": drug_name,
                "error": f"Drug '{drug_name}' not found in database",
                "suggestion": "Try common drugs like: metformin, atorvastatin, lisinopril"
            }
        
        drug_info = self.DRUG_DATABASE[drug_lower]
        
        result = {
            "drug_name": drug_name,
            "ipa_transcription": drug_info["ipa"],
            "syllable_breakdown": drug_info["syllables"],
            "emphasis": drug_info["emphasis"],
            "common_errors": drug_info["common_errors"]
        }
        
        if format_type == "simple":
            return {
                "drug_name": drug_name,
                "pronunciation": "-".join(drug_info["syllables"]),
                "emphasis": drug_info["emphasis"]
            }
        elif format_type == "ipa":
            return {
                "drug_name": drug_name,
                "ipa": drug_info["ipa"]
            }
        
        # detailed format - add SSML
        result["audio_ssml"] = self._generate_ssml(drug_name, drug_info)
        return result
    
    def _generate_ssml(self, drug_name: str, drug_info: Dict) -> str:
        """Generate SSML for audio synthesis."""
        syllables = " ".join(drug_info["syllables"])
        emphasis = drug_info["emphasis"]
        return f'<speak><emphasis level="strong">{emphasis}</emphasis> in {syllables}</speak>'
    
    def list_drugs(self) -> List[str]:
        """Return list of available drugs."""
        return list(self.DRUG_DATABASE.keys())


def main():
    parser = argparse.ArgumentParser(
        description="Drug Pronunciation - Medical drug name pronunciation assistant",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get detailed pronunciation for metformin
  python main.py --drug metformin
  
  # Get simple pronunciation
  python main.py --drug atorvastatin --format simple
  
  # Get IPA only
  python main.py --drug lisinopril --format ipa
  
  # List all available drugs
  python main.py --list
        """
    )
    
    parser.add_argument(
        "--drug", "-d",
        type=str,
        help="Drug name (generic or brand)"
    )
    
    parser.add_argument(
        "--format", "-f",
        type=str,
        choices=["ipa", "simple", "detailed"],
        default="detailed",
        help="Output format (default: detailed)"
    )
    
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List all available drugs"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output JSON file path (optional)"
    )
    
    args = parser.parse_args()
    
    pronouncer = DrugPronunciation()
    
    # Handle list command
    if args.list:
        drugs = pronouncer.list_drugs()
        print("Available drugs in database:")
        for drug in drugs:
            print(f"  - {drug}")
        return
    
    # Require drug name if not listing
    if not args.drug:
        parser.print_help()
        sys.exit(1)
    
    # Get pronunciation
    result = pronouncer.get_pronunciation(args.drug, args.format)
    
    # Check for errors
    if "error" in result:
        print(f"Error: {result['error']}", file=sys.stderr)
        if "suggestion" in result:
            print(f"Suggestion: {result['suggestion']}", file=sys.stderr)
        sys.exit(1)
    
    # Output
    output = json.dumps(result, indent=2)
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Pronunciation guide saved to: {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
