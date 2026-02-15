#!/usr/bin/env python3
"""Concept Explainer - Medical concept analogies for accessible explanations.

Uses everyday analogies to explain complex medical concepts to different audiences
(patients, children, students).
"""

import argparse
import json
import sys
from typing import Dict, Optional


class ConceptExplainer:
    """Explains medical concepts using analogies."""
    
    ANALOGIES = {
        "thrombosis": {
            "analogy": "Like a traffic jam in a highway",
            "explanation": "Blood clots form when blood flow slows or stops, similar to how traffic jams occur when cars can't move freely."
        },
        "immune system": {
            "analogy": "Like a security system in a building",
            "explanation": "The immune system patrols the body looking for intruders (pathogens), just like security guards check for unauthorized visitors."
        },
        "antibiotic resistance": {
            "analogy": "Like weeds becoming resistant to weed killer",
            "explanation": "Bacteria evolve to survive antibiotics, similar to how weeds adapt to survive repeated herbicide use."
        },
        "inflammation": {
            "analogy": "Like a fire alarm system",
            "explanation": "Inflammation is the body's alarm response to injury or infection, calling emergency responders (immune cells) to the scene."
        },
        "blood pressure": {
            "analogy": "Like water pressure in garden hoses",
            "explanation": "Blood pressure is the force of blood pushing against artery walls, similar to water pressure in hoses. Too high can damage the pipes."
        }
    }
    
    def explain(self, concept: str, audience: str = "patient") -> Dict:
        """Explain concept with analogy.
        
        Args:
            concept: Medical concept to explain
            audience: Target audience (child, patient, student)
            
        Returns:
            Dictionary with explanation and analogy
        """
        concept_lower = concept.lower()
        
        if concept_lower in self.ANALOGIES:
            data = self.ANALOGIES[concept_lower]
            
            # Adjust explanation based on audience
            explanation = data["explanation"]
            if audience == "child":
                explanation = self._simplify_for_child(explanation)
            elif audience == "student":
                explanation = self._add_detail_for_student(explanation)
            
            return {
                "concept": concept,
                "explanation": explanation,
                "analogy": data["analogy"],
                "audience": audience
            }
        
        return {
            "concept": concept,
            "explanation": f"{concept} is a medical condition or physiological process.",
            "analogy": "Analogy not available in database.",
            "audience": audience
        }
    
    def _simplify_for_child(self, explanation: str) -> str:
        """Simplify explanation for children."""
        # Simple simplification - in practice would be more sophisticated
        return explanation.split(".")[0] + "."
    
    def _add_detail_for_student(self, explanation: str) -> str:
        """Add technical detail for students."""
        return explanation + " This involves complex physiological mechanisms at the cellular level."
    
    def list_concepts(self) -> list:
        """Return list of available concepts."""
        return list(self.ANALOGIES.keys())


def main():
    parser = argparse.ArgumentParser(
        description="Concept Explainer - Explain medical concepts using analogies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Explain thrombosis to a patient
  python main.py --concept "thrombosis"
  
  # Explain to a child
  python main.py --concept "immune system" --audience child
  
  # Explain to a medical student
  python main.py --concept "antibiotic resistance" --audience student
  
  # List all available concepts
  python main.py --list
        """
    )
    
    parser.add_argument(
        "--concept", "-c",
        type=str,
        help="Medical concept to explain (e.g., 'thrombosis', 'immune system')"
    )
    
    parser.add_argument(
        "--audience", "-a",
        type=str,
        choices=["child", "patient", "student"],
        default="patient",
        help="Target audience (default: patient)"
    )
    
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List all available concepts"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output JSON file path (optional)"
    )
    
    args = parser.parse_args()
    
    explainer = ConceptExplainer()
    
    # Handle list command
    if args.list:
        concepts = explainer.list_concepts()
        print("Available medical concepts:")
        for concept in concepts:
            print(f"  - {concept}")
        return
    
    # Require concept if not listing
    if not args.concept:
        parser.print_help()
        sys.exit(1)
    
    # Generate explanation
    result = explainer.explain(args.concept, args.audience)
    
    # Output
    output = json.dumps(result, indent=2)
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Explanation saved to: {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
