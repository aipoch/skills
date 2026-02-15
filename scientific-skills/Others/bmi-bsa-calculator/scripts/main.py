#!/usr/bin/env python3
"""BMI & BSA Calculator - Body metrics for clinical dosing.

Calculates Body Mass Index (BMI) and Body Surface Area (BSA) using standard formulas.
BSA is commonly used for drug dosing in oncology and pediatrics.
"""

import argparse
import json
import math
import sys
from typing import Dict


class BMIBSACalculator:
    """Calculates BMI and BSA using standard medical formulas."""
    
    def calculate(self, weight_kg: float, height_cm: float) -> Dict:
        """Calculate BMI and BSA.
        
        Args:
            weight_kg: Weight in kilograms
            height_cm: Height in centimeters
            
        Returns:
            Dictionary with BMI, BSA, and BMI category
        """
        if weight_kg <= 0 or height_cm <= 0:
            raise ValueError("Weight and height must be positive values")
        
        height_m = height_cm / 100
        
        # BMI calculation (kg/m²)
        bmi = weight_kg / (height_m ** 2)
        
        # BSA calculation (DuBois formula)
        # BSA (m²) = 0.007184 × weight^0.425 × height^0.725
        bsa = 0.007184 * (weight_kg ** 0.425) * (height_cm ** 0.725)
        
        # BMI categories (WHO standards)
        if bmi < 18.5:
            category = "Underweight"
        elif bmi < 25:
            category = "Normal"
        elif bmi < 30:
            category = "Overweight"
        else:
            category = "Obese"
        
        return {
            "bmi": round(bmi, 1),
            "bsa_m2": round(bsa, 2),
            "bmi_category": category,
            "weight_kg": weight_kg,
            "height_cm": height_cm
        }
    
    def calculate_dose(self, bsa: float, standard_dose_mg: float) -> Dict:
        """Calculate drug dose based on BSA.
        
        Args:
            bsa: Body surface area in m²
            standard_dose_mg: Standard dose per m² in mg
            
        Returns:
            Dictionary with calculated dose
        """
        dose = bsa * standard_dose_mg
        return {
            "bsa": bsa,
            "dose_mg": round(dose, 1),
            "dose_per_m2": standard_dose_mg
        }


def main():
    parser = argparse.ArgumentParser(
        description="BMI & BSA Calculator - Body metrics for clinical dosing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate BMI and BSA
  python main.py --weight 70 --height 175
  
  # Calculate with drug dosing
  python main.py --weight 70 --height 175 --dose 100
  
  # Output as JSON
  python main.py --weight 70 --height 175 --format json
        """
    )
    
    parser.add_argument(
        "--weight", "-w",
        type=float,
        required=True,
        help="Weight in kilograms"
    )
    
    parser.add_argument(
        "--height", "-H",
        type=float,
        required=True,
        help="Height in centimeters"
    )
    
    parser.add_argument(
        "--dose", "-d",
        type=float,
        help="Standard drug dose per m² in mg (optional)"
    )
    
    parser.add_argument(
        "--format", "-f",
        type=str,
        choices=["json", "text"],
        default="text",
        help="Output format (default: text)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output file path (optional)"
    )
    
    args = parser.parse_args()
    
    try:
        calc = BMIBSACalculator()
        result = calc.calculate(args.weight, args.height)
        
        # Add dose calculation if provided
        if args.dose:
            dose_info = calc.calculate_dose(result["bsa_m2"], args.dose)
            result["dosing"] = dose_info
        
        # Format output
        if args.format == "json":
            output = json.dumps(result, indent=2)
        else:
            output = f"""
BMI & BSA Calculation
====================
Weight: {result['weight_kg']} kg
Height: {result['height_cm']} cm

BMI: {result['bmi']} kg/m²
Category: {result['bmi_category']}

BSA: {result['bsa_m2']} m²
(DuBois formula)"""
            
            if args.dose:
                output += f"""

Drug Dosing
-----------
Standard dose: {args.dose} mg/m²
Calculated dose: {result['dosing']['dose_mg']} mg
"""
        
        # Output
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Results saved to: {args.output}")
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
