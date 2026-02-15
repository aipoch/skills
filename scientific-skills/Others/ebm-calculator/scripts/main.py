#!/usr/bin/env python3
"""EBM Calculator - Evidence-Based Medicine diagnostic statistics calculator.

Calculates sensitivity, specificity, PPV, NPV, likelihood ratios, and NNT
for clinical decision making and biostatistics education.
"""

import argparse
import json
import sys
from typing import Dict, Optional


class EBMCalculator:
    """Calculates EBM diagnostic test statistics."""
    
    def calculate(self, tp: int, fn: int, tn: int, fp: int, prevalence: float = None) -> Dict:
        """Calculate all EBM metrics from confusion matrix."""
        
        total = tp + fn + tn + fp
        disease_total = tp + fn
        healthy_total = tn + fp
        
        # Sensitivity and Specificity
        sensitivity = tp / disease_total if disease_total > 0 else 0
        specificity = tn / healthy_total if healthy_total > 0 else 0
        
        # PPV and NPV (using prevalence if provided)
        if prevalence is not None:
            # Bayes' theorem for population-adjusted PPV/NPV
            p_disease = prevalence
            p_healthy = 1 - prevalence
            
            ppv = (sensitivity * p_disease) / ((sensitivity * p_disease) + ((1 - specificity) * p_healthy))
            npv = (specificity * p_healthy) / ((specificity * p_healthy) + ((1 - sensitivity) * p_disease))
        else:
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        
        # Likelihood Ratios
        lr_positive = sensitivity / (1 - specificity) if specificity < 1 else float('inf')
        lr_negative = (1 - sensitivity) / specificity if specificity > 0 else float('inf')
        
        # Diagnostic Accuracy
        accuracy = (tp + tn) / total if total > 0 else 0
        
        return {
            "sensitivity": round(sensitivity, 4),
            "specificity": round(specificity, 4),
            "ppv": round(ppv, 4),
            "npv": round(npv, 4),
            "lr_positive": round(lr_positive, 4) if lr_positive != float('inf') else "Infinity",
            "lr_negative": round(lr_negative, 4) if lr_negative != float('inf') else "Infinity",
            "accuracy": round(accuracy, 4),
            "interpretation": self._interpret(sensitivity, specificity),
            "sample_size": total
        }
    
    def calculate_nnt(self, control_event_rate: float, experimental_event_rate: float) -> Dict:
        """Calculate Number Needed to Treat."""
        arr = control_event_rate - experimental_event_rate
        nnt = 1 / arr if arr > 0 else float('inf')
        
        return {
            "absolute_risk_reduction": round(arr, 4),
            "relative_risk_reduction": round(arr / control_event_rate, 4) if control_event_rate > 0 else 0,
            "nnt": round(nnt, 1) if nnt != float('inf') else "Infinity",
            "interpretation": f"Need to treat {int(nnt)} patients to prevent 1 event" if nnt != float('inf') else "No benefit"
        }
    
    def _interpret(self, sens: float, spec: float) -> str:
        """Interpret diagnostic performance."""
        if sens >= 0.95 and spec >= 0.95:
            return "Excellent diagnostic test"
        elif sens >= 0.90 and spec >= 0.90:
            return "Good diagnostic test"
        elif sens < 0.70 or spec < 0.70:
            return "Poor diagnostic test - consider alternatives"
        else:
            return "Moderate diagnostic test"
    
    def pretest_to_posttest(self, pretest_prob: float, lr: float) -> Dict:
        """Convert pre-test to post-test probability using likelihood ratio."""
        pretest_odds = pretest_prob / (1 - pretest_prob)
        posttest_odds = pretest_odds * lr
        posttest_prob = posttest_odds / (1 + posttest_odds)
        
        return {
            "pretest_probability": round(pretest_prob, 4),
            "likelihood_ratio": round(lr, 4),
            "posttest_probability": round(posttest_prob, 4),
            "probability_change": round(posttest_prob - pretest_prob, 4)
        }


def main():
    parser = argparse.ArgumentParser(
        description="EBM Calculator - Evidence-Based Medicine diagnostic statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate diagnostic test statistics
  python main.py --tp 90 --fn 10 --tn 80 --fp 20
  
  # With prevalence adjustment
  python main.py --tp 90 --fn 10 --tn 80 --fp 20 --prevalence 0.1
  
  # Calculate NNT
  python main.py --mode nnt --control-rate 0.25 --experimental-rate 0.15
  
  # Pre-test to post-test probability
  python main.py --mode probability --pretest 0.3 --lr 4.5
        """
    )
    
    parser.add_argument(
        "--mode", "-m",
        type=str,
        choices=["diagnostic", "nnt", "probability"],
        default="diagnostic",
        help="Calculation mode (default: diagnostic)"
    )
    
    # Diagnostic test parameters
    parser.add_argument("--tp", "--true-pos", type=int, help="True positives")
    parser.add_argument("--fn", "--false-neg", type=int, help="False negatives")
    parser.add_argument("--tn", "--true-neg", type=int, help="True negatives")
    parser.add_argument("--fp", "--false-pos", type=int, help="False positives")
    parser.add_argument("--prevalence", "-p", type=float, help="Disease prevalence (0-1)")
    
    # NNT parameters
    parser.add_argument("--control-rate", type=float, help="Control event rate (0-1)")
    parser.add_argument("--experimental-rate", type=float, help="Experimental event rate (0-1)")
    
    # Probability parameters
    parser.add_argument("--pretest", type=float, help="Pre-test probability (0-1)")
    parser.add_argument("--lr", type=float, help="Likelihood ratio")
    
    parser.add_argument("--output", "-o", type=str, help="Output JSON file path (optional)")
    
    args = parser.parse_args()
    
    calc = EBMCalculator()
    
    # Execute based on mode
    if args.mode == "diagnostic":
        # Check required parameters
        if None in [args.tp, args.fn, args.tn, args.fp]:
            print("Error: Diagnostic mode requires --tp, --fn, --tn, --fp", file=sys.stderr)
            sys.exit(1)
        result = calc.calculate(args.tp, args.fn, args.tn, args.fp, args.prevalence)
    
    elif args.mode == "nnt":
        if None in [args.control_rate, args.experimental_rate]:
            print("Error: NNT mode requires --control-rate and --experimental-rate", file=sys.stderr)
            sys.exit(1)
        result = calc.calculate_nnt(args.control_rate, args.experimental_rate)
    
    elif args.mode == "probability":
        if None in [args.pretest, args.lr]:
            print("Error: Probability mode requires --pretest and --lr", file=sys.stderr)
            sys.exit(1)
        result = calc.pretest_to_posttest(args.pretest, args.lr)
    
    # Output
    output = json.dumps(result, indent=2)
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Results saved to: {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
