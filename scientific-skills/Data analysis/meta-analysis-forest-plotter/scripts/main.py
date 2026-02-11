#!/usr/bin/env python3
"""
Meta-analysis Forest Plotter (Pure Python)

Perform meta-analysis calculations and generate publication-quality forest plots.
Implements fixed-effect and random-effects models with heterogeneity statistics.

Features:
- Fixed-effect model (inverse variance weighting)
- Random-effects model (DerSimonian-Laird method)
- Heterogeneity statistics (Q, I², τ²)
- Publication-quality forest plots (PNG/PDF/SVG)
- Support for CSV/JSON input formats
"""

import argparse
import json
import csv
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
import numpy as np
from scipy import stats

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("⚠️  Warning: matplotlib not installed. Install with: pip install matplotlib")


@dataclass
class Study:
    """Represents a single study in the meta-analysis."""
    name: str
    effect: float  # Effect size (e.g., log odds ratio, mean difference)
    lower: float   # Lower confidence interval
    upper: float   # Upper confidence interval
    weight: Optional[float] = None  # Optional weight
    n: Optional[int] = None  # Sample size (optional)
    
    @property
    def se(self) -> float:
        """Calculate standard error from confidence interval."""
        # Assuming 95% CI: CI = effect ± 1.96*SE
        return (self.upper - self.lower) / (2 * 1.96)
    
    @property
    def variance(self) -> float:
        """Calculate variance."""
        return self.se ** 2
    
    @property
    def precision(self) -> float:
        """Calculate precision (1/variance)."""
        return 1 / self.variance


@dataclass
class MetaAnalysisResult:
    """Results from meta-analysis."""
    # Fixed-effect model
    fe_effect: float
    fe_lower: float
    fe_upper: float
    fe_se: float
    fe_pvalue: float
    
    # Random-effects model
    re_effect: float
    re_lower: float
    re_upper: float
    re_se: float
    re_pvalue: float
    
    # Heterogeneity statistics
    Q: float          # Cochran's Q statistic
    I2: float         # I-squared (%)
    tau2: float       # Between-study variance
    Q_pvalue: float   # P-value for Q statistic
    
    # Study weights
    fe_weights: List[float]
    re_weights: List[float]
    
    # Input studies
    studies: List[Study]


class MetaAnalyzer:
    """Perform meta-analysis calculations."""
    
    @staticmethod
    def fixed_effect(studies: List[Study]) -> Tuple[float, float, float, List[float]]:
        """
        Fixed-effect model using inverse variance weighting.
        
        Returns:
            overall_effect, se, z, weights
        """
        # Calculate weights (inverse variance)
        weights = [1 / s.variance for s in studies]
        total_weight = sum(weights)
        
        # Calculate weighted mean
        weighted_sum = sum(s.effect * w for s, w in zip(studies, weights))
        overall_effect = weighted_sum / total_weight
        
        # Calculate SE of overall effect
        se = np.sqrt(1 / total_weight)
        
        # Calculate z-score
        z = overall_effect / se
        
        # Normalize weights to sum to 100%
        normalized_weights = [100 * w / total_weight for w in weights]
        
        return overall_effect, se, z, normalized_weights
    
    @staticmethod
    def calculate_heterogeneity(studies: List[Study], fe_effect: float) -> Tuple[float, float, float, float]:
        """
        Calculate heterogeneity statistics.
        
        Returns:
            Q, I2, tau2, Q_pvalue
        """
        k = len(studies)
        
        # Cochran's Q statistic
        Q = sum((s.effect - fe_effect) ** 2 / s.variance for s in studies)
        
        # Degrees of freedom
        df = k - 1
        
        # P-value for Q (chi-square distribution)
        Q_pvalue = 1 - stats.chi2.cdf(Q, df)
        
        # I-squared (percentage of variation due to heterogeneity)
        if Q <= df:
            I2 = 0
        else:
            I2 = max(0, 100 * (Q - df) / Q)
        
        # Tau-squared (between-study variance)
        # DerSimonian-Laird estimator
        weights = [1 / s.variance for s in studies]
        total_weight = sum(weights)
        sum_w_squared = sum(w ** 2 for w in weights)
        sum_w = sum(w for w in weights)
        
        # Calculate U
        U = sum_w - sum_w_squared / sum_w
        
        if U > 0 and Q > df:
            tau2 = max(0, (Q - df) / U)
        else:
            tau2 = 0
        
        return Q, I2, tau2, Q_pvalue
    
    @staticmethod
    def random_effects(studies: List[Study], tau2: float) -> Tuple[float, float, float, List[float]]:
        """
        Random-effects model (DerSimonian-Laird method).
        
        Returns:
            overall_effect, se, z, weights
        """
        # Calculate random-effects weights
        weights = [1 / (s.variance + tau2) for s in studies]
        total_weight = sum(weights)
        
        # Calculate weighted mean
        weighted_sum = sum(s.effect * w for s, w in zip(studies, weights))
        overall_effect = weighted_sum / total_weight
        
        # Calculate SE
        se = np.sqrt(1 / total_weight)
        
        # Calculate z-score
        z = overall_effect / se
        
        # Normalize weights
        normalized_weights = [100 * w / total_weight for w in weights]
        
        return overall_effect, se, z, normalized_weights
    
    @classmethod
    def analyze(cls, studies: List[Study]) -> MetaAnalysisResult:
        """Perform complete meta-analysis."""
        if len(studies) < 2:
            raise ValueError("At least 2 studies required for meta-analysis")
        
        # Fixed-effect model
        fe_effect, fe_se, fe_z, fe_weights = cls.fixed_effect(studies)
        fe_lower = fe_effect - 1.96 * fe_se
        fe_upper = fe_effect + 1.96 * fe_se
        fe_pvalue = 2 * (1 - stats.norm.cdf(abs(fe_z)))
        
        # Heterogeneity
        Q, I2, tau2, Q_pvalue = cls.calculate_heterogeneity(studies, fe_effect)
        
        # Random-effects model
        re_effect, re_se, re_z, re_weights = cls.random_effects(studies, tau2)
        re_lower = re_effect - 1.96 * re_se
        re_upper = re_effect + 1.96 * re_se
        re_pvalue = 2 * (1 - stats.norm.cdf(abs(re_z)))
        
        return MetaAnalysisResult(
            fe_effect=fe_effect,
            fe_lower=fe_lower,
            fe_upper=fe_upper,
            fe_se=fe_se,
            fe_pvalue=fe_pvalue,
            re_effect=re_effect,
            re_lower=re_lower,
            re_upper=re_upper,
            re_se=re_se,
            re_pvalue=re_pvalue,
            Q=Q,
            I2=I2,
            tau2=tau2,
            Q_pvalue=Q_pvalue,
            fe_weights=fe_weights,
            re_weights=re_weights,
            studies=studies
        )


class ForestPlotter:
    """Generate forest plots."""
    
    def __init__(self, result: MetaAnalysisResult):
        self.result = result
    
    def plot(self, 
             output_file: str = "forest_plot.png",
             title: str = "Meta-analysis Forest Plot",
             effect_label: str = "Effect Size",
             null_value: float = 0,
             show_weights: bool = True,
             show_model: str = "both") -> None:
        """
        Generate forest plot.
        
        Args:
            output_file: Output filename (png/pdf/svg)
            title: Plot title
            effect_label: Label for x-axis
            null_value: Value representing no effect (0 for MD, 1 for log(OR))
            show_weights: Whether to show study weights
            show_model: "fixed", "random", or "both"
        """
        if not MATPLOTLIB_AVAILABLE:
            raise ImportError("matplotlib is required for plotting")
        
        studies = self.result.studies
        n_studies = len(studies)
        
        # Calculate figure size
        fig_height = max(6, n_studies * 0.5 + 4)
        
        # Create figure with multiple subplots for alignment
        fig = plt.figure(figsize=(12, fig_height))
        
        # Create grid spec
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(1, 3, width_ratios=[3, 1, 4], wspace=0.05)
        
        # Left panel: Study names
        ax_left = fig.add_subplot(gs[0])
        ax_left.axis('off')
        
        # Middle panel: Weights (optional)
        if show_weights:
            ax_mid = fig.add_subplot(gs[1])
            ax_mid.axis('off')
        
        # Right panel: Forest plot
        ax_right = fig.add_subplot(gs[2])
        
        # Y positions (reverse order so Study 1 is at top)
        y_pos = np.arange(n_studies + 2)[::-1]  # +2 for models
        
        # Determine which model to show
        if show_model == "fixed":
            summary_effect = self.result.fe_effect
            summary_lower = self.result.fe_lower
            summary_upper = self.result.fe_upper
            summary_weights = self.result.fe_weights
        elif show_model == "random":
            summary_effect = self.result.re_effect
            summary_lower = self.result.re_lower
            summary_upper = self.result.re_upper
            summary_weights = self.result.re_weights
        else:  # both - use random effects as primary
            summary_effect = self.result.re_effect
            summary_lower = self.result.re_lower
            summary_upper = self.result.re_upper
            summary_weights = self.result.re_weights
        
        # Determine x-axis limits
        all_effects = [s.effect for s in studies] + [summary_effect]
        all_lowers = [s.lower for s in studies] + [summary_lower]
        all_uppers = [s.upper for s in studies] + [summary_upper]
        
        x_min = min(all_lowers) - 0.1 * abs(min(all_lowers))
        x_max = max(all_uppers) + 0.1 * abs(max(all_uppers))
        
        # Add margin
        margin = (x_max - x_min) * 0.1
        x_min -= margin
        x_max += margin
        
        # Plot individual studies
        for i, study in enumerate(studies):
            y = y_pos[i + 1]  # Skip top position
            
            # Plot confidence interval
            ax_right.plot([study.lower, study.upper], [y, y], 
                         'k-', linewidth=1.5, alpha=0.7)
            
            # Plot effect size marker
            ax_right.plot(study.effect, y, 'ks', markersize=8)
            
            # Add study name to left panel
            ax_left.text(1, y, study.name, 
                        ha='right', va='center', fontsize=10)
            
            # Add weight to middle panel
            if show_weights:
                weight_str = f"{summary_weights[i]:.1f}%"
                ax_mid.text(0.5, y, weight_str,
                           ha='center', va='center', fontsize=9)
        
        # Plot summary effect (diamond)
        y_summary = y_pos[-1]  # Bottom position
        
        # Create diamond
        diamond = plt.Polygon([
            [summary_lower, y_summary],
            [summary_effect, y_summary + 0.3],
            [summary_upper, y_summary],
            [summary_effect, y_summary - 0.3]
        ], fill=True, facecolor='red', edgecolor='darkred', linewidth=2)
        ax_right.add_patch(diamond)
        
        # Add summary label
        if show_model == "both":
            summary_label = "Random-effects model"
        elif show_model == "fixed":
            summary_label = "Fixed-effect model"
        else:
            summary_label = "Random-effects model"
        
        ax_left.text(1, y_summary, summary_label,
                    ha='right', va='center', fontsize=10, fontweight='bold')
        
        if show_weights:
            ax_mid.text(0.5, y_summary, "100%",
                       ha='center', va='center', fontsize=9, fontweight='bold')
        
        # Add separator line
        ax_right.axhline(y=y_pos[0] - 0.5, color='gray', linestyle='-', linewidth=0.5)
        
        # Configure right panel
        ax_right.set_ylim(-0.5, n_studies + 2.5)
        ax_right.set_xlim(x_min, x_max)
        ax_right.set_yticks([])
        ax_right.set_xlabel(effect_label, fontsize=11)
        ax_right.grid(axis='x', alpha=0.3)
        
        # Add vertical line at null value
        ax_right.axvline(x=null_value, color='red', linestyle='--', 
                        linewidth=1, alpha=0.5)
        
        # Add effect size values and CIs on the right
        for i, study in enumerate(studies):
            y = y_pos[i + 1]
            text = f"{study.effect:.2f} [{study.lower:.2f}, {study.upper:.2f}]"
            ax_right.text(x_max + margin * 0.1, y, text,
                         ha='left', va='center', fontsize=9)
        
        # Add summary effect text
        summary_text = f"{summary_effect:.2f} [{summary_lower:.2f}, {summary_upper:.2f}]"
        ax_right.text(x_max + margin * 0.1, y_summary, summary_text,
                     ha='left', va='center', fontsize=9, fontweight='bold')
        
        # Configure middle panel
        if show_weights:
            ax_mid.set_ylim(-0.5, n_studies + 2.5)
            ax_mid.set_xlim(0, 1)
            ax_mid.text(0.5, y_pos[0], "Weight",
                       ha='center', va='center', fontsize=9, fontweight='bold')
        
        # Configure left panel
        ax_left.set_ylim(-0.5, n_studies + 2.5)
        ax_left.set_xlim(0, 1)
        ax_left.text(1, y_pos[0], "Study",
                    ha='right', va='center', fontsize=10, fontweight='bold')
        
        # Add title
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
        
        # Add heterogeneity statistics as text box
        hetero_text = (
            f"Heterogeneity: Q={self.result.Q:.2f}, "
            f"I²={self.result.I2:.1f}%, "
            f"τ²={self.result.tau2:.3f}"
        )
        fig.text(0.5, 0.02, hetero_text, ha='center', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Add model statistics
        if show_model in ["random", "both"]:
            model_text = (
                f"Random-effects: ES={self.result.re_effect:.2f} "
                f"(95% CI: {self.result.re_lower:.2f} to {self.result.re_upper:.2f}), "
                f"p={self.result.re_pvalue:.3f}"
            )
        else:
            model_text = (
                f"Fixed-effect: ES={self.result.fe_effect:.2f} "
                f"(95% CI: {self.result.fe_lower:.2f} to {self.result.fe_upper:.2f}), "
                f"p={self.result.fe_pvalue:.3f}"
            )
        fig.text(0.5, 0.06, model_text, ha='center', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
        
        plt.tight_layout(rect=[0, 0.08, 1, 0.96])
        
        # Save figure
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Forest plot saved to: {output_file}")


def read_studies_from_csv(filepath: str) -> List[Study]:
    """Read studies from CSV file."""
    studies = []
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            study = Study(
                name=row.get('name', row.get('study', 'Unknown')),
                effect=float(row['effect']),
                lower=float(row['lower']),
                upper=float(row['upper']),
                weight=float(row.get('weight', 0)) if row.get('weight') else None,
                n=int(row.get('n', 0)) if row.get('n') else None
            )
            studies.append(study)
    return studies


def read_studies_from_json(filepath: str) -> List[Study]:
    """Read studies from JSON file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    studies = []
    for item in data:
        study = Study(
            name=item.get('name', item.get('study', 'Unknown')),
            effect=float(item['effect']),
            lower=float(item['lower']),
            upper=float(item['upper']),
            weight=float(item.get('weight', 0)) if item.get('weight') else None,
            n=int(item.get('n', 0)) if item.get('n') else None
        )
        studies.append(study)
    return studies


def print_results(result: MetaAnalysisResult):
    """Print meta-analysis results to console."""
    print("\n" + "="*70)
    print("META-ANALYSIS RESULTS")
    print("="*70)
    
    print(f"\nNumber of studies: {len(result.studies)}")
    
    print("\n" + "-"*70)
    print("STUDY-LEVEL DATA")
    print("-"*70)
    print(f"{'Study':<20} {'Effect':>10} {'95% CI':>20} {'Weight':>10}")
    print("-"*70)
    
    for i, study in enumerate(result.studies):
        ci = f"[{study.lower:.2f}, {study.upper:.2f}]"
        weight = f"{result.re_weights[i]:.1f}%"
        print(f"{study.name:<20} {study.effect:>10.2f} {ci:>20} {weight:>10}")
    
    print("\n" + "-"*70)
    print("SUMMARY EFFECTS")
    print("-"*70)
    
    print(f"\nFixed-effect model:")
    print(f"  Effect size: {result.fe_effect:.3f}")
    print(f"  95% CI: [{result.fe_lower:.3f}, {result.fe_upper:.3f}]")
    print(f"  SE: {result.fe_se:.3f}")
    print(f"  Z: {result.fe_effect/result.fe_se:.3f}")
    print(f"  P-value: {result.fe_pvalue:.4f}")
    
    print(f"\nRandom-effects model:")
    print(f"  Effect size: {result.re_effect:.3f}")
    print(f"  95% CI: [{result.re_lower:.3f}, {result.re_upper:.3f}]")
    print(f"  SE: {result.re_se:.3f}")
    print(f"  Z: {result.re_effect/result.re_se:.3f}")
    print(f"  P-value: {result.re_pvalue:.4f}")
    
    print("\n" + "-"*70)
    print("HETEROGENEITY STATISTICS")
    print("-"*70)
    print(f"  Cochran's Q: {result.Q:.3f}")
    print(f"  Degrees of freedom: {len(result.studies) - 1}")
    print(f"  P-value: {result.Q_pvalue:.4f}")
    print(f"  I²: {result.I2:.1f}%")
    print(f"  τ²: {result.tau2:.4f}")
    
    interpretation = []
    if result.I2 < 25:
        interpretation.append("Low heterogeneity")
    elif result.I2 < 50:
        interpretation.append("Moderate heterogeneity")
    elif result.I2 < 75:
        interpretation.append("Substantial heterogeneity")
    else:
        interpretation.append("Considerable heterogeneity")
    
    if result.Q_pvalue < 0.05:
        interpretation.append("(significant)")
    else:
        interpretation.append("(not significant)")
    
    print(f"  Interpretation: {' '.join(interpretation)}")
    print("="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Meta-analysis Forest Plotter (Pure Python)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode
  python main.py --demo
  
  # Analyze from CSV file
  python main.py --input studies.csv --output forest.png
  
  # Analyze from JSON file
  python main.py --input studies.json --output forest.pdf
  
  # Custom settings
  python main.py --input studies.csv --model random --title "My Meta-analysis"
  
CSV Format:
  name,effect,lower,upper,weight,n
  Study 1,0.5,0.3,0.7,25,100
  Study 2,0.6,0.4,0.8,30,120
        """
    )
    
    parser.add_argument("--input", "-i", type=str,
                       help="Input file (CSV or JSON)")
    parser.add_argument("--output", "-o", type=str, default="forest_plot.png",
                       help="Output file (png/pdf/svg)")
    parser.add_argument("--title", "-t", type=str, default="Meta-analysis Forest Plot",
                       help="Plot title")
    parser.add_argument("--model", "-m", type=str, default="both",
                       choices=["fixed", "random", "both"],
                       help="Model to display (default: both)")
    parser.add_argument("--effect-label", type=str, default="Effect Size",
                       help="Label for effect size axis")
    parser.add_argument("--null-value", type=float, default=0,
                       help="Null effect value (default: 0)")
    parser.add_argument("--no-weights", action="store_true",
                       help="Hide weights column")
    parser.add_argument("--demo", action="store_true",
                       help="Run demo with example data")
    
    args = parser.parse_args()
    
    # Check matplotlib
    if not MATPLOTLIB_AVAILABLE:
        print("❌ Error: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)
    
    # Load studies
    if args.demo:
        print("🎮 Running DEMO mode with example data...")
        studies = [
            Study(name="Smith et al. 2020", effect=0.45, lower=0.20, upper=0.70, n=120),
            Study(name="Jones et al. 2019", effect=0.62, lower=0.35, upper=0.89, n=150),
            Study(name="Wang et al. 2021", effect=0.38, lower=0.12, upper=0.64, n=98),
            Study(name="Chen et al. 2018", effect=0.71, lower=0.48, upper=0.94, n=175),
            Study(name="Johnson et al. 2022", effect=0.55, lower=0.28, upper=0.82, n=134),
        ]
    elif args.input:
        print(f"📁 Loading studies from: {args.input}")
        if args.input.endswith('.json'):
            studies = read_studies_from_json(args.input)
        elif args.input.endswith('.csv'):
            studies = read_studies_from_csv(args.input)
        else:
            print("❌ Error: Input file must be .csv or .json")
            sys.exit(1)
        print(f"✓ Loaded {len(studies)} studies")
    else:
        parser.print_help()
        print("\n❌ Error: Please provide --input file or use --demo")
        sys.exit(1)
    
    # Perform meta-analysis
    print("\n🔍 Performing meta-analysis...")
    result = MetaAnalyzer.analyze(studies)
    
    # Print results
    print_results(result)
    
    # Generate forest plot
    print("🎨 Generating forest plot...")
    plotter = ForestPlotter(result)
    plotter.plot(
        output_file=args.output,
        title=args.title,
        effect_label=args.effect_label,
        null_value=args.null_value,
        show_weights=not args.no_weights,
        show_model=args.model
    )
    
    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\n✓ Results summary:")
    print(f"  Studies analyzed: {len(studies)}")
    print(f"  Random-effects ES: {result.re_effect:.3f} (95% CI: {result.re_lower:.3f}-{result.re_upper:.3f})")
    print(f"  Heterogeneity I²: {result.I2:.1f}%")
    print(f"\n✓ Output saved to: {args.output}")


if __name__ == "__main__":
    main()
