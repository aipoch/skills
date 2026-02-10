#!/usr/bin/env python3
"""
Mouse Colony Planner - Transgenic Mouse Breeding Planning Tool

Features:
- Calculate breeding timelines
- Estimate cage requirements
- Predict breeding costs
"""

import argparse
import math
from dataclasses import dataclass
from typing import List, Dict, Tuple
from enum import Enum


class BreedingScheme(Enum):
    """Breeding scheme types"""
    HETEROZYGOTE = "heterozygote"  # Heterozygote x Wild type
    HOMOZYGOTE = "homozygote"      # Heterozygote x Heterozygote
    CONDITIONAL = "conditional"    # Conditional knockout (Cre/loxp)


@dataclass
class BreedingParams:
    """Breeding parameters"""
    gestation_days: int = 21       # Gestation period
    weaning_days: int = 21         # Weaning age
    sexual_maturity_days: int = 42 # Sexual maturity age
    litter_size: int = 8           # Average litter size
    female_puberty: int = 35       # Female puberty (breeding age)
    male_puberty: int = 35         # Male puberty (breeding age)
    
    # Genotype ratios
    het_ratio: float = 0.5         # Heterozygote ratio
    homo_ratio: float = 0.25       # Homozygote ratio (het x het)
    

@dataclass
class CostParams:
    """Cost parameters"""
    cage_cost_per_day: float = 3.0      # Cost per cage per day (USD)
    genotyping_cost: float = 15.0       # Genotyping cost per mouse (USD)
    mouse_purchase_cost: float = 50.0   # Mouse purchase cost (USD)


@dataclass
class Phase:
    """Breeding phase"""
    name: str
    duration_days: int
    cages_needed: int
    description: str


@dataclass
class ColonyPlan:
    """Breeding plan results"""
    scheme: BreedingScheme
    phases: List[Phase]
    total_days: int
    total_cages: int
    total_cost: float
    breeding_pairs: int
    expected_target_genotype: int


def calculate_breeding_plan(
    scheme: BreedingScheme,
    initial_females: int,
    initial_males: int,
    target_pups: int,
    breeding_params: BreedingParams = None,
    cost_params: CostParams = None,
    cage_capacity: int = 5
) -> ColonyPlan:
    """
    Calculate breeding plan
    
    Args:
        scheme: Breeding scheme
        initial_females: Starting number of females
        initial_males: Starting number of males
        target_pups: Target number of specific genotype mice
        breeding_params: Breeding parameters
        cost_params: Cost parameters
        cage_capacity: Maximum cage capacity
    
    Returns:
        ColonyPlan: Breeding plan
    """
    if breeding_params is None:
        breeding_params = BreedingParams()
    if cost_params is None:
        cost_params = CostParams()
    
    phases = []
    
    # Calculate required breeding pairs
    if scheme == BreedingScheme.HETEROZYGOTE:
        # Heterozygote x Wild type → 50% heterozygotes
        pups_per_pair = breeding_params.litter_size * breeding_params.het_ratio
        breeding_pairs_needed = math.ceil(target_pups / pups_per_pair)
        
    elif scheme == BreedingScheme.HOMOZYGOTE:
        # Heterozygote x Heterozygote → 25% homozygotes
        pups_per_pair = breeding_params.litter_size * breeding_params.homo_ratio
        breeding_pairs_needed = math.ceil(target_pups / pups_per_pair)
        
    else:  # CONDITIONAL
        # Two-step breeding: first get flox/+, then mate with Cre
        # Simplified to require more breeding pairs
        pups_per_pair = breeding_params.litter_size * breeding_params.homo_ratio * 0.5
        breeding_pairs_needed = math.ceil(target_pups / pups_per_pair)
    
    # Ensure sufficient mice
    breeding_pairs = min(
        initial_females,
        initial_males,
        breeding_pairs_needed
    )
    
    # Phase 1: Acclimation (3-7 days)
    adapt_duration = 7
    adapt_cages = math.ceil((initial_females + initial_males) / cage_capacity)
    phases.append(Phase(
        name="Acclimation",
        duration_days=adapt_duration,
        cages_needed=adapt_cages,
        description="New mice acclimate to environment, health observation"
    ))
    
    # Phase 2: Breeding (gestation + lactation)
    # One litter time = gestation + weaning = 21 + 21 = 42 days
    breed_duration = breeding_params.gestation_days + breeding_params.weaning_days
    breed_cages = math.ceil((breeding_pairs * 2 + breeding_pairs * breeding_params.litter_size) / cage_capacity)
    
    # May need multiple litters
    litters_needed = math.ceil(breeding_pairs_needed / breeding_pairs)
    actual_breed_duration = breed_duration * litters_needed
    
    phases.append(Phase(
        name="Breeding Phase",
        duration_days=actual_breed_duration,
        cages_needed=max(breed_cages, adapt_cages),
        description=f"Mating, gestation, delivery, lactation (estimated {litters_needed} litters)"
    ))
    
    # Phase 3: Post-weaning housing
    wean_duration = 21  # Post-weaning to genotyping/separation
    
    if scheme == BreedingScheme.CONDITIONAL:
        # Conditional knockout requires additional step
        wean_duration += breeding_params.sexual_maturity_days  # Need to reach sexual maturity before mating
        
        # Second phase breeding cages
        cond_breed_cages = math.ceil(target_pups / cage_capacity)
        phases.append(Phase(
            name="Conditional Knockout Phase 2",
            duration_days=breeding_params.gestation_days + breeding_params.weaning_days,
            cages_needed=cond_breed_cages,
            description="Flox mice mated with Cre driver to obtain conditional knockout mice"
        ))
    
    # Post-weaning separation
    pups_per_litter = breeding_pairs * breeding_params.litter_size * litters_needed
    wean_cages = math.ceil(pups_per_litter / cage_capacity)
    phases.append(Phase(
        name="Post-weaning Housing",
        duration_days=wean_duration,
        cages_needed=wean_cages,
        description="Post-weaning separation, preparation for genotyping"
    ))
    
    # Phase 4: Genotyping
    geno_duration = 3  # Sampling + testing time
    geno_cages = wean_cages  # Still need cages during genotyping
    phases.append(Phase(
        name="Genotyping",
        duration_days=geno_duration,
        cages_needed=geno_cages,
        description="DNA extraction and PCR genotyping"
    ))
    
    # Calculate total time and max cages
    total_days = sum(p.duration_days for p in phases)
    max_cages = max(p.cages_needed for p in phases)
    
    # Calculate costs
    # Cage cost = sum of (cage days per phase) × unit price
    cage_days = sum(p.duration_days * p.cages_needed for p in phases)
    cage_cost = cage_days * cost_params.cage_cost_per_day
    
    # Genotyping cost
    total_pups = breeding_pairs * breeding_params.litter_size * litters_needed
    if scheme == BreedingScheme.CONDITIONAL:
        genotyping_cost = total_pups * cost_params.genotyping_cost * 2  # Need two rounds of genotyping
    else:
        genotyping_cost = total_pups * cost_params.genotyping_cost
    
    # Mouse purchase cost (initial mice)
    purchase_cost = (initial_females + initial_males) * cost_params.mouse_purchase_cost
    
    total_cost = cage_cost + genotyping_cost + purchase_cost
    
    # Calculate expected target genotype mice
    if scheme == BreedingScheme.HETEROZYGOTE:
        expected_target = int(total_pups * breeding_params.het_ratio)
    elif scheme == BreedingScheme.HOMOZYGOTE:
        expected_target = int(total_pups * breeding_params.homo_ratio)
    else:
        expected_target = int(total_pups * breeding_params.homo_ratio * 0.5)
    
    return ColonyPlan(
        scheme=scheme,
        phases=phases,
        total_days=total_days,
        total_cages=max_cages,
        total_cost=total_cost,
        breeding_pairs=breeding_pairs,
        expected_target_genotype=expected_target
    )


def format_plan_output(plan: ColonyPlan) -> str:
    """Format breeding plan output"""
    lines = []
    lines.append("=" * 60)
    lines.append(f"🐭 Transgenic Mouse Breeding Plan - {plan.scheme.value}")
    lines.append("=" * 60)
    lines.append("")
    
    lines.append("📋 Breeding Phases:")
    lines.append("-" * 40)
    for i, phase in enumerate(plan.phases, 1):
        lines.append(f"  Phase {i}: {phase.name}")
        lines.append(f"    Duration: {phase.duration_days} days")
        lines.append(f"    Cages needed: {phase.cages_needed}")
        lines.append(f"    Description: {phase.description}")
        lines.append("")
    
    lines.append("-" * 40)
    lines.append(f"⏱️  Estimated Total Time: {plan.total_days} days ({plan.total_days/30:.1f} months)")
    lines.append(f"🏠 Maximum Cages: {plan.total_cages}")
    lines.append(f"🧬 Breeding Pairs: {plan.breeding_pairs}")
    lines.append("")
    
    lines.append("💰 Cost Estimate:")
    lines.append("-" * 40)
    lines.append(f"  Expected target genotype mice: {plan.expected_target_genotype}")
    lines.append(f"  Cost per mouse: ${plan.total_cost/max(plan.expected_target_genotype, 1):.1f}")
    lines.append(f"  💵 Total Cost: ${plan.total_cost:.2f}")
    lines.append("")
    
    lines.append("=" * 60)
    
    return "\n".join(lines)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Mouse Colony Planner - Transgenic Mouse Breeding Planning Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --scheme heterozygote --females 10 --males 5 --target-pups 10
  python main.py --scheme homozygote --females 20 --males 10 --target-pups 20
  python main.py --scheme conditional --females 15 --males 15 --target-pups 15
        """
    )
    
    parser.add_argument(
        "--scheme",
        type=str,
        required=True,
        choices=[s.value for s in BreedingScheme],
        help="Breeding scheme type"
    )
    
    parser.add_argument(
        "--females",
        type=int,
        required=True,
        help="Starting number of females"
    )
    
    parser.add_argument(
        "--males",
        type=int,
        required=True,
        help="Starting number of males"
    )
    
    parser.add_argument(
        "--target-pups",
        type=int,
        default=10,
        help="Target number of specific genotype mice (default: 10)"
    )
    
    parser.add_argument(
        "--gestation",
        type=int,
        default=21,
        help="Gestation period (days) (default: 21)"
    )
    
    parser.add_argument(
        "--weaning",
        type=int,
        default=21,
        help="Weaning age (days) (default: 21)"
    )
    
    parser.add_argument(
        "--sexual-maturity",
        type=int,
        default=42,
        help="Sexual maturity age (days) (default: 42)"
    )
    
    parser.add_argument(
        "--cage-capacity",
        type=int,
        default=5,
        help="Maximum cage capacity (default: 5)"
    )
    
    parser.add_argument(
        "--cage-cost",
        type=float,
        default=3.0,
        help="Cost per cage per day (USD) (default: 3.0)"
    )
    
    parser.add_argument(
        "--genotyping-cost",
        type=float,
        default=15.0,
        help="Genotyping cost per mouse (USD) (default: 15.0)"
    )
    
    return parser.parse_args()


def main():
    """Main function"""
    args = parse_arguments()
    
    # Create breeding parameters
    breeding_params = BreedingParams(
        gestation_days=args.gestation,
        weaning_days=args.weaning,
        sexual_maturity_days=args.sexual_maturity
    )
    
    # Create cost parameters
    cost_params = CostParams(
        cage_cost_per_day=args.cage_cost,
        genotyping_cost=args.genotyping_cost
    )
    
    # Parse breeding scheme
    scheme = BreedingScheme(args.scheme)
    
    # Calculate breeding plan
    plan = calculate_breeding_plan(
        scheme=scheme,
        initial_females=args.females,
        initial_males=args.males,
        target_pups=args.target_pups,
        breeding_params=breeding_params,
        cost_params=cost_params,
        cage_capacity=args.cage_capacity
    )
    
    # Output results
    print(format_plan_output(plan))


if __name__ == "__main__":
    main()
