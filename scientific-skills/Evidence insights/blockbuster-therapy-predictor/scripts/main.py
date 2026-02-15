#!/usr/bin/env python3
"""
Blockbuster Therapy Predictor
Predicts the potential of early-stage technology platforms to become blockbuster therapies.

Data dimensions: Clinical trials + Patent landscape + VC funding
"""

import json
import argparse
from datetime import datetime
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
from enum import Enum
import random


class InvestmentRecommendation(Enum):
    """Investment recommendation levels"""
    STRONG_BUY = "Strongly Recommended"
    BUY = "Recommended"
    HOLD = "Watch"
    CAUTION = "Cautious"


@dataclass
class TechnologyRoute:
    """Technology route data model"""
    name: str
    category: str
    description: str
    
    # Clinical trial data
    clinical_trials_phase1: int = 0
    clinical_trials_phase2: int = 0
    clinical_trials_phase3: int = 0
    clinical_success_rate: float = 0.0
    avg_time_to_market: float = 0.0  # years
    indications_count: int = 0
    
    # Patent data
    patent_count: int = 0
    patent_growth_rate: float = 0.0  # annual growth rate
    core_patents: int = 0
    geographic_coverage: int = 0  # number of countries covered
    
    # Funding data
    total_funding_usd: float = 0.0  # millions USD
    funding_rounds: int = 0
    top_vc_backed: bool = False
    last_valuation_usd: float = 0.0  # millions USD
    companies_count: int = 0


@dataclass
class PredictionResult:
    """Prediction result model"""
    tech_name: str
    maturity_score: float
    market_potential_score: float
    momentum_score: float
    blockbuster_index: float
    recommendation: str
    key_drivers: List[str]
    risk_factors: List[str]
    timeline_prediction: str


class ClinicalDataAnalyzer:
    """Clinical trial data analyzer"""
    
    @staticmethod
    def calculate_clinical_score(tech: TechnologyRoute) -> float:
        """Calculate clinical stage score (0-100)"""
        # Trial phase weights
        phase_weights = {1: 0.2, 2: 0.5, 3: 0.8}
        total_trials = (tech.clinical_trials_phase1 + 
                       tech.clinical_trials_phase2 + 
                       tech.clinical_trials_phase3)
        
        if total_trials == 0:
            return 10.0  # base score
        
        weighted_score = (
            tech.clinical_trials_phase1 * phase_weights[1] +
            tech.clinical_trials_phase2 * phase_weights[2] +
            tech.clinical_trials_phase3 * phase_weights[3]
        ) / total_trials * 100
        
        # Success rate adjustment
        weighted_score *= (0.5 + tech.clinical_success_rate)
        
        # Indication diversity bonus
        indication_bonus = min(tech.indications_count * 2, 15)
        
        return min(weighted_score + indication_bonus, 100)


class PatentAnalyzer:
    """Patent landscape analyzer"""
    
    @staticmethod
    def calculate_patent_depth_score(tech: TechnologyRoute) -> float:
        """Calculate patent depth score (0-100)"""
        if tech.patent_count == 0:
            return 5.0
        
        # Base patent quantity score
        patent_base = min(tech.patent_count / 10, 40)
        
        # Core patent quality score
        core_quality = min(tech.core_patents * 3, 30)
        
        # Growth rate score
        growth_score = min(tech.patent_growth_rate * 2, 20)
        
        # Geographic coverage score
        geo_score = min(tech.geographic_coverage * 2, 10)
        
        return min(patent_base + core_quality + growth_score + geo_score, 100)


class FundingAnalyzer:
    """Funding data analyzer"""
    
    @staticmethod
    def calculate_funding_score(tech: TechnologyRoute) -> float:
        """Calculate funding stage score (0-100)"""
        if tech.total_funding_usd == 0:
            return 5.0
        
        # Funding scale score
        funding_score = min(tech.total_funding_usd / 50, 40)
        
        # Round maturity score
        round_score = min(tech.funding_rounds * 8, 25)
        
        # Top-tier VC endorsement score
        vc_bonus = 20 if tech.top_vc_backed else 0
        
        # Company count score (ecosystem activity)
        ecosystem_score = min(tech.companies_count * 3, 15)
        
        return min(funding_score + round_score + vc_bonus + ecosystem_score, 100)


class MarketPotentialEvaluator:
    """Market potential evaluator"""
    
    MARKET_SIZE_ESTIMATES = {
        "PROTAC": 35.0,  # 2030 estimate in billions USD
        "mRNA": 45.0,
        "CRISPR": 25.0,
        "CAR-T": 20.0,
        "Bispecific": 30.0,
        "ADC": 28.0,
        "Cell Therapy": 22.0,
        "Gene Therapy": 18.0,
        "RNAi": 15.0,
        "Allogeneic": 12.0,
    }
    
    UNMET_NEED_SCORES = {
        "PROTAC": 85,  # breakthrough for undruggable targets
        "mRNA": 90,    # vaccines + tumor immunology
        "CRISPR": 88,  # genetic disease cure
        "CAR-T": 75,   # validated in hematological tumors, solid tumors pending
        "Bispecific": 80,
        "ADC": 78,
        "Cell Therapy": 72,
        "Gene Therapy": 85,
        "RNAi": 70,
        "Allogeneic": 82,
    }
    
    COMPETITIVE_LANDSCORE = {
        "PROTAC": 75,   # moderate competition, large differentiation space
        "mRNA": 65,     # intense competition but large market
        "CRISPR": 70,   # high technical barrier
        "CAR-T": 60,    # intense competition
        "Bispecific": 65,
        "ADC": 70,
        "Cell Therapy": 68,
        "Gene Therapy": 72,
        "RNAi": 75,
        "Allogeneic": 70,
    }
    
    @classmethod
    def calculate_market_potential(cls, tech_name: str) -> float:
        """Calculate market potential score (0-100)"""
        market_size = cls.MARKET_SIZE_ESTIMATES.get(tech_name, 10.0)
        unmet_need = cls.UNMET_NEED_SCORES.get(tech_name, 60)
        competitive = cls.COMPETITIVE_LANDSCORE.get(tech_name, 60)
        
        # Market size standardization (max 50 points)
        size_score = min(market_size / 50 * 50, 50)
        
        # Unmet need (35 points)
        need_score = unmet_need * 0.35
        
        # Competitive landscape (15 points)
        comp_score = competitive * 0.15
        
        return min(size_score + need_score + comp_score, 100)


class BlockbusterPredictor:
    """Blockbuster therapy prediction engine"""
    
    def __init__(self):
        self.clinical_analyzer = ClinicalDataAnalyzer()
        self.patent_analyzer = PatentAnalyzer()
        self.funding_analyzer = FundingAnalyzer()
        self.market_evaluator = MarketPotentialEvaluator()
    
    def calculate_maturity_score(self, tech: TechnologyRoute) -> float:
        """Calculate technology maturity score"""
        clinical = self.clinical_analyzer.calculate_clinical_score(tech)
        patent = self.patent_analyzer.calculate_patent_depth_score(tech)
        funding = self.funding_analyzer.calculate_funding_score(tech)
        
        return clinical * 0.4 + patent * 0.3 + funding * 0.3
    
    def calculate_momentum_score(self, tech: TechnologyRoute) -> float:
        """Calculate development momentum score"""
        factors = []
        
        # Patent growth momentum
        if tech.patent_growth_rate > 30:
            factors.append(25)
        elif tech.patent_growth_rate > 15:
            factors.append(15)
        else:
            factors.append(5)
        
        # Funding activity
        if tech.funding_rounds >= 3:
            factors.append(25)
        elif tech.funding_rounds >= 2:
            factors.append(15)
        else:
            factors.append(5)
        
        # Clinical progress
        if tech.clinical_trials_phase3 > 0:
            factors.append(30)
        elif tech.clinical_trials_phase2 > 2:
            factors.append(20)
        elif tech.clinical_trials_phase2 > 0:
            factors.append(10)
        else:
            factors.append(5)
        
        # Ecosystem activity
        eco_score = min(tech.companies_count * 5, 20)
        factors.append(eco_score)
        
        return sum(factors)
    
    def get_recommendation(self, index: float) -> str:
        """Generate investment recommendation based on index"""
        if index >= 80:
            return InvestmentRecommendation.STRONG_BUY.value
        elif index >= 60:
            return InvestmentRecommendation.BUY.value
        elif index >= 40:
            return InvestmentRecommendation.HOLD.value
        else:
            return InvestmentRecommendation.CAUTION.value
    
    def identify_key_drivers(self, tech: TechnologyRoute, scores: Dict) -> List[str]:
        """Identify key driving factors"""
        drivers = []
        
        if tech.clinical_trials_phase3 > 0:
            drivers.append("Phase III clinical trials ongoing, approaching commercialization")
        elif tech.clinical_trials_phase2 > 3:
            drivers.append("Multiple Phase II clinical trials in progress")
        
        if tech.patent_growth_rate > 25:
            drivers.append("Rapid patent growth, strengthening technology moat")
        
        if tech.top_vc_backed:
            drivers.append("Backed by top-tier VCs with sufficient funding")
        
        if tech.indications_count > 3:
            drivers.append("Multi-indication portfolio with broad market potential")
        
        if tech.core_patents > 5:
            drivers.append("Leading core patent portfolio")
        
        return drivers if drivers else ["Emerging technology platform, worth continuous monitoring"]
    
    def identify_risks(self, tech: TechnologyRoute) -> List[str]:
        """Identify risk factors"""
        risks = []
        
        if tech.clinical_trials_phase3 == 0 and tech.clinical_trials_phase2 == 0:
            risks.append("Early-stage with insufficient clinical validation")
        
        if tech.clinical_success_rate < 0.3:
            risks.append("Low historical success rate")
        
        if tech.patent_count < 10:
            risks.append("Relatively weak patent portfolio")
        
        if tech.total_funding_usd < 100:
            risks.append("Limited funding scale may affect R&D progress")
        
        if tech.avg_time_to_market > 8:
            risks.append("Long development cycle with high uncertainty")
        
        return risks if risks else ["Standard R&D risks"]
    
    def predict_timeline(self, tech: TechnologyRoute) -> str:
        """Predict commercialization timeline"""
        if tech.clinical_trials_phase3 > 0:
            return "First product expected to launch in 2-4 years"
        elif tech.clinical_trials_phase2 > 2:
            return "Expected to enter Phase III within 4-6 years"
        elif tech.clinical_trials_phase2 > 0:
            return "Key clinical data expected in 5-7 years"
        else:
            return "Expected to reach commercialization stage in 7-10 years"
    
    def predict(self, tech: TechnologyRoute) -> PredictionResult:
        """Execute complete prediction"""
        maturity = self.calculate_maturity_score(tech)
        market_potential = self.market_evaluator.calculate_market_potential(tech.name)
        momentum = self.calculate_momentum_score(tech)
        
        # Blockbuster index calculation
        blockbuster_index = (
            market_potential * 0.5 +
            maturity * 0.3 +
            momentum * 0.2
        )
        
        return PredictionResult(
            tech_name=tech.name,
            maturity_score=round(maturity, 1),
            market_potential_score=round(market_potential, 1),
            momentum_score=round(momentum, 1),
            blockbuster_index=round(blockbuster_index, 1),
            recommendation=self.get_recommendation(blockbuster_index),
            key_drivers=self.identify_key_drivers(tech, {}),
            risk_factors=self.identify_risks(tech),
            timeline_prediction=self.predict_timeline(tech)
        )


class DataLoader:
    """Data loader - mock/real data sources"""
    
    @staticmethod
    def load_mock_data() -> List[TechnologyRoute]:
        """Load mock data for demonstration"""
        technologies = [
            TechnologyRoute(
                name="PROTAC",
                category="Protein Degradation",
                description="Proteolysis Targeting Chimera technology",
                clinical_trials_phase1=15,
                clinical_trials_phase2=8,
                clinical_trials_phase3=2,
                clinical_success_rate=0.35,
                avg_time_to_market=6.5,
                indications_count=5,
                patent_count=180,
                patent_growth_rate=35,
                core_patents=25,
                geographic_coverage=12,
                total_funding_usd=2800,
                funding_rounds=4,
                top_vc_backed=True,
                last_valuation_usd=450,
                companies_count=18
            ),
            TechnologyRoute(
                name="mRNA",
                category="Nucleic Acid Drugs",
                description="Messenger RNA therapy platform",
                clinical_trials_phase1=25,
                clinical_trials_phase2=12,
                clinical_trials_phase3=5,
                clinical_success_rate=0.45,
                avg_time_to_market=5.0,
                indications_count=8,
                patent_count=450,
                patent_growth_rate=28,
                core_patents=40,
                geographic_coverage=15,
                total_funding_usd=5200,
                funding_rounds=5,
                top_vc_backed=True,
                last_valuation_usd=800,
                companies_count=25
            ),
            TechnologyRoute(
                name="CRISPR",
                category="Gene Editing",
                description="CRISPR-Cas gene editing technology",
                clinical_trials_phase1=12,
                clinical_trials_phase2=6,
                clinical_trials_phase3=2,
                clinical_success_rate=0.40,
                avg_time_to_market=7.0,
                indications_count=6,
                patent_count=320,
                patent_growth_rate=22,
                core_patents=35,
                geographic_coverage=14,
                total_funding_usd=3800,
                funding_rounds=4,
                top_vc_backed=True,
                last_valuation_usd=600,
                companies_count=15
            ),
            TechnologyRoute(
                name="CAR-T",
                category="Cell Therapy",
                description="Chimeric Antigen Receptor T-cell therapy",
                clinical_trials_phase1=30,
                clinical_trials_phase2=15,
                clinical_trials_phase3=8,
                clinical_success_rate=0.50,
                avg_time_to_market=4.5,
                indications_count=4,
                patent_count=520,
                patent_growth_rate=18,
                core_patents=45,
                geographic_coverage=16,
                total_funding_usd=6500,
                funding_rounds=6,
                top_vc_backed=True,
                last_valuation_usd=1200,
                companies_count=32
            ),
            TechnologyRoute(
                name="Bispecific",
                category="Antibody Drugs",
                description="Bispecific antibody technology",
                clinical_trials_phase1=20,
                clinical_trials_phase2=10,
                clinical_trials_phase3=4,
                clinical_success_rate=0.42,
                avg_time_to_market=5.5,
                indications_count=6,
                patent_count=380,
                patent_growth_rate=20,
                core_patents=30,
                geographic_coverage=14,
                total_funding_usd=3200,
                funding_rounds=5,
                top_vc_backed=True,
                last_valuation_usd=550,
                companies_count=22
            ),
            TechnologyRoute(
                name="ADC",
                category="Antibody Drugs",
                description="Antibody-Drug Conjugate",
                clinical_trials_phase1=22,
                clinical_trials_phase2=12,
                clinical_trials_phase3=6,
                clinical_success_rate=0.48,
                avg_time_to_market=5.0,
                indications_count=5,
                patent_count=420,
                patent_growth_rate=25,
                core_patents=38,
                geographic_coverage=15,
                total_funding_usd=4100,
                funding_rounds=5,
                top_vc_backed=True,
                last_valuation_usd=700,
                companies_count=28
            ),
            TechnologyRoute(
                name="RNAi",
                category="Nucleic Acid Drugs",
                description="RNA interference therapy",
                clinical_trials_phase1=10,
                clinical_trials_phase2=5,
                clinical_trials_phase3=2,
                clinical_success_rate=0.38,
                avg_time_to_market=6.0,
                indications_count=4,
                patent_count=280,
                patent_growth_rate=15,
                core_patents=22,
                geographic_coverage=12,
                total_funding_usd=2100,
                funding_rounds=4,
                top_vc_backed=False,
                last_valuation_usd=350,
                companies_count=12
            ),
            TechnologyRoute(
                name="Gene Therapy",
                category="Gene Therapy",
                description="AAV vector gene therapy",
                clinical_trials_phase1=14,
                clinical_trials_phase2=7,
                clinical_trials_phase3=3,
                clinical_success_rate=0.35,
                avg_time_to_market=6.5,
                indications_count=5,
                patent_count=350,
                patent_growth_rate=18,
                core_patents=28,
                geographic_coverage=13,
                total_funding_usd=2900,
                funding_rounds=4,
                top_vc_backed=True,
                last_valuation_usd=480,
                companies_count=16
            ),
            TechnologyRoute(
                name="Allogeneic",
                category="Cell Therapy",
                description="Universal/Allogeneic cell therapy",
                clinical_trials_phase1=8,
                clinical_trials_phase2=4,
                clinical_trials_phase3=1,
                clinical_success_rate=0.30,
                avg_time_to_market=7.5,
                indications_count=3,
                patent_count=150,
                patent_growth_rate=40,
                core_patents=15,
                geographic_coverage=10,
                total_funding_usd=1200,
                funding_rounds=3,
                top_vc_backed=False,
                last_valuation_usd=220,
                companies_count=10
            ),
            TechnologyRoute(
                name="Cell Therapy",
                category="Cell Therapy",
                description="General cell therapy platform",
                clinical_trials_phase1=12,
                clinical_trials_phase2=6,
                clinical_trials_phase3=2,
                clinical_success_rate=0.33,
                avg_time_to_market=6.8,
                indications_count=4,
                patent_count=220,
                patent_growth_rate=24,
                core_patents=20,
                geographic_coverage=11,
                total_funding_usd=2400,
                funding_rounds=4,
                top_vc_backed=True,
                last_valuation_usd=380,
                companies_count=14
            ),
        ]
        return technologies


class ReportGenerator:
    """Report generator"""
    
    @staticmethod
    def generate_console_report(results: List[PredictionResult], threshold: float = 0):
        """Generate console report"""
        # Filter and sort
        filtered = [r for r in results if r.blockbuster_index >= threshold]
        sorted_results = sorted(filtered, key=lambda x: x.blockbuster_index, reverse=True)
        
        print("\n" + "=" * 100)
        print("🏆 BLOCKBUSTER THERAPY PREDICTOR Report".center(100))
        print("=" * 100)
        print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Technologies analyzed: {len(sorted_results)}")
        print("-" * 100)
        
        # Ranking table
        print("\n📊 Technology Rankings")
        print("-" * 100)
        print(f"{'Rank':<6}{'Technology':<15}{'Blockbuster Index':<20}{'Maturity':<12}{'Market Potential':<18}{'Momentum':<12}{'Recommendation':<15}")
        print("-" * 100)
        
        for i, r in enumerate(sorted_results, 1):
            emoji = "🥇" if i == 1 else "🥈" if i == 2 else "🥉" if i == 3 else "  "
            print(f"{emoji} {i:<4}{r.tech_name:<15}{r.blockbuster_index:<20.1f}"
                  f"{r.maturity_score:<12.1f}{r.market_potential_score:<18.1f}"
                  f"{r.momentum_score:<12.1f}{r.recommendation:<15}")
        
        print("-" * 100)
        
        # Detailed analysis
        print("\n📋 Detailed Assessment Report")
        print("=" * 100)
        
        for i, r in enumerate(sorted_results[:5], 1):  # Top 5 detailed report
            print(f"\n[{i}] {r.tech_name} - Blockbuster Index: {r.blockbuster_index:.1f}")
            print("-" * 80)
            print(f"  📈 Score Details: Maturity({r.maturity_score:.1f}) | Market Potential({r.market_potential_score:.1f}) | Momentum({r.momentum_score:.1f})")
            print(f"  💡 Key Drivers: {', '.join(r.key_drivers[:3])}")
            print(f"  ⚠️  Risk Factors: {', '.join(r.risk_factors[:2])}")
            print(f"  ⏰ Timeline Prediction: {r.timeline_prediction}")
            print(f"  🎯 Investment Recommendation: {r.recommendation}")
        
        print("\n" + "=" * 100)
        print("💡 Investment Recommendation Levels: Strongly Recommended(≥80) | Recommended(60-79) | Watch(40-59) | Cautious(<40)")
        print("=" * 100)
    
    @staticmethod
    def generate_json_report(results: List[PredictionResult], threshold: float = 0) -> str:
        """Generate JSON report"""
        filtered = [r for r in results if r.blockbuster_index >= threshold]
        sorted_results = sorted(filtered, key=lambda x: x.blockbuster_index, reverse=True)
        
        report = {
            "generated_at": datetime.now().isoformat(),
            "total_routes": len(sorted_results),
            "rankings": [
                {
                    "rank": i + 1,
                    "tech_name": r.tech_name,
                    "blockbuster_index": r.blockbuster_index,
                    "maturity_score": r.maturity_score,
                    "market_potential_score": r.market_potential_score,
                    "momentum_score": r.momentum_score,
                    "recommendation": r.recommendation,
                    "key_drivers": r.key_drivers,
                    "risk_factors": r.risk_factors,
                    "timeline_prediction": r.timeline_prediction
                }
                for i, r in enumerate(sorted_results)
            ]
        }
        
        return json.dumps(report, ensure_ascii=False, indent=2)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Blockbuster Therapy Predictor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py                          # Run complete analysis
  python main.py --tech PROTAC,mRNA       # Analyze specified technologies only
  python main.py --output json            # Output in JSON format
  python main.py --threshold 70           # Show only technologies with index ≥70
        """
    )
    
    parser.add_argument(
        "--mode",
        choices=["full", "quick"],
        default="full",
        help="Analysis mode: full=complete analysis, quick=quick analysis"
    )
    parser.add_argument(
        "--tech",
        type=str,
        help="Specify technologies to analyze, comma-separated, e.g.: PROTAC,mRNA,CRISPR"
    )
    parser.add_argument(
        "--output",
        choices=["console", "json"],
        default="console",
        help="Output format"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0,
        help="Minimum blockbuster index threshold (0-100)"
    )
    parser.add_argument(
        "--save",
        type=str,
        help="Save report to file path"
    )
    
    args = parser.parse_args()
    
    # Load data
    print("📥 Loading data...")
    all_techs = DataLoader.load_mock_data()
    
    # Filter specified technologies
    if args.tech:
        target_techs = [t.strip() for t in args.tech.split(",")]
        all_techs = [t for t in all_techs if t.name in target_techs]
        if not all_techs:
            print(f"❌ Specified technologies not found: {args.tech}")
            return
    
    # Execute prediction
    print(f"🔬 Analyzing {len(all_techs)} technology routes...")
    predictor = BlockbusterPredictor()
    results = [predictor.predict(tech) for tech in all_techs]
    
    # Generate report
    if args.output == "json":
        report = ReportGenerator.generate_json_report(results, args.threshold)
        print(report)
        if args.save:
            with open(args.save, "w", encoding="utf-8") as f:
                f.write(report)
            print(f"\n✅ Report saved to: {args.save}")
    else:
        ReportGenerator.generate_console_report(results, args.threshold)
        if args.save:
            json_report = ReportGenerator.generate_json_report(results, args.threshold)
            with open(args.save, "w", encoding="utf-8") as f:
                f.write(json_report)
            print(f"\n✅ JSON report saved to: {args.save}")


if __name__ == "__main__":
    main()
