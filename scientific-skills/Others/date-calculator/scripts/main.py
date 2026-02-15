#!/usr/bin/env python3
"""Date Calculator - Medical date calculations for gestational age and follow-up windows."""

import argparse
import json
import sys
from datetime import datetime, timedelta


class DateCalculator:
    """Calculates medical dates including gestational age and follow-up windows."""
    
    def calculate_gestational_age(self, lmp_date: str) -> dict:
        """Calculate gestational age from last menstrual period.
        
        Args:
            lmp_date: Last menstrual period date (YYYY-MM-DD)
            
        Returns:
            Dictionary with gestational age and estimated delivery date
        """
        try:
            lmp = datetime.strptime(lmp_date, "%Y-%m-%d")
            today = datetime.now()
            
            days = (today - lmp).days
            weeks = days // 7
            remaining_days = days % 7
            
            # Estimated due date (40 weeks from LMP)
            edd = lmp + timedelta(days=280)
            
            return {
                "lmp_date": lmp_date,
                "gestational_age": f"{weeks} weeks {remaining_days} days",
                "gestational_age_days": days,
                "estimated_delivery_date": edd.strftime("%Y-%m-%d"),
                "calculation_date": today.strftime("%Y-%m-%d")
            }
        except ValueError:
            return {"error": "Invalid date format. Use YYYY-MM-DD"}
    
    def calculate_followup_window(self, start_date: str, weeks: int = 4, window_days: int = 7) -> dict:
        """Calculate follow-up date window.
        
        Args:
            start_date: Start date (YYYY-MM-DD)
            weeks: Number of weeks for follow-up (default: 4)
            window_days: Window size in days (default: 7)
            
        Returns:
            Dictionary with follow-up window dates
        """
        try:
            start = datetime.strptime(start_date, "%Y-%m-%d")
            
            window_start = start + timedelta(weeks=weeks)
            window_end = window_start + timedelta(days=window_days)
            
            return {
                "start_date": start_date,
                "followup_weeks": weeks,
                "window_start": window_start.strftime("%Y-%m-%d"),
                "window_end": window_end.strftime("%Y-%m-%d"),
                "window_range": f"{window_start.strftime('%Y-%m-%d')} to {window_end.strftime('%Y-%m-%d')}"
            }
        except ValueError:
            return {"error": "Invalid date format. Use YYYY-MM-DD"}


def main():
    parser = argparse.ArgumentParser(
        description="Date Calculator - Medical date calculations for gestational age and follow-up windows",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate gestational age
  python main.py --type gestational --date 2024-01-15
  
  # Calculate 4-week follow-up window
  python main.py --type followup --date 2024-03-01
  
  # Calculate custom follow-up (6 weeks)
  python main.py --type followup --date 2024-03-01 --weeks 6
        """
    )
    
    parser.add_argument(
        "--type", "-t",
        type=str,
        choices=["gestational", "followup"],
        required=True,
        help="Calculation type (gestational or followup)"
    )
    
    parser.add_argument(
        "--date", "-d",
        type=str,
        required=True,
        help="Date in YYYY-MM-DD format (LMP for gestational, start date for followup)"
    )
    
    parser.add_argument(
        "--weeks",
        type=int,
        default=4,
        help="Number of weeks for follow-up (default: 4)"
    )
    
    parser.add_argument(
        "--window-days",
        type=int,
        default=7,
        help="Follow-up window size in days (default: 7)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output JSON file path (optional)"
    )
    
    args = parser.parse_args()
    
    calc = DateCalculator()
    
    if args.type == "gestational":
        result = calc.calculate_gestational_age(args.date)
    elif args.type == "followup":
        result = calc.calculate_followup_window(args.date, args.weeks, args.window_days)
    else:
        print(f"Error: Unknown calculation type: {args.type}", file=sys.stderr)
        sys.exit(1)
    
    # Check for errors
    if "error" in result:
        print(f"Error: {result['error']}", file=sys.stderr)
        sys.exit(1)
    
    # Output
    output = json.dumps(result, indent=2)
    
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Result saved to: {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
