#!/usr/bin/env python3
"""
CO2 Tank Monitor - CO2 gas cylinder monitoring system
Simulates IoT monitoring to predict tank depletion time and prevent weekend outages.
"""

import argparse
import random
import sys
from datetime import datetime, timedelta


def get_current_time() -> datetime:
    """Get current time (easy to mock for testing)"""
    return datetime.now()


def simulate_sensor_data():
    """Simulate sensor data reading"""
    # Simulate 40L cylinder, full pressure ~15MPa, working pressure 8-10MPa, alarm pressure ~2MPa
    pressure = round(random.uniform(2.5, 12.0), 2)
    capacity = random.choice([10, 40])
    daily_consumption = round(random.uniform(0.5, 3.0), 2)
    return pressure, capacity, daily_consumption


def calculate_remaining_days(pressure: float, daily_consumption: float) -> float:
    """Calculate remaining days"""
    if daily_consumption <= 0:
        return float('inf')
    return pressure / daily_consumption


def calculate_depletion_time(remaining_days: float) -> datetime:
    """Calculate estimated depletion time"""
    return get_current_time() + timedelta(days=remaining_days)


def is_weekend(dt: datetime) -> bool:
    """Check if weekend (Saturday or Sunday)"""
    return dt.weekday() >= 5  # 5=Saturday, 6=Sunday


def will_deplete_on_weekend(depletion_time: datetime, alert_days: int) -> bool:
    """Check if depletion occurs during weekend"""
    now = get_current_time()
    
    # Calculate alert start time
    alert_start = now + timedelta(days=alert_days)
    
    # If depletion time is within alert period and on weekend
    if depletion_time <= alert_start:
        return is_weekend(depletion_time)
    
    # Check if it spans across weekend
    days_until_depletion = (depletion_time - now).days
    for i in range(int(days_until_depletion) + 1):
        check_day = now + timedelta(days=i)
        if is_weekend(check_day) and check_day <= depletion_time:
            # If depletes during weekend
            weekend_start = check_day.replace(hour=0, minute=0, second=0)
            weekend_end = weekend_start + timedelta(days=2)
            if weekend_start <= depletion_time <= weekend_end:
                return True
    
    return False


def get_status(remaining_days: float, alert_days: int, depletion_time: datetime) -> tuple:
    """
    Get status code and description
    Returns: (status_code, status_icon, status_text, recommendations)
    """
    if remaining_days <= 0:
        return 2, "🔴", "DEPLETED", ["⚠️  Tank is empty! Replace immediately!"]
    
    weekend_risk = will_deplete_on_weekend(depletion_time, alert_days)
    
    if remaining_days <= alert_days or weekend_risk:
        recommendations = []
        if remaining_days <= alert_days:
            recommendations.append(f"⚠️  Tank will deplete in {remaining_days:.1f} days")
        if weekend_risk:
            recommendations.append("⚠️  Tank will deplete on weekend!")
        recommendations.append("💡 Recommendation: Replace tank immediately or arrange weekend monitoring")
        return 2, "🔴", "DANGER", recommendations
    
    if remaining_days <= alert_days + 2:
        return 1, "🟡", "WARNING", [
            f"⏰ Remaining time: {remaining_days:.1f} days",
            "💡 Recommendation: Monitor tank status, prepare for replacement"
        ]
    
    return 0, "🟢", "OK", [
        f"✅ Sufficient remaining time: {remaining_days:.1f} days",
        "💡 No immediate action required"
    ]


def format_report(
    pressure: float,
    capacity: int,
    daily_consumption: float,
    remaining_days: float,
    depletion_time: datetime,
    status_code: int,
    status_icon: str,
    status_text: str,
    recommendations: list
) -> str:
    """Format report"""
    now = get_current_time()
    
    # Format depletion time with weekday
    weekdays = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    weekday_str = weekdays[depletion_time.weekday()]
    depletion_str = depletion_time.strftime(f"%Y-%m-%d %H:%M ({weekday_str})")
    
    report_lines = [
        "=" * 40,
        "       CO2 Tank Monitor Report",
        "=" * 40,
        f"📅 Current Time: {now.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "📊 Sensor Data:",
        f"   Current Pressure: {pressure:.2f} MPa",
        f"   Tank Capacity: {capacity} L",
        f"   Daily Consumption: {daily_consumption:.2f} MPa/day",
        "",
        "⏱️  Prediction Analysis:",
        f"   Estimated Remaining Days: {remaining_days:.1f} days",
        f"   Estimated Depletion Time: {depletion_str}",
        "",
        f"🚨 Alert Status: {status_icon} {status_text}",
    ]
    
    for rec in recommendations:
        report_lines.append(f"   {rec}")
    
    report_lines.append("=" * 40)
    
    return "\n".join(report_lines)


def main():
    parser = argparse.ArgumentParser(
        description="CO2 Tank Monitor - Monitor CO2 gas cylinder status to prevent weekend outages"
    )
    parser.add_argument(
        "--pressure", "-p",
        type=float,
        default=8.0,
        help="Current tank pressure (MPa), default 8.0"
    )
    parser.add_argument(
        "--capacity", "-c",
        type=int,
        default=40,
        choices=[10, 40],
        help="Tank capacity (L), default 40"
    )
    parser.add_argument(
        "--daily-consumption", "-d",
        type=float,
        default=1.5,
        help="Daily consumption rate (MPa/day), default 1.5"
    )
    parser.add_argument(
        "--alert-days", "-a",
        type=int,
        default=2,
        help="Alert threshold in days, default 2"
    )
    parser.add_argument(
        "--simulate", "-s",
        action="store_true",
        help="Enable simulation mode (random data generation)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Quiet mode, output only return code"
    )
    
    args = parser.parse_args()
    
    # Get data
    if args.simulate:
        pressure, capacity, daily_consumption = simulate_sensor_data()
    else:
        pressure = args.pressure
        capacity = args.capacity
        daily_consumption = args.daily_consumption
    
    # Calculate
    remaining_days = calculate_remaining_days(pressure, daily_consumption)
    depletion_time = calculate_depletion_time(remaining_days)
    
    # Get status
    status_code, status_icon, status_text, recommendations = get_status(
        remaining_days, args.alert_days, depletion_time
    )
    
    # Output report
    if not args.quiet:
        report = format_report(
            pressure, capacity, daily_consumption,
            remaining_days, depletion_time,
            status_code, status_icon, status_text, recommendations
        )
        print(report)
    
    # Return status code
    sys.exit(status_code)


if __name__ == "__main__":
    main()
