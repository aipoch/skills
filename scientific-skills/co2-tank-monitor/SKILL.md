---
name: co2-tank-monitor
description: IoT monitoring simulation to predict CO2 tank depletion and prevent weekend gas outages in cell culture facilities. Monitors cylinder pressure, calculates consumption rates, and provides early warnings for timely replacement.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
  skill-author: AIPOCH
---

# CO2 Tank Monitor

Monitor CO2 cylinder pressure and predict depletion times to prevent gas outages in cell culture incubators, particularly during weekends when laboratories are unmanned. Provides automated alerts and consumption tracking for proactive cylinder management.

**Key Capabilities:**
- **Pressure-Based Depletion Prediction**: Calculate remaining cylinder life based on current pressure and consumption rates
- **Weekend Risk Detection**: Identify if depletion will occur during weekends when no staff is available
- **Multi-Cylinder Support**: Handle different cylinder sizes (10L, 40L) and specifications
- **Automated Alert System**: Generate status reports with actionable recommendations
- **Simulation Mode**: Test monitoring scenarios with randomized data for training

---

## When to Use

**✅ Use this skill when:**
- Monitoring **cell culture CO2 incubators** that require continuous gas supply
- Setting up **automated daily checks** for cylinder status (e.g., via cron jobs)
- Planning **cylinder replacement schedules** to avoid service interruptions
- **Training new lab members** on gas monitoring procedures (use simulation mode)
- Managing **multiple incubators** with different consumption rates
- **Pre-holiday assessments** to ensure adequate gas supply during lab closures
- Creating **standard operating procedures** (SOPs) for gas monitoring

**❌ Do NOT use when:**
- Monitoring **liquid nitrogen tanks** for cryogenic storage → Use `freezer-sample-locator` or specialized LN2 monitors
- Tracking **compressed air or N2 cylinders** for non-CO2 applications → Adapt parameters but verify compatibility
- Real-time **IoT sensor integration** is available → Use actual sensor API instead of simulation
- Cylinder uses **weight-based measurement** instead of pressure → Modify calculation logic accordingly
- **Emergency gas leak detection** is needed → Use gas detection safety systems with alarms
- Managing **bulk tank installations** with automatic switchover → These systems have built-in monitoring

**Related Skills:**
- **上游 (Upstream)**: `equipment-maintenance-log`, `lab-inventory-tracker`
- **下游 (Downstream)**: `eln-template-creator`, `lab-budget-forecaster`

---

## Integration with Other Skills

**Upstream Skills:**
- `equipment-maintenance-log`: Track CO2 incubator maintenance history that affects consumption rates
- `lab-inventory-tracker`: Monitor cylinder inventory and replacement schedules
- `equipment-maintenance-log`: Record calibration dates for pressure gauges

**Downstream Skills:**
- `eln-template-creator`: Generate experiment templates that include gas monitoring checklists
- `lab-budget-forecaster`: Forecast CO2 gas costs based on consumption trends
- `waste-disposal-guide`: Coordinate empty cylinder disposal with replacement

**Complete Workflow:**
```
Daily Morning Check → co2-tank-monitor → equipment-maintenance-log → lab-inventory-tracker → Replacement Scheduling
```

---

## Core Capabilities

### 1. Pressure-Based Depletion Prediction

Calculate remaining cylinder lifetime based on current pressure readings and historical consumption patterns.

```python
from scripts.main import calculate_remaining_days, calculate_depletion_time
from datetime import datetime

# Current cylinder status
current_pressure = 8.0  # MPa
daily_consumption = 1.5  # MPa/day

# Calculate remaining time
remaining_days = calculate_remaining_days(current_pressure, daily_consumption)
print(f"Remaining days: {remaining_days:.1f}")

# Calculate depletion datetime
depletion_time = calculate_depletion_time(remaining_days)
print(f"Estimated depletion: {depletion_time.strftime('%Y-%m-%d %H:%M')}")

# Formula: remaining_days = pressure / daily_consumption
```

**Parameters:**

| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `pressure` | float | Yes | Current cylinder pressure in MPa | 8.0 |
| `daily_consumption` | float | Yes | Average daily consumption rate (MPa/day) | 1.5 |
| `capacity` | int | No | Cylinder capacity in liters (10 or 40) | 40 |

**Calculation Method:**
```
remaining_days = current_pressure / daily_consumption
```

**Best Practices:**
- ✅ **Calibrate consumption rate** based on historical data (monitor for 1-2 weeks)
- ✅ **Account for usage variations** - higher on heavy culture days, lower on weekends
- ✅ **Check pressure at consistent time** (e.g., 9 AM daily) for accurate trends
- ✅ **Update consumption rate seasonally** - incubator usage may vary

**Common Issues and Solutions:**

**Issue: Inaccurate consumption estimates**
- Symptom: Predicted depletion date consistently off by several days
- Solution: Monitor actual consumption over 2-4 weeks; adjust daily_consumption parameter

**Issue: Pressure fluctuations causing false alerts**
- Symptom: Inconsistent pressure readings due to temperature changes
- Solution: Take readings at same time daily; allow cylinder to equilibrate after use

### 2. Weekend Depletion Risk Detection

Identify if cylinder depletion will occur during weekends when laboratory staff is unavailable.

```python
from scripts.main import is_weekend, will_deplete_on_weekend, calculate_depletion_time
from datetime import datetime

# Calculate depletion time
remaining_days = 3.5
depletion_time = calculate_depletion_time(remaining_days)

# Check if depletion is on weekend
if is_weekend(depletion_time):
    print(f"⚠️  Warning: Depletion on {depletion_time.strftime('%A')}")
else:
    print(f"✅ Depletion on weekday: {depletion_time.strftime('%A')}")

# Check weekend risk with alert threshold
alert_days = 2
weekend_risk = will_deplete_on_weekend(depletion_time, alert_days)

if weekend_risk:
    print("🔴 CRITICAL: Cylinder will deplete during weekend!")
    print("   Action: Replace before Friday or arrange weekend duty")
```

**Weekend Risk Scenarios:**

| Scenario | Risk Level | Action Required |
|----------|------------|-----------------|
| Depletion on Saturday/Sunday | 🔴 High | Immediate replacement or weekend duty |
| Depletion Monday morning | 🟡 Medium | Replace Friday afternoon |
| Depletion Friday evening | 🟡 Medium | Monitor closely; replace if possible |
| Depletion mid-week | 🟢 Low | Schedule routine replacement |

**Best Practices:**
- ✅ **Set alert_days = 2** for standard work weeks (alert by Wednesday for Friday replacement)
- ✅ **Plan Friday replacements** for any cylinder with <4 days remaining
- ✅ **Create weekend duty roster** for critical experiments during high-risk periods
- ✅ **Consider holiday schedules** - extend alert_days before long weekends

**Common Issues and Solutions:**

**Issue: False weekend alerts for Monday depletion**
- Symptom: Alerts triggered for Monday morning depletion when lab reopens
- Solution: Adjust logic or manually review; Monday AM depletion may be acceptable

**Issue: Missing Friday evening depletion risk**
- Symptom: Cylinder runs out Friday night after staff leaves
- Solution: Check Friday 5 PM status specifically; consider Friday evening as weekend

### 3. Multi-Cylinder Size Support

Support different CO2 cylinder specifications commonly used in laboratory settings.

```python
from scripts.main import simulate_sensor_data

# Cylinder specifications
CYLINDER_SPECS = {
    10: {
        "name": "Portable (10L)",
        "full_pressure": 15.0,  # MPa
        "typical_usage": "Small incubators, backup",
        "duration_1_5_mpa": "~10 days"
    },
    40: {
        "name": "Standard (40L)",
        "full_pressure": 15.0,  # MPa  
        "typical_usage": "Main incubators, high usage",
        "duration_1_5_mpa": "~40 days"
    }
}

# Select appropriate cylinder
capacity = 40  # Liters
specs = CYLINDER_SPECS[capacity]
print(f"Cylinder: {specs['name']}")
print(f"Full pressure: {specs['full_pressure']} MPa")
print(f"Typical usage: {specs['typical_usage']}")

# Calculate capacity-based estimates
full_pressure = 15.0  # MPa (typical full cylinder)
daily_consumption = 1.5  # MPa/day
max_days = full_pressure / daily_consumption
print(f"Maximum duration at {daily_consumption} MPa/day: {max_days:.1f} days")
```

**Cylinder Specifications:**

| Capacity | Typical Use | Full Pressure | Duration (@1.5 MPa/day) |
|----------|-------------|---------------|------------------------|
| **10L** | Small incubators, backup | ~15 MPa | ~10 days |
| **40L** | Main incubators, heavy use | ~15 MPa | ~40 days |

**Best Practices:**
- ✅ **Use 40L cylinders** for main incubators to reduce replacement frequency
- ✅ **Keep 10L cylinders** as emergency backup for critical experiments
- ✅ **Record cylinder capacity** in monitoring logs for accurate predictions
- ✅ **Standardize on one size** per facility when possible for simplified inventory

**Common Issues and Solutions:**

**Issue: Wrong capacity setting causing prediction errors**
- Symptom: Predictions consistently off by 4x (10L vs 40L confusion)
- Solution: Verify cylinder label; update capacity parameter in monitoring script

**Issue: Mixing cylinder sizes in same calculation**
- Symptom: Different incubators on different cylinder sizes
- Solution: Run separate monitoring instances for each cylinder/incubator pair

### 4. Automated Alert System with Status Levels

Generate color-coded status reports with specific recommendations based on urgency.

```python
from scripts.main import get_status, calculate_remaining_days, calculate_depletion_time

# Example scenarios
scenarios = [
    {"pressure": 12.0, "consumption": 1.5, "days": 2},  # Normal
    {"pressure": 5.0, "consumption": 1.5, "days": 2},   # Caution  
    {"pressure": 2.5, "consumption": 1.5, "days": 2},  # Danger
]

for scenario in scenarios:
    remaining = calculate_remaining_days(
        scenario["pressure"], 
        scenario["consumption"]
    )
    depletion = calculate_depletion_time(remaining)
    
    status_code, icon, status_text, recommendations = get_status(
        remaining, scenario["days"], depletion
    )
    
    print(f"\nPressure: {scenario['pressure']} MPa")
    print(f"Status: {icon} {status_text} (Code: {status_code})")
    print("Recommendations:")
    for rec in recommendations:
        print(f"  {rec}")
```

**Status Levels:**

| Code | Icon | Status | Condition | Action |
|------|------|--------|-----------|--------|
| **0** | 🟢 | Normal | Days > alert_days + 2 | No action needed |
| **1** | 🟡 | Caution | Days within alert_days + 2 | Monitor closely |
| **2** | 🔴 | Danger | Days ≤ alert_days or weekend depletion | Replace immediately |

**Return Codes for Automation:**

| Exit Code | Meaning | Automation Action |
|-----------|---------|-------------------|
| **0** | Normal | Continue routine monitoring |
| **1** | Caution | Send email notification |
| **2** | Danger | Send urgent alert; page on-call staff |

**Best Practices:**
- ✅ **Integrate with alerting systems** using return codes for automated notifications
- ✅ **Customize recommendations** for your lab's specific procedures
- ✅ **Escalate based on status** - email for caution, SMS/call for danger
- ✅ **Log all status changes** for trend analysis and audit trails

**Common Issues and Solutions:**

**Issue: Alert fatigue from frequent caution alerts**
- Symptom: Too many caution-level notifications causing ignored alerts
- Solution: Adjust alert_days threshold; batch daily reports instead of immediate alerts

**Issue: Missing critical weekend alerts**
- Symptom: Weekend depletions not detected when they span multiple days
- Solution: Ensure will_deplete_on_weekend() logic covers multi-day weekend periods

### 5. Simulation Mode for Training

Generate randomized cylinder data for training staff without affecting real monitoring systems.

```python
from scripts.main import simulate_sensor_data, calculate_remaining_days

# Generate 5 random scenarios for training
print("Training Scenarios (Simulated Data):\n")

for i in range(5):
    pressure, capacity, consumption = simulate_sensor_data()
    remaining = calculate_remaining_days(pressure, consumption)
    
    print(f"Scenario {i+1}:")
    print(f"  Pressure: {pressure:.2f} MPa")
    print(f"  Capacity: {capacity} L")
    print(f"  Daily consumption: {consumption:.2f} MPa/day")
    print(f"  Estimated remaining: {remaining:.1f} days")
    
    # Training question
    if remaining < 3:
        print(f"  🚨 ACTION: Replace immediately!")
    elif remaining < 7:
        print(f"  ⚠️  ACTION: Schedule replacement this week")
    else:
        print(f"  ✅ STATUS: Monitor normally")
    print()

# Command line usage for simulation mode:
# python scripts/main.py --simulate
```

**Training Applications:**

| Use Case | Simulation Parameters | Learning Objective |
|----------|----------------------|-------------------|
| **New staff training** | Various random scenarios | Recognize different alert levels |
| **Weekend risk awareness** | Force weekend depletion scenarios | Understand critical timing |
| **Consumption calculation** | Different pressure/consumption combos | Practice manual calculations |
| **Emergency response** | Low pressure scenarios (<3 MPa) | Learn urgent procedures |

**Best Practices:**
- ✅ **Run simulations weekly** to keep staff familiar with alert formats
- ✅ **Create scenario library** with specific training cases (weekend depletion, low pressure, etc.)
- ✅ **Include simulation in onboarding** for new lab members
- ✅ **Test alert pathways** using simulation to verify notification systems

**Common Issues and Solutions:**

**Issue: Simulation data too random for consistent training**
- Symptom: Trainees see different scenarios each session
- Solution: Set random seed for reproducible scenarios; create predefined training datasets

**Issue: Confusion between simulation and real data**
- Symptom: Staff mistakes simulation for actual cylinder status
- Solution: Clearly label all simulation outputs; use distinct formatting

### 6. Automated Scheduling Integration

Integrate with cron jobs or task schedulers for automated daily monitoring.

```python
# Example cron setup (documented in comments)
"""
Cron job configuration for daily monitoring:

# Daily check at 9:00 AM
0 9 * * * cd /lab/scripts && python co2_tank_monitor.py --pressure $(cat /sensor/pressure.log | tail -1) --quiet

# Check before weekends (Friday at 5 PM)
0 17 * * 5 cd /lab/scripts && python co2_tank_monitor.py --pressure $(cat /sensor/pressure.log | tail -1)

# Log results
0 9 * * * cd /lab/scripts && python co2_tank_monitor.py >> /var/log/co2_monitor.log 2>&1
"""

# Integration with sensor reading
def read_sensor_data(sensor_log_path: str) -> float:
    """Read current pressure from sensor log file."""
    try:
        with open(sensor_log_path, 'r') as f:
            lines = f.readlines()
            # Get last line and extract pressure value
            last_line = lines[-1].strip()
            pressure = float(last_line.split()[0])  # Assumes format: "8.5 MPa"
            return pressure
    except Exception as e:
        print(f"Error reading sensor: {e}")
        return None

# Usage in automated script
sensor_pressure = read_sensor_data("/path/to/pressure_sensor.log")
if sensor_pressure:
    # Run monitoring with actual sensor data
    import subprocess
    result = subprocess.run([
        "python", "scripts/main.py",
        "--pressure", str(sensor_pressure),
        "--quiet"
    ], capture_output=True)
    
    # Handle return code
    if result.returncode == 2:
        send_urgent_alert("CO2 cylinder requires immediate replacement!")
    elif result.returncode == 1:
        send_notification("CO2 cylinder needs attention soon")
```

**Integration Patterns:**

| Integration Type | Trigger | Action | Return Code Handling |
|-----------------|---------|--------|---------------------|
| **Daily cron** | 9:00 AM daily | Check status | Log results; alert if ≠0 |
| **Pre-weekend** | Friday 5 PM | Weekend risk check | Force alert if weekend risk |
| **Real-time** | Sensor threshold | Immediate check | Page on-call if danger |
| **Manual** | Lab staff | Ad-hoc check | Display full report |

**Best Practices:**
- ✅ **Use --quiet mode** for automated runs to reduce log noise
- ✅ **Capture return codes** for alerting and logging
- ✅ **Test cron jobs** manually before deploying to production
- ✅ **Monitor the monitor** - ensure scheduled checks are running
- ✅ **Implement heartbeat** - alert if monitoring script fails to run

**Common Issues and Solutions:**

**Issue: Cron job not executing**
- Symptom: No monitoring logs or alerts
- Causes: Path issues, permissions, environment variables
- Solution: Use full paths; test manually first; check cron logs

**Issue: Sensor data unavailable**
- Symptom: Cannot read pressure from sensor log
- Solution: Implement fallback to manual input; alert on sensor failure

---

## Complete Workflow Example

**From daily monitoring to replacement scheduling:**

```bash
# Step 1: Manual morning check with full report
python scripts/main.py --pressure 8.5 --daily-consumption 1.2

# Step 2: Automated daily check (in cron)
0 9 * * * python scripts/main.py --pressure $(cat sensor.log | tail -1) --quiet

# Step 3: Pre-weekend check
python scripts/main.py --pressure 5.5 --alert-days 3

# Step 4: Simulation for training
python scripts/main.py --simulate
```

**Python API Usage:**

```python
from scripts.main import (
    calculate_remaining_days,
    calculate_depletion_time,
    will_deplete_on_weekend,
    get_status
)
from datetime import datetime

def monitor_co2_cylinder(
    pressure: float,
    daily_consumption: float,
    capacity: int = 40,
    alert_days: int = 2
) -> dict:
    """
    Complete CO2 cylinder monitoring workflow.
    
    Returns:
        Dictionary with status, predictions, and recommendations
    """
    # Calculate predictions
    remaining_days = calculate_remaining_days(pressure, daily_consumption)
    depletion_time = calculate_depletion_time(remaining_days)
    
    # Assess weekend risk
    weekend_risk = will_deplete_on_weekend(depletion_time, alert_days)
    
    # Get status and recommendations
    status_code, icon, status_text, recommendations = get_status(
        remaining_days, alert_days, depletion_time
    )
    
    # Compile report
    report = {
        "timestamp": datetime.now().isoformat(),
        "cylinder": {
            "pressure_mpa": pressure,
            "capacity_l": capacity,
            "daily_consumption_mpa": daily_consumption
        },
        "prediction": {
            "remaining_days": round(remaining_days, 1),
            "depletion_time": depletion_time.isoformat(),
            "weekend_risk": weekend_risk
        },
        "status": {
            "code": status_code,
            "level": status_text,
            "icon": icon
        },
        "recommendations": recommendations,
        "action_required": status_code >= 1
    }
    
    return report

# Execute monitoring
result = monitor_co2_cylinder(
    pressure=6.5,
    daily_consumption=1.5,
    capacity=40,
    alert_days=2
)

# Display formatted report
print(f"CO2 Cylinder Monitor Report")
print(f"{'='*40}")
print(f"Status: {result['status']['icon']} {result['status']['level']}")
print(f"Remaining: {result['prediction']['remaining_days']} days")
print(f"Depletion: {result['prediction']['depletion_time']}")
print(f"Weekend Risk: {'Yes' if result['prediction']['weekend_risk'] else 'No'}")
print(f"\nRecommendations:")
for rec in result['recommendations']:
    print(f"  {rec}")
```

**Expected Output Files:**

```
monitoring_logs/
├── co2_daily_checks.log      # Daily monitoring results
├── co2_weekend_alerts.log    # Weekend-specific alerts
└── co2_replacement_history.json  # Cylinder change tracking
```

---

## Common Patterns

### Pattern 1: Daily Morning Monitoring Routine

**Scenario**: Lab technician checks all CO2 cylinders every morning at 9 AM.

```json
{
  "monitoring_type": "daily_routine",
  "schedule": "9:00 AM daily",
  "cylinders": [
    {"location": "Incubator A", "capacity": 40, "typical_consumption": 1.5},
    {"location": "Incubator B", "capacity": 40, "typical_consumption": 1.2},
    {"location": "Backup", "capacity": 10, "typical_consumption": 0.3}
  ],
  "alert_threshold": 2,
  "actions": {
    "normal": "Log and continue",
    "caution": "Schedule replacement within 3 days",
    "danger": "Replace immediately or arrange weekend coverage"
  }
}
```

**Workflow:**
1. Read pressure gauges on all cylinders at 9 AM
2. Run monitoring script for each cylinder
3. Log status codes and any alerts
4. If any cylinder shows status ≥1, notify lab manager
5. Schedule replacements for caution-level cylinders
6. For danger-level, initiate immediate replacement procedure

**Output Example:**
```
Morning CO2 Check - 2026-02-09 09:00
=====================================
Incubator A (40L): 🟢 Normal - 12.5 days remaining
Incubator B (40L): 🟡 Caution - 4.2 days remaining  
Backup (10L): 🟢 Normal - 8.1 days remaining

Action Required:
  - Schedule Incubator B cylinder replacement for this week
```

### Pattern 2: Pre-Holiday Weekend Assessment

**Scenario**: Before a 3-day weekend (e.g., Memorial Day), ensure adequate gas supply.

```json
{
  "monitoring_type": "pre_holiday",
  "holiday": "Memorial Day Weekend",
  "lab_closure": "3 days (Sat-Mon)",
  "alert_days": 4,
  "special_considerations": [
    "No staff on-site for 72 hours",
    "Critical experiments running",
    "Extended alert threshold"
  ]
}
```

**Workflow:**
1. On Friday before holiday, run assessment with extended alert_days (4)
2. Check all cylinders with weekend risk detection
3. For any cylinder with depletion during closure period:
   - Prioritize immediate replacement
   - Consider emergency weekend delivery
   - Transfer critical cultures to backup incubator
4. Document gas status in holiday handoff notes
5. Set up emergency contact for gas supplier

**Output Example:**
```
Pre-Holiday Assessment - Memorial Day Weekend
==============================================
🔴 CRITICAL: Incubator A cylinder will deplete on Sunday
   Current pressure: 4.5 MPa
   Depletion: Sunday 11:30 PM
   
⚠️  Lab closure: 3 days with no staff

IMMEDIATE ACTIONS:
  1. Replace Incubator A cylinder TODAY (Friday)
  2. Verify backup incubator is operational
  3. Contact gas supplier for emergency weekend number
  4. Transfer critical samples to Incubator B if concerned
```

### Pattern 3: New Lab Member Training

**Scenario**: Train new technician on CO2 monitoring using simulation mode.

```json
{
  "training_type": "new_staff_onboarding",
  "mode": "simulation",
  "scenarios": [
    {"name": "Normal operation", "pressure_range": [10, 15]},
    {"name": "Approaching replacement", "pressure_range": [4, 6]},
    {"name": "Critical weekend risk", "pressure_range": [2, 4]},
    {"name": "Emergency low pressure", "pressure_range": [0.5, 2]}
  ],
  "learning_objectives": [
    "Recognize different alert levels",
    "Calculate remaining days manually",
    "Understand weekend risk",
    "Learn replacement procedures"
  ]
}
```

**Workflow:**
1. Run simulation mode to generate 10 random scenarios
2. For each scenario, trainee manually calculates remaining days
3. Compare trainee calculation with script output
4. Discuss appropriate actions for each status level
5. Practice emergency procedures with "danger" scenarios
6. Review actual monitoring logs from past month
7. Shadow experienced technician for first week

**Output Example:**
```
Training Session - CO2 Monitoring
==================================

Scenario 1/5: Normal Operation
  Simulated pressure: 12.3 MPa
  Daily consumption: 1.4 MPa/day
  
  Trainee calculation: 8.8 days remaining
  Script result: 8.8 days remaining ✅
  
  Discussion: What actions needed?
  ✓ Correct: Continue routine monitoring

Scenario 2/5: Weekend Risk
  Simulated pressure: 3.8 MPa
  Daily consumption: 1.5 MPa/day
  Depletion: Sunday 2 PM
  
  Status: 🔴 Danger
  
  Discussion: What if this were Friday morning?
  ✓ Correct: Immediate replacement required
```

### Pattern 4: Multi-Incubator Facility Management

**Scenario**: Large facility with 6 incubators on different CO2 cylinders.

```json
{
  "facility_type": "multi_incubator",
  "total_incubators": 6,
  "cylinder_configuration": {
    "main_incubators": {"count": 4, "capacity": 40, "consumption": "1.5-2.0 MPa/day"},
    "backup_incubators": {"count": 2, "capacity": 10, "consumption": "0.3-0.5 MPa/day"}
  },
  "monitoring_strategy": "centralized_dashboard",
  "rotation_schedule": "staggered_replacement"
}
```

**Workflow:**
1. Create monitoring dashboard tracking all 6 cylinders
2. Run monitoring script for each cylinder daily
3. Implement staggered replacement schedule to avoid all cylinders needing replacement simultaneously
4. Track consumption trends to optimize delivery schedule
5. Maintain 1-2 backup cylinders in inventory
6. Negotiate bulk pricing with gas supplier based on usage patterns
7. Monthly review of consumption trends and cost optimization

**Output Example:**
```
Facility CO2 Dashboard - 2026-02-09
====================================

Main Incubators (40L):
  Inc-01: 🟢 18.2 days | Last replaced: 2026-01-22
  Inc-02: 🟡 4.5 days | REPLACE THIS WEEK
  Inc-03: 🟢 22.1 days | Last replaced: 2026-01-15
  Inc-04: 🟢 15.8 days | Last replaced: 2026-01-28

Backup Incubators (10L):
  Back-01: 🟢 6.2 days
  Back-02: 🟢 8.4 days

This Week Actions:
  - Replace Inc-02 cylinder by Wednesday
  - Order replacement for delivery next week
  
Monthly Consumption: 42 MPa (28% under budget)
```

---

## Quality Checklist

**Pre-Monitoring Setup:**
- [ ] **CRITICAL**: Verify pressure gauge is calibrated and functional
- [ ] Confirm cylinder size (10L or 40L) matches monitoring parameters
- [ ] Establish baseline daily consumption rate through 1-2 weeks observation
- [ ] Set appropriate alert_days threshold (typically 2 for standard weeks)
- [ ] Configure automated scheduling (cron) for daily checks
- [ ] Set up alerting pathway (email, SMS) for danger-level status
- [ ] Create emergency contact list for gas supplier
- [ ] Ensure backup cylinder availability for critical incubators

**During Daily Monitoring:**
- [ ] Take pressure readings at consistent time (e.g., 9:00 AM)
- [ ] **CRITICAL**: Record actual pressure value, not estimates
- [ ] Note any unusual consumption patterns (door openings, new cultures)
- [ ] Verify status code and understand implications
- [ ] Check for weekend depletion risk (especially on Thursdays/Fridays)
- [ ] Log all readings for trend analysis
- [ ] Update consumption rate if usage patterns change significantly
- [ ] Inspect cylinder and regulator for leaks or damage

**Alert Response:**
- [ ] **CRITICAL**: Acknowledge danger-level alerts immediately
- [ ] For weekend depletion risk, arrange replacement before Friday 5 PM
- [ ] Notify lab manager of caution-level cylinders needing replacement
- [ ] Document all actions taken in response to alerts
- [ ] Verify replacement cylinder is available before removing empty one
- [ ] After replacement, verify incubator maintains proper CO2 concentration
- [ ] Update monitoring logs with replacement date and new cylinder info
- [ ] Return empty cylinder to storage area for pickup

**Post-Replacement Verification:**
- [ ] **CRITICAL**: Verify new cylinder valve is fully open
- [ ] Check pressure gauge reading on new cylinder (should be ~15 MPa)
- [ ] Monitor incubator CO2 concentration for 30 minutes to verify stability
- [ ] Verify no gas leaks at regulator connections (soap bubble test)
- [ ] Update monitoring parameters with new cylinder start pressure
- [ ] Log replacement in inventory management system
- [ ] Notify relevant lab members of cylinder change
- [ ] Schedule next replacement based on predicted depletion date

---

## Common Pitfalls

**Monitoring Setup Issues:**
- ❌ **Inconsistent reading times** → Daily variations affect trend accuracy
  - ✅ Take readings at same time each day (e.g., 9:00 AM ± 30 min)
  
- ❌ **Wrong consumption estimates** → Inaccurate depletion predictions
  - ✅ Calculate from actual usage over 2+ weeks; update seasonally
  
- ❌ **Ignoring temperature effects** → Pressure varies with ambient temperature
  - ✅ Allow cylinder to equilibrate to room temperature before reading
  
- ❌ **Not accounting for usage variations** → Weekend vs weekday consumption differs
  - ✅ Monitor separately for different usage patterns; use average

**Alert and Response Issues:**
- ❌ **Alert fatigue** → Too many notifications cause important alerts to be missed
  - ✅ Batch daily reports; only escalate urgent alerts immediately
  
- ❌ **Missing weekend alerts** → Friday check doesn't catch Saturday depletion
  - ✅ Always check weekend risk explicitly; use Friday 5 PM check
  
- ❌ **Delayed replacement** → Waiting too long leads to actual depletion
  - ✅ Replace when caution status triggered; don't wait for danger
  
- ❌ **No backup plan** → Cylinder fails with no replacement available
  - ✅ Maintain 1-2 backup cylinders; know emergency supplier contacts

**Calculation and Data Issues:**
- ❌ **Wrong cylinder capacity** → 10L vs 40L confusion causes 4x error
  - ✅ Always verify cylinder label; double-check capacity parameter
  
- ❌ **Pressure unit confusion** → PSI vs MPa or Bar mixing
  - ✅ Standardize on MPa; convert if gauge shows different units
  
- ❌ **Not tracking multiple cylinders** → Mixing readings from different tanks
  - ✅ Label cylinders clearly; use location-specific monitoring
  
- ❌ **Ignoring cumulative error** → Small daily errors compound over weeks
  - ✅ Verify predictions against actual depletion; calibrate regularly

**Operational Issues:**
- ❌ **Gauges not calibrated** → Inaccurate readings lead to wrong predictions
  - ✅ Calibrate pressure gauges annually or per manufacturer recommendation
  
- ❌ **Leaks undetected** → Higher consumption than predicted
  - ✅ Regular leak checks; sudden consumption increases indicate leaks
  
- ❌ **Regulator problems** → Pressure drops even with adequate gas
  - ✅ Replace regulators per schedule; check for ice buildup in flow
  
- ❌ **No documentation** → Cannot track trends or troubleshoot issues
  - ✅ Maintain detailed logs; review monthly for optimization

---

## Troubleshooting

**Problem: Predicted depletion date consistently wrong**
- Symptoms: Actual depletion 3-5 days earlier or later than predicted
- Causes:
  - Incorrect consumption rate estimate
  - Changing usage patterns (more/fewer cultures)
  - Temperature variations affecting pressure readings
  - Leaks in system
- Solutions:
  - Recalculate consumption rate from recent actual usage
  - Monitor for 2+ weeks to establish accurate baseline
  - Take readings at consistent temperature conditions
  - Check for leaks using soapy water at connections

**Problem: Frequent false alerts (alert fatigue)**
- Symptoms: Receiving caution/danger alerts that don't require action
- Causes:
  - Alert threshold too sensitive
  - Normal fluctuations triggering alerts
  - Multiple notifications for same cylinder
- Solutions:
  - Adjust alert_days parameter (increase from 2 to 3)
  - Implement alert aggregation (daily summary vs immediate)
  - Use quiet mode for automated checks; manual for full reports

**Problem: Weekend depletion not detected**
- Symptoms: Cylinder runs out Saturday/Sunday despite monitoring
- Causes:
  - Friday check missed weekend risk
  - Depletion spans multiple days
  - Holiday weekends not considered
- Solutions:
  - Add explicit Friday 5 PM check for weekend risk
  - Extend alert_days before long weekends
  - Use will_deplete_on_weekend() logic that checks multi-day periods

**Problem: Cannot read pressure gauge accurately**
- Symptoms: Inconsistent or unclear pressure readings
- Causes:
  - Gauge needle stuck or damaged
  - Condensation or dirt on gauge face
  - Parallax error reading needle position
- Solutions:
  - Replace faulty gauges immediately
  - Clean gauge face regularly
  - View needle straight-on to avoid parallax
  - Consider digital pressure sensors for better accuracy

**Problem: Cylinder depletes faster than expected**
- Symptoms: Predicted 10 days, actually depleted in 6 days
- Causes:
  - Increased incubator usage (more cultures, frequent door openings)
  - Gas leak in system
  - Regulator malfunction
  - Temperature increase (higher consumption)
- Solutions:
  - Check for leaks at all connections
  - Verify incubator door seals
  - Inspect regulator for proper function
  - Update consumption rate to reflect actual usage

**Problem: Monitoring script not running automatically**
- Symptoms: No daily logs or alerts despite cron setup
- Causes:
  - Cron job misconfigured
  - Path issues in script
  - Permissions problems
  - Script errors not logged
- Solutions:
  - Test cron job manually first
  - Use full absolute paths in cron and script
  - Check cron logs: `grep CRON /var/log/syslog`
  - Redirect stderr to log file for debugging

---

## References

Available in `references/` directory:

- (No reference files currently available for this skill)

**External Resources:**
- Cell Culture CO2 Guidelines: https://www.thermofisher.com/cellculture
- Gas Cylinder Safety: https://www.osha.gov/gascylinders
- CO2 Incubator Best Practices: https://www.sigmaaldrich.com/incubators

---

## Scripts

Located in `scripts/` directory:

- `main.py` - CO2 tank monitoring engine with prediction and alerting logic

---

## Pressure Conversion Reference

| Unit | MPa | PSI | Bar |
|------|-----|-----|-----|
| **MPa** | 1.0 | 145.0 | 10.0 |
| **PSI** | 0.0069 | 1.0 | 0.069 |
| **Bar** | 0.1 | 14.5 | 1.0 |

**Typical Cylinder Pressures:**
- **Full cylinder**: ~15 MPa (~2200 PSI)
- **Working pressure**: 8-10 MPa
- **Replace threshold**: ~3-5 MPa
- **Empty**: <1 MPa

---

**Last Updated**: 2026-02-09  
**Skill ID**: 182  
**Version**: 2.0 (K-Dense Standard)
