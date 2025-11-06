# Estimation Techniques

Comprehensive methods for realistic time and effort estimation.

## Estimation Methods

### Bottom-Up Estimation

**Most accurate method**: Estimate each small task, sum to get total.

**Process**:
1. Break project into small tasks (1-3 days each)
2. Estimate each task individually
3. Sum task estimates
4. Add buffers at task, phase, and project levels

**Example**:
```
Task A: 2 days
Task B: 3 days
Task C: 1 day
Subtotal: 6 days
Buffer (20%): 1.2 days
Total: 7.2 days ≈ 8 days
```

**Best for**: Well-understood work, detailed planning phase

**Pros**: Most accurate, catches hidden work
**Cons**: Time-consuming, requires detail

### Top-Down Estimation

**Faster method**: Estimate overall, allocate to parts.

**Process**:
1. Estimate total project duration
2. Allocate percentage to each phase
3. Validate against constraints

**Example**:
```
Total project: 40 days

Allocation:
- Planning: 10% (4 days)
- Design: 20% (8 days)
- Build: 50% (20 days)
- Test: 15% (6 days)
- Deploy: 5% (2 days)
Total: 40 days
```

**Best for**: Early planning, similar to past projects

**Pros**: Quick, good for high-level planning
**Cons**: Less accurate, may miss details

### Three-Point Estimation

**Accounts for uncertainty** using optimistic, likely, and pessimistic estimates.

**Formula**: Expected = (O + 4M + P) / 6

Where:
- O = Optimistic (best case)
- M = Most Likely (realistic)
- P = Pessimistic (worst case)

**Example**:
```
Task: Build user authentication

O (best case): 3 days
M (likely): 5 days
P (worst case): 9 days

Expected = (3 + 4×5 + 9) / 6 = 32 / 6 = 5.3 days ≈ 6 days
```

**Best for**: Uncertain work, novel tasks, risk-aware planning

**Standard Deviation** (uncertainty measure): (P - O) / 6
Example: (9 - 3) / 6 = 1 day uncertainty

### Analogous Estimation

**Compare to similar past projects**.

**Process**:
1. Find similar completed project
2. Identify differences (size, complexity, team)
3. Adjust estimate based on differences

**Example**:
```
Past project: E-commerce site (120 days)
New project: E-commerce site with mobile app

Adjustments:
- Mobile app adds 30%: +36 days
- Experienced team: -10%: -12 days
- Estimate: 120 + 36 - 12 = 144 days
```

**Best for**: Similar work, have historical data

**Pros**: Quick, leverages experience
**Cons**: Requires comparable projects, may miss unique aspects

### Parametric Estimation

**Use statistical relationships** between variables.

**Examples**:
- Lines of code × time per line
- Square footage × time per sqft
- Features × time per feature

**Example**:
```
Historical data: Mobile apps average 2 days per screen

New app: 40 screens
Base estimate: 40 × 2 = 80 days
Complexity adjustment: +20% = 96 days
```

**Best for**: Repetitive work, have metrics

**Pros**: Data-driven, scalable
**Cons**: Requires metrics, may oversimplify

## Buffer Strategies

### Task-Level Buffers

Add contingency to individual estimates:
- **Routine work**: 10-15% buffer
- **Familiar work**: 15-20% buffer
- **Uncertain work**: 20-30% buffer
- **Novel work**: 30-50% buffer

### Project-Level Buffers

**Management Reserve**: 10-20% of total for unknowns

**Contingency Reserve**: For known risks (based on risk assessment)

**Buffer Placement**:
- At milestones (not scattered)
- At integration points
- At external dependencies
- At end of phases

## Effort vs Duration

**Key Distinction**:
- **Effort**: Total work hours (person-hours)
- **Duration**: Calendar time elapsed

**Example**:
- Effort: 40 hours of work
- Duration: 2 weeks (person working 4 hours/day)

**Factors affecting duration**:
- Team availability (full-time vs part-time)
- Productivity rate (not 8 productive hours/day)
- Dependencies (waiting time)
- Context switching
- Meetings and interruptions

**Typical productivity**:
- 6 productive hours per 8-hour day
- 5 productive hours when shared across projects
- 4 productive hours with high meeting load

## Common Estimation Biases

### Planning Fallacy

**Tendency to underestimate** time, costs, risks.

**Causes**:
- Optimism bias
- Anchoring on best-case
- Ignoring past performance
- Not accounting for unknowns

**Mitigation**:
- Use reference class forecasting (historical data)
- Three-point estimation
- Add explicit buffers
- Review past estimates vs actuals

### Anchoring Bias

**First number** heavily influences subsequent estimates.

**Example**: "Should take about 2 weeks..." then all estimates cluster near 2 weeks

**Mitigation**:
- Estimate independently first
- Use multiple methods
- Question initial estimates
- Compare to historical data

### Student Syndrome

**Work expands to fill time available**.

**Mitigation**:
- Use timeboxing
- Regular check-ins
- Milestone-based tracking
- Reduce artificial padding

## Velocity and Capacity

### Team Velocity

**Average amount of work completed per iteration**.

**Calculation**:
```
Sprint 1: 23 points
Sprint 2: 27 points
Sprint 3: 25 points

Average velocity: (23 + 27 + 25) / 3 = 25 points
```

**Use**: Plan future sprints based on actual velocity

**Stabilization**: Velocity stabilizes after 3-4 sprints

### Capacity Planning

**Available hours per iteration**:
```
Team: 6 people
Sprint: 2 weeks = 10 days
Hours/day: 6 productive hours

Gross capacity: 6 × 10 × 6 = 360 hours

Subtract:
- Sprint ceremonies: 12 hours
- Other meetings: 20 hours
- Training: 8 hours
- Support: 10 hours
- PTO: 16 hours

Net capacity: 360 - 66 = 294 hours
```

**Plan conservatively**: Use 80-90% of net capacity

## Estimation Techniques by Project Type

### Software Development

**Methods**: Story points, t-shirt sizing, planning poker

**Factors**:
- Complexity
- Uncertainty
- Dependencies
- Technical debt

**Typical velocity**: 15-40 story points per sprint (2 weeks)

### Construction/Physical Projects

**Methods**: Parametric (sqft, units), analogous

**Factors**:
- Materials
- Weather delays
- Inspections/permits
- Contractor availability

**Buffers**: 20-30% for weather and delays

### Research Projects

**Methods**: Phased approach with reviews

**Factors**:
- Exploratory nature
- Literature review time
- IRB/approval processes
- Data collection challenges

**Buffers**: 30-50% for unknown unknowns

### Event Planning

**Methods**: Backward planning from fixed date

**Factors**:
- Vendor lead times
- Approval processes
- Rehearsal time
- Setup/teardown

**Critical**: Account for procurement delays (30-90 days for some vendors)

## Improving Estimation Accuracy

### Track Actuals vs Estimates

**For each task, record**:
- Original estimate
- Actual time spent
- Variance
- Reason for variance

**Analyze patterns**:
- Consistently over/under?
- Which types of tasks?
- Common missed factors?

### Estimation Workshops

**Team estimation**:
- Planning poker
- T-shirt sizing
- Delphi method (anonymous rounds)

**Benefits**:
- Diverse perspectives
- Knowledge sharing
- Better buy-in
- Catches blind spots

### Reference Class Forecasting

**Use historical data** from similar projects:
1. Identify reference class (similar projects)
2. Collect actual durations
3. Use distribution for estimation
4. Adjust for unique factors

**Example**:
Past 10 mobile apps: 60-180 days, average 110 days
New app estimate: 100-120 days (toward lower end due to experience)

---

## Quick Reference

| Situation | Best Method | Buffer |
|-----------|------------|---------|
| Detailed plan | Bottom-up | 20-30% |
| Early planning | Top-down | 30-40% |
| Uncertain work | Three-point | 30-50% |
| Similar to past | Analogous | 15-25% |
| Repetitive work | Parametric | 10-20% |
| Novel/research | Phased + reviews | 40-60% |
| Fixed deadline | Backward planning | 25-35% |

**Remember**: 
- All estimates are guesses - build in buffers
- Accuracy improves with detail and experience
- Track actuals to improve future estimates
- Include team in estimation process
