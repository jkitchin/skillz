# Planning Frameworks - Detailed Guide

Comprehensive step-by-step guides for each planning methodology.

## Table of Contents

- [Work Breakdown Structure (WBS)](#work-breakdown-structure-wbs)
- [Backward Planning](#backward-planning)
- [Critical Path Method (CPM)](#critical-path-method-cpm)
- [Agile/Iterative Planning](#agileiterative-planning)
- [Phased/Milestone Planning](#phasedmilestone-planning)
- [Risk Management](#risk-management)
- [Resource Planning](#resource-planning)

---

## Work Breakdown Structure (WBS)

### Overview

Hierarchical decomposition of project into progressively smaller, manageable components.

### Step-by-Step Process

#### 1. Define Project Scope (10 minutes)

**Write clear project statement**:
- What: The deliverable or outcome
- Why: The purpose and value
- Boundaries: What's included and excluded

Example: "Website redesign project to improve user experience and increase conversions by 20%, excluding backend infrastructure changes"

#### 2. Identify Major Phases (15 minutes)

**Break project into 3-7 top-level phases**:

Typical phases might include:
- Initiation/Planning
- Design
- Development/Build
- Testing/QA
- Deployment/Launch
- Closure/Handoff

**Questions to ask**:
- What are the natural stages of work?
- What major milestones exist?
- How does work flow chronologically?

#### 3. Decompose Phases into Deliverables (30 minutes)

**For each phase, identify concrete deliverables**:

Example for "Design" phase:
- User research findings
- Wireframes
- Visual mockups
- Design system documentation
- Prototype

**Use nouns** for deliverables (things produced)

**100% Rule**: All deliverables should represent 100% of phase output

#### 4. Break Deliverables into Tasks (45 minutes)

**For each deliverable, list tasks required**:

Example for "Visual mockups" deliverable:
- Review brand guidelines
- Create homepage mockup
- Create product page mockups
- Create checkout flow mockups
- Review with stakeholders
- Incorporate feedback
- Finalize designs

**Use verbs** for tasks (actions taken)

**Size guideline**: Tasks should be 1-3 days of work

#### 5. Add Sub-tasks if Needed (optional, 20 minutes)

**Only if tasks are still too large**:

Example: "Create homepage mockup" might become:
- Sketch layout options
- Select best layout
- Design header section
- Design hero section
- Design features section
- Design footer section

**Stop when** tasks are manageable (1-3 days)

**Avoid** going more than 4-5 levels deep

#### 6. Validate Completeness (15 minutes)

**Check each level**:
- Does this represent 100% of parent's work?
- Are sibling items mutually exclusive?
- Is anything missing?
- Is anything duplicated?

**Cross-check** with stakeholders

### WBS Formats

#### Hierarchical List
```
1.0 Website Redesign
  1.1 Planning Phase
    1.1.1 Define requirements
    1.1.2 Create project charter
    1.1.3 Identify stakeholders
  1.2 Design Phase
    1.2.1 User Research
      1.2.1.1 Conduct user interviews
      1.2.1.2 Analyze current analytics
      1.2.1.3 Create user personas
    1.2.2 Wireframes
      1.2.2.1 Sketch initial concepts
      1.2.2.2 Create digital wireframes
      1.2.2.3 User testing of wireframes
```

#### Outline Format
```
Website Redesign
├── Planning
│   ├── Requirements definition
│   ├── Project charter
│   └── Stakeholder identification
├── Design
│   ├── User Research
│   │   ├── User interviews
│   │   ├── Analytics analysis
│   │   └── Persona creation
│   └── Wireframes
│       ├── Concept sketches
│       ├── Digital wireframes
│       └── User testing
```

### Tips for Effective WBS

**Do**:
- Start top-down (big picture first)
- Use consistent level of detail
- Make items mutually exclusive
- Ensure completeness at each level
- Use team input

**Don't**:
- Mix tasks and deliverables at same level
- Go too deep (usually 3-4 levels enough)
- Include dependencies in WBS (save for scheduling)
- Assign durations yet (comes later)

---

## Backward Planning

### Overview

Start from goal/deadline and work backwards to identify all prerequisite steps.

### Step-by-Step Process

#### 1. Define End State (10 minutes)

**Be specific about goal**:
- What: Exactly what's achieved
- When: Specific date and time
- Success criteria: How you'll know it's done
- Deliverable: What exists at completion

Example: "Product launch event on June 15, 2024 at 2pm ET with 500 registered attendees, demo ready, press kit available"

#### 2. Identify Immediate Prerequisites (10 minutes)

**Ask**: "What must be true/completed the moment before this happens?"

For product launch:
- Venue confirmed and set up
- Attendee registrations complete
- Product demo functioning
- Marketing materials printed
- Press kit distributed
- Team rehearsed and ready

#### 3. Work Backwards for Each Prerequisite (30-60 minutes)

**For each item, ask**: "What must happen before this?"

Example: "Product demo functioning" requires:
- Demo script finalized
- Demo environment configured
- Demo data prepared
- Demo tested successfully

Continue backwards:
- "Demo script finalized" requires:
  - Product features identified
  - Story arc created
  - Script written
  - Script reviewed

Keep going until you reach present day or tasks you can start now.

#### 4. Build Dependency Network (20 minutes)

**Identify which prerequisites**:
- Must be sequential
- Can run in parallel
- Have cross-dependencies

Example:
- "Venue setup" depends on "Venue confirmed"
- "Attendee registrations" can run parallel to "Demo preparation"
- "Press kit printing" depends on "Marketing copy finalized"

#### 5. Estimate Durations (30 minutes)

**For each task, estimate time needed**:
- Include buffer (20-30%)
- Consider dependencies
- Account for approvals/reviews
- Factor in procurement lead times

Working backwards from June 15:
- Venue setup: 1 day before (June 14)
- Venue confirmed: 60 days before (April 16)
- Venue options identified: 75 days before (April 1)

#### 6. Calculate Start Dates (15 minutes)

**Work backwards from deadline**:
- Latest completion date for each task
- Buffer for risks
- Identify earliest start date

If June 15 is launch and venue needs 60 days notice:
- Latest venue confirmation: April 16
- Add buffer (7 days): April 9
- Start venue search: March 1 (allowing 6 weeks to choose)

#### 7. Create Forward Timeline (20 minutes)

**Reverse the sequence**:
- Convert "days before deadline" to actual dates
- Create forward-looking schedule
- Add milestones
- Assign owners

### Backward Planning Template

```
GOAL: [Specific end state]
DATE: [Exact deadline]

WORKING BACKWARDS:

T-0 days (Deadline):
- [Final state achieved]

T-1 day:
- [What must be ready day before]

T-1 week:
- [What must be ready week before]

T-2 weeks:
- [Prerequisites for T-1 week]

T-1 month:
- [Prerequisites for T-2 weeks]

Continue backwards to present...

FORWARD TIMELINE:
[Today] → [Milestone 1] → [Milestone 2] → [Goal]
```

### Common Pitfalls

**Missing prerequisites**:
- Review/approval time
- Procurement/shipping delays
- Learning/training time
- Setup/configuration time
- Handoff/transition time

**Insufficient buffers**:
- External dependencies (vendors, approvals)
- Novel/uncertain work
- Critical path tasks

---

## Critical Path Method (CPM)

### Overview

Identify sequence of tasks that determines minimum project duration and focus management attention.

### Key Concepts

**Critical Path**: Longest sequence of dependent tasks; determines minimum project duration

**Float/Slack**: Amount of time a task can be delayed without delaying project

**Critical Tasks**: Tasks on critical path (zero float); any delay delays entire project

**Non-Critical Tasks**: Tasks with float; can be delayed without affecting project

### Calculation Steps

#### 1. List All Tasks with Durations (15 minutes)

```
Task A: Design wireframes (5 days)
Task B: Client approval (2 days)
Task C: Create visual designs (7 days)
Task D: Develop frontend (10 days)
Task E: Develop backend (8 days)
Task F: Integration testing (3 days)
Task G: User acceptance testing (4 days)
Task H: Deploy (1 day)
```

#### 2. Identify Dependencies (15 minutes)

```
Task B depends on Task A (must have wireframes to approve)
Task C depends on Task B (need approval before designs)
Task D depends on Task C (need designs for frontend)
Task E depends on Task C (need designs for backend)
Task F depends on Task D and E (both must be complete)
Task G depends on Task F (integration must pass first)
Task H depends on Task G (UAT must pass first)
```

#### 3. Calculate Early Start/Finish (Forward Pass) (20 minutes)

**Early Start (ES)**: Earliest a task can start
**Early Finish (EF)**: ES + Duration

```
Task A: ES=0, EF=5
Task B: ES=5 (after A), EF=7
Task C: ES=7 (after B), EF=14
Task D: ES=14 (after C), EF=24
Task E: ES=14 (after C), EF=22
Task F: ES=24 (after D and E, take latest), EF=27
Task G: ES=27 (after F), EF=31
Task H: ES=31 (after G), EF=32
```

**Project Duration**: 32 days

#### 4. Calculate Late Start/Finish (Backward Pass) (20 minutes)

**Late Finish (LF)**: Latest a task can finish without delaying project
**Late Start (LS)**: LF - Duration

Working backwards from day 32:

```
Task H: LF=32, LS=31
Task G: LF=31, LS=27
Task F: LF=27, LS=24
Task E: LF=24, LS=16
Task D: LF=24, LS=14
Task C: LF=14, LS=7
Task B: LF=7, LS=5
Task A: LF=5, LS=0
```

#### 5. Calculate Float/Slack (10 minutes)

**Float = LS - ES** (or LF - EF)

```
Task A: 0-0 = 0 (CRITICAL)
Task B: 5-5 = 0 (CRITICAL)
Task C: 7-7 = 0 (CRITICAL)
Task D: 14-14 = 0 (CRITICAL)
Task E: 16-14 = 2 days float
Task F: 24-24 = 0 (CRITICAL)
Task G: 27-27 = 0 (CRITICAL)
Task H: 31-31 = 0 (CRITICAL)
```

#### 6. Identify Critical Path (5 minutes)

Tasks with zero float form the critical path:

**Critical Path**: A → B → C → D → F → G → H (32 days)

**Non-Critical**: Task E (has 2 days float)

### Using Critical Path Information

**Focus Management**:
- Monitor critical tasks closely
- Prioritize resources for critical path
- Can't delay critical tasks without delaying project

**Optimization**:
- Shorten critical tasks to shorten project
- Shift resources from non-critical to critical tasks
- Add parallelization where possible

**Risk Management**:
- Buffer critical tasks
- Have backup plans for critical dependencies
- Monitor critical path changes

### Automated Calculation

Use `scripts/critical_path.py` for automatic calculation:

```bash
python scripts/critical_path.py tasks.json
```

Input format (JSON):
```json
{
  "tasks": [
    {"id": "A", "name": "Design wireframes", "duration": 5, "dependencies": []},
    {"id": "B", "name": "Client approval", "duration": 2, "dependencies": ["A"]},
    {"id": "C", "name": "Visual designs", "duration": 7, "dependencies": ["B"]},
    ...
  ]
}
```

---

## Agile/Iterative Planning

### Sprint Planning

#### Preparation (Before Meeting)

**Product Owner**:
- Prioritize backlog
- Ensure top items are "ready" (clear, estimated, valuable)
- Define sprint goal candidate

**Team**:
- Review backlog
- Identify questions/concerns
- Understand capacity

#### Sprint Planning Meeting (2-4 hours for 2-week sprint)

**Part 1: What will we deliver? (1-2 hours)**

1. Review sprint goal
2. Review team capacity
3. Select backlog items
4. Confirm commitment

**Part 2: How will we do the work? (1-2 hours)**

1. Break items into tasks
2. Identify dependencies
3. Assign initial ownership
4. Identify risks/blockers

**Output**:
- Sprint goal
- Selected backlog items
- Task breakdown
- Team commitment

#### Capacity Planning

**Calculate available hours**:
```
Team size: 5 people
Sprint length: 10 days
Hours per day: 6 (accounting for meetings, breaks)

Gross capacity: 5 × 10 × 6 = 300 hours

Subtract:
- Ceremonies: 8 hours
- Other meetings: 10 hours
- PTO/holidays: 16 hours

Net capacity: 266 hours
```

**Apply velocity**:
If past 3 sprints completed 50, 48, 52 story points:
- Average velocity: 50 points
- Plan for ~50 points this sprint

#### Backlog Refinement

**Regular grooming (weekly)**:
- Clarify requirements
- Add acceptance criteria
- Estimate effort
- Split large items
- Remove obsolete items

**Story Splitting**:
Large stories split by:
- User roles (admin vs regular user)
- CRUD operations (create first, then read/update/delete)
- Happy path vs edge cases
- Simple vs complex variations

#### Daily Standup (15 minutes)

**Each person answers**:
1. What did I complete yesterday?
2. What will I do today?
3. Any blockers?

**Purpose**: Coordination, not status reporting

**Scrum Master**: Remove blockers, facilitate

#### Sprint Review (1-2 hours)

**Demonstrate completed work**:
- Show working software
- Get feedback
- Discuss what's done vs not done

**Not a status meeting**: Live demo, working software

#### Sprint Retrospective (1-1.5 hours)

**Discuss**:
- What went well?
- What could be improved?
- What will we commit to changing?

**Output**: 1-3 concrete improvements for next sprint

### Kanban Planning

**Columns** (typical):
- Backlog → To Do → In Progress → Review → Done

**WIP Limits**: Maximum items in each column
- Prevents overload
- Highlights bottlenecks
- Improves flow

**Continuous planning**:
- Pull from backlog when capacity available
- Reprioritize anytime
- No fixed iterations

---

## Phased/Milestone Planning

### Phase Definition

**Phases** represent distinct stages with clear boundaries:

**Typical phases**:
1. Initiation
2. Planning
3. Execution
4. Monitoring & Control (continuous)
5. Closure

**Or domain-specific**:
- Research: Literature review → Methodology → Data collection → Analysis → Writing
- Construction: Design → Permits → Foundation → Structure → Finishes → Handover
- Software: Discovery → Design → Development → Testing → Deployment → Maintenance

**Phase gates**: Review/approval points between phases

### Milestone Planning

**Good milestones**:
- SMART (Specific, Measurable, Achievable, Relevant, Time-bound)
- Binary (clearly done or not done)
- Meaningful (represent significant progress)
- Visible (stakeholders can see)

**Milestone types**:
- **Deliverable**: "Website design approved"
- **Decision**: "Go/no-go decision made"
- **Event**: "Product launched"
- **Gate**: "Phase 1 completed"

**Spacing guidelines**:
- Short projects (<1 month): Weekly
- Medium projects (1-6 months): Bi-weekly or monthly
- Long projects (>6 months): Monthly

### Stage-Gate Process

**At each gate, review**:
- Deliverables complete?
- Quality criteria met?
- Risks acceptable?
- Budget on track?
- Continue to next phase?

**Possible outcomes**:
- Proceed to next phase
- Conditional approval (with fixes)
- Recycle to current phase
- Hold for more information
- Cancel project

---

## Risk Management

### Risk Register

Track all identified risks:

| Risk | Category | Probability | Impact | Score | Response | Owner | Status |
|------|----------|-------------|--------|-------|----------|-------|--------|
| Key developer leaves | Resource | 3 | 4 | 12 | Knowledge sharing, backup | PM | Active |
| Vendor delays | External | 4 | 3 | 12 | Early ordering, backup vendor | Ops | Active |
| Scope creep | Scope | 4 | 4 | 16 | Change control process | PM | Mitigated |

**Probability & Impact**: 1-5 scale
**Risk Score**: Probability × Impact
**Priority**: High (12+), Medium (6-11), Low (1-5)

### Response Planning

**For each high-priority risk**:
1. Prevention actions (reduce probability)
2. Mitigation actions (reduce impact)
3. Contingency plan (if it occurs)
4. Trigger indicators (how to know it's happening)
5. Owner (who monitors and responds)

See `references/templates.md` for risk register template.

---

## Resource Planning

### RACI Matrix

Define who does what:

| Task/Decision | Team Lead | Designer | Developer | QA | Client |
|--------------|-----------|----------|-----------|----|----|
| Requirements | C | C | I | I | A |
| Design | A | R | C | I | C |
| Development | A | I | R | C | I |
| Testing | C | I | C | R | I |
| Approval | I | I | I | I | A |

**R**: Responsible (does the work) - can be multiple
**A**: Accountable (ultimately answers) - only one
**C**: Consulted (provides input)
**I**: Informed (kept in loop)

See `references/templates.md` for RACI template.

---

