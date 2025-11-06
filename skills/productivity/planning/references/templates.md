# Planning Templates

Ready-to-use templates for common planning artifacts. Copy, customize, and adapt to your needs.

## Table of Contents

- [Project Charter](#project-charter)
- [Work Breakdown Structure (WBS)](#work-breakdown-structure-wbs)
- [Risk Register](#risk-register)
- [RACI Matrix](#raci-matrix)
- [Project Plan (One-Pager)](#project-plan-one-pager)
- [Sprint Planning Template](#sprint-planning-template)
- [Retrospective Template](#retrospective-template)
- [Status Report Template](#status-report-template)
- [Milestone Tracker](#milestone-tracker)
- [Gantt Chart Structure](#gantt-chart-structure)

---

## Project Charter

Use this to formally initiate a project and align stakeholders.

```markdown
# Project Charter: [Project Name]

**Date:** [Date]
**Project Manager:** [Name]
**Sponsor:** [Name]

## Project Purpose

Why this project exists and what problem it solves:

[2-3 sentences on business justification]

## Objectives

What we aim to achieve (SMART goals):

1. [Objective 1 - Specific, Measurable, Achievable, Relevant, Time-bound]
2. [Objective 2]
3. [Objective 3]

## Scope

**In Scope:**
- [Deliverable/feature 1]
- [Deliverable/feature 2]
- [Deliverable/feature 3]

**Out of Scope:**
- [What we're explicitly NOT doing]
- [Boundary clarification]

## Success Criteria

How we'll measure success:

- [ ] [Criterion 1 - quantifiable]
- [ ] [Criterion 2 - quantifiable]
- [ ] [Criterion 3 - quantifiable]

## Key Stakeholders

| Name | Role | Responsibility | Contact |
|------|------|----------------|---------|
| [Name] | Sponsor | Approvals, resources | [email] |
| [Name] | Project Manager | Delivery | [email] |
| [Name] | Technical Lead | Architecture | [email] |
| [Name] | End User Rep | Requirements | [email] |

## Timeline

**Start Date:** [Date]
**End Date:** [Date]
**Duration:** [X weeks/months]

**Key Milestones:**
- [Date]: [Milestone 1]
- [Date]: [Milestone 2]
- [Date]: [Milestone 3]

## Budget

**Total Budget:** $[Amount]

**Breakdown:**
- Personnel: $[Amount]
- Tools/Software: $[Amount]
- External Services: $[Amount]
- Contingency (15%): $[Amount]

## High-Level Risks

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| [Risk 1] | High/Med/Low | High/Med/Low | [Strategy] |
| [Risk 2] | High/Med/Low | High/Med/Low | [Strategy] |

## Assumptions

- [Assumption 1]
- [Assumption 2]
- [Assumption 3]

## Constraints

- [Constraint 1 - time, budget, resource, technical]
- [Constraint 2]

## Approvals

| Name | Role | Signature | Date |
|------|------|-----------|------|
| [Name] | Sponsor | | |
| [Name] | Stakeholder | | |

---

**Notes:**

[Any additional context or information]
```

---

## Work Breakdown Structure (WBS)

Hierarchical decomposition of project work.

### Outline Format

```
1.0 [Project Name]
  1.1 [Phase 1 Name]
    1.1.1 [Deliverable 1]
      1.1.1.1 [Task 1]
      1.1.1.2 [Task 2]
      1.1.1.3 [Task 3]
    1.1.2 [Deliverable 2]
      1.1.2.1 [Task 1]
      1.1.2.2 [Task 2]
  1.2 [Phase 2 Name]
    1.2.1 [Deliverable 1]
      1.2.1.1 [Task 1]
      1.2.1.2 [Task 2]
    1.2.2 [Deliverable 2]
      1.2.2.1 [Task 1]
      1.2.2.2 [Task 2]
      1.2.2.3 [Task 3]
  1.3 [Phase 3 Name]
    1.3.1 [Deliverable 1]
    1.3.2 [Deliverable 2]
```

### Task Detail Template

```markdown
**Task ID:** 1.2.1.1
**Task Name:** [Descriptive name]
**Description:** [What needs to be done]
**Owner:** [Person responsible]
**Duration:** [X days/hours]
**Dependencies:** [Task IDs this depends on]
**Deliverable:** [Output/artifact produced]
**Status:** Not Started / In Progress / Complete
**Notes:** [Additional context]
```

---

## Risk Register

Track and manage project risks.

### Markdown Table Format

```markdown
# Risk Register

**Project:** [Project Name]
**Last Updated:** [Date]

## Risk Scoring
- **Probability:** 1=Very Low, 2=Low, 3=Medium, 4=High, 5=Very High
- **Impact:** 1=Negligible, 2=Minor, 3=Moderate, 4=Significant, 5=Severe
- **Risk Score:** Probability × Impact
- **Priority:** High (15-25), Medium (6-14), Low (1-5)

| ID | Risk Description | Category | Prob | Impact | Score | Priority | Response Strategy | Mitigation Actions | Owner | Status | Trigger Indicators |
|----|-----------------|----------|------|--------|-------|----------|-------------------|-------------------|-------|--------|-------------------|
| R01 | [Describe the risk] | Schedule | 4 | 4 | 16 | High | Mitigate | [Actions to reduce probability/impact] | [Name] | Active | [How to know it's happening] |
| R02 | [Describe the risk] | Technical | 3 | 5 | 15 | High | Mitigate | [Actions] | [Name] | Active | [Triggers] |
| R03 | [Describe the risk] | Resource | 3 | 3 | 9 | Medium | Accept | [Contingency plan] | [Name] | Monitoring | [Triggers] |
| R04 | [Describe the risk] | External | 2 | 4 | 8 | Medium | Transfer | [Contract terms, insurance] | [Name] | Active | [Triggers] |

**Categories:** Schedule, Technical, Resource, External, Scope, Quality, Budget

**Response Strategies:**
- **Avoid:** Eliminate the risk by changing approach
- **Mitigate:** Reduce probability or impact
- **Transfer:** Shift to another party (vendor, insurance)
- **Accept:** Monitor with contingency plan

**Status:** Active (mitigating), Monitoring (watching), Occurred (contingency activated), Closed
```

### CSV Format

For spreadsheet import:

```csv
ID,Risk Description,Category,Probability,Impact,Score,Priority,Response Strategy,Mitigation Actions,Owner,Status,Trigger Indicators
R01,Key developer leaves,Resource,3,4,12,High,Mitigate,"Knowledge sharing sessions, documentation, pair programming",PM,Active,"Decreased engagement, job search signals"
R02,Third-party API changes,External,4,3,12,Medium,Mitigate,"API versioning, fallback service, monitoring",Tech Lead,Active,"Deprecation notices, vendor announcements"
R03,Scope creep,Scope,4,4,16,High,Avoid,"Change control process, stakeholder alignment",PM,Active,"Unplanned requests, expanding requirements"
```

---

## RACI Matrix

Define roles and responsibilities.

```markdown
# RACI Matrix

**Project:** [Project Name]

**Legend:**
- **R** = Responsible (does the work)
- **A** = Accountable (ultimately answerable, only one per task)
- **C** = Consulted (provides input)
- **I** = Informed (kept in the loop)

| Task/Decision | Project Manager | Tech Lead | Designer | Developer | QA | Stakeholder | Client |
|--------------|----------------|-----------|----------|-----------|-------|-------------|--------|
| Project charter | A | C | I | I | I | C | C |
| Requirements gathering | A | C | C | C | I | R | A |
| Technical design | C | A | I | R | C | I | I |
| UI/UX design | C | I | A/R | C | I | C | C |
| Development | C | A | I | R | C | I | I |
| Code review | I | A | I | C | I | I | I |
| Testing plan | C | C | I | C | A/R | I | I |
| Testing execution | C | C | I | C | A/R | I | I |
| Bug fixes | C | A | C | R | C | I | I |
| Deployment | A | R | I | R | C | I | I |
| User training | R | I | C | I | I | C | C |
| Go-live approval | C | C | I | I | I | C | A |
| Post-launch support | A | R | I | R | C | I | C |

**Notes:**
- Each row should have exactly one "A" (Accountable)
- "R" (Responsible) can be multiple people
- Avoid having too many "C" - consultation should be selective
```

---

## Project Plan (One-Pager)

High-level overview for stakeholders.

```markdown
# Project Plan: [Project Name]

**Date:** [Date] | **Version:** 1.0 | **Status:** [Planning/Active/On Hold]

---

## Overview

**Goal:** [One sentence - what are we achieving]

**Why:** [One sentence - business value]

**Timeline:** [Start Date] → [End Date] ([X weeks/months])

**Budget:** $[Amount]

**Owner:** [Name]

---

## Key Deliverables

1. **[Deliverable 1]** - [Brief description]
2. **[Deliverable 2]** - [Brief description]
3. **[Deliverable 3]** - [Brief description]

---

## Phases & Milestones

```
Phase 1: [Name] (Weeks 1-3)
├─ [Milestone 1] - [Date]
└─ [Deliverable]

Phase 2: [Name] (Weeks 4-7)
├─ [Milestone 2] - [Date]
└─ [Deliverable]

Phase 3: [Name] (Weeks 8-10)
├─ [Milestone 3] - [Date]
└─ [Deliverable]
```

---

## Success Metrics

- [ ] [Quantifiable metric 1]
- [ ] [Quantifiable metric 2]
- [ ] [Quantifiable metric 3]

---

## Top Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| [Risk 1] | High | [How we're addressing it] |
| [Risk 2] | Medium | [How we're addressing it] |

---

## Team

- **Project Manager:** [Name]
- **Tech Lead:** [Name]
- **Team Members:** [Names]
- **Stakeholders:** [Names]

---

## Next Steps

1. [ ] [Action 1] - [Owner] - [Due Date]
2. [ ] [Action 2] - [Owner] - [Due Date]
3. [ ] [Action 3] - [Owner] - [Due Date]

---

**Questions or concerns?** Contact [Name] at [email]
```

---

## Sprint Planning Template

For agile/iterative planning (2-week sprint example).

```markdown
# Sprint Planning: Sprint [#]

**Date:** [Date]
**Sprint Duration:** [Start Date] → [End Date] (2 weeks)
**Team Capacity:** [X hours] (accounting for PTO, meetings)

---

## Sprint Goal

[One sentence describing what this sprint will achieve]

**Success Criteria:**
- [ ] [Criterion 1]
- [ ] [Criterion 2]

---

## Team Capacity

| Team Member | Available Days | Hours/Day | Total Hours | Notes |
|-------------|---------------|-----------|-------------|-------|
| [Name] | 10 | 6 | 60 | |
| [Name] | 8 | 6 | 48 | 2 days PTO |
| [Name] | 10 | 6 | 60 | |
| **Total** | | | **168** | |

**Capacity Adjustments:**
- Ceremonies: -8 hours
- Other meetings: -10 hours
- **Net Capacity: 150 hours**

---

## Sprint Backlog

### Committed Items

| ID | Story/Task | Priority | Estimate | Owner | Status |
|----|-----------|----------|----------|-------|--------|
| US-101 | [User story description] | High | 8h | [Name] | Not Started |
| US-102 | [User story description] | High | 13h | [Name] | Not Started |
| US-103 | [User story description] | Medium | 5h | [Name] | Not Started |
| TASK-201 | [Technical task] | High | 8h | [Name] | Not Started |

**Total Committed:** [X hours / Y story points]

### Stretch Goals

If we finish early:

| ID | Story/Task | Estimate | Notes |
|----|-----------|----------|-------|
| US-104 | [Lower priority item] | 5h | Nice to have |

---

## Task Breakdown

### US-101: [Story Name]

**Acceptance Criteria:**
- [ ] [Criterion 1]
- [ ] [Criterion 2]
- [ ] [Criterion 3]

**Tasks:**
- [ ] [Task 1] - [Estimate] - [Owner]
- [ ] [Task 2] - [Estimate] - [Owner]
- [ ] [Task 3] - [Estimate] - [Owner]

**Dependencies:** [Any blockers or prerequisites]

---

## Definition of Done

For this sprint, work is done when:

- [ ] Code complete and committed
- [ ] Unit tests written and passing
- [ ] Code reviewed and approved
- [ ] Integration tested
- [ ] Documentation updated
- [ ] Acceptance criteria met
- [ ] Deployed to staging
- [ ] Demo-ready

---

## Risks & Blockers

| Risk/Blocker | Impact | Mitigation | Owner |
|--------------|--------|------------|-------|
| [Issue 1] | High | [Plan] | [Name] |
| [Issue 2] | Medium | [Plan] | [Name] |

---

## Daily Standup Schedule

**Time:** [Time] daily
**Duration:** 15 minutes
**Format:** What I did, what I'll do, blockers

---

## Sprint Ceremonies

- **Sprint Planning:** [Date/Time] - 2 hours
- **Daily Standup:** [Time] daily - 15 min
- **Sprint Review:** [Date/Time] - 1 hour
- **Retrospective:** [Date/Time] - 1 hour

---

**Notes:**

[Any additional context or decisions made during planning]
```

---

## Retrospective Template

Continuous improvement session.

```markdown
# Sprint Retrospective: Sprint [#]

**Date:** [Date]
**Participants:** [Names]
**Facilitator:** [Name]

---

## Sprint Metrics

**Committed:** [X story points / hours]
**Completed:** [Y story points / hours]
**Completion Rate:** [Y/X %]

**Velocity Trend:**
- Sprint N-2: [Z points]
- Sprint N-1: [Z points]
- This Sprint: [Y points]

---

## What Went Well? 👍

Things to celebrate and continue:

1. [Positive observation 1]
   - Impact: [How this helped]
   - Why: [Root cause of success]

2. [Positive observation 2]
   - Impact: [How this helped]
   - Why: [Root cause of success]

3. [Positive observation 3]
   - Impact: [How this helped]
   - Why: [Root cause of success]

---

## What Could Be Improved? 🤔

Things that slowed us down or caused frustration:

1. [Challenge 1]
   - Impact: [How this hurt]
   - Frequency: [How often this happened]
   - Root cause: [Why this happened]

2. [Challenge 2]
   - Impact: [How this hurt]
   - Frequency: [How often this happened]
   - Root cause: [Why this happened]

3. [Challenge 3]
   - Impact: [How this hurt]
   - Frequency: [How often this happened]
   - Root cause: [Why this happened]

---

## Action Items

**Commit to 1-3 concrete improvements for next sprint:**

| Action | Owner | Deadline | Success Measure |
|--------|-------|----------|----------------|
| [Specific action 1] | [Name] | [Date] | [How we'll know it worked] |
| [Specific action 2] | [Name] | [Date] | [How we'll know it worked] |
| [Specific action 3] | [Name] | [Date] | [How we'll know it worked] |

**Follow-up from last retro:**

| Previous Action | Status | Notes |
|----------------|--------|-------|
| [Action from last sprint] | ✅ Done / 🔄 In Progress / ❌ Not Done | [Outcome] |

---

## Appreciations 💙

Shout-outs to team members:

- **[Name]:** [What they did and why it mattered]
- **[Name]:** [What they did and why it mattered]

---

## Insights & Themes

Recurring patterns or observations:

- [Pattern 1]
- [Pattern 2]

---

**Retro Format Used:** [Start/Stop/Continue, Mad/Sad/Glad, 4Ls, etc.]

**Next Retro:** [Date]
```

---

## Status Report Template

Regular progress updates for stakeholders.

```markdown
# Project Status Report

**Project:** [Project Name]
**Reporting Period:** [Date Range]
**Prepared By:** [Name]
**Date:** [Date]

---

## Executive Summary

**Overall Status:** 🟢 On Track / 🟡 At Risk / 🔴 Off Track

[2-3 sentence summary of current state and key developments]

---

## Progress This Period

### Completed

✅ [Milestone/deliverable 1]
✅ [Milestone/deliverable 2]
✅ [Milestone/deliverable 3]

### In Progress

🔄 [Active work item 1] - [% complete]
🔄 [Active work item 2] - [% complete]

### Planned Next Period

📅 [What's coming up next]
📅 [What's coming up next]

---

## Metrics

| Metric | Target | Actual | Status |
|--------|--------|--------|--------|
| % Complete | [X%] | [Y%] | 🟢/🟡/🔴 |
| Budget Used | [X%] | [Y%] | 🟢/🟡/🔴 |
| On Schedule | On time | [Days ahead/behind] | 🟢/🟡/🔴 |
| Quality (Defects) | < [N] | [Actual] | 🟢/🟡/🔴 |

---

## Milestones

| Milestone | Target Date | Forecast Date | Status |
|-----------|-------------|---------------|--------|
| [Milestone 1] | [Date] | [Date] | ✅ Complete |
| [Milestone 2] | [Date] | [Date] | 🔄 In Progress |
| [Milestone 3] | [Date] | [Date] | 📅 Upcoming |
| [Milestone 4] | [Date] | [Date] | ⏸️ Not Started |

---

## Risks & Issues

### New This Period

🆕 **[Risk/Issue]** - Priority: High/Medium/Low
- Impact: [Description]
- Mitigation: [Plan]
- Owner: [Name]

### Ongoing

🔴 **[High priority issue]**
- Status: [Update]
- Action: [Next steps]

🟡 **[Medium priority issue]**
- Status: [Update]
- Action: [Next steps]

### Resolved

✅ **[Resolved issue]** - [How it was resolved]

---

## Budget

**Total Budget:** $[Amount]
**Spent to Date:** $[Amount] ([X%])
**Forecast to Complete:** $[Amount]
**Variance:** $[Amount] (Over/Under)

---

## Resource Status

| Resource | Allocation | Availability | Issues |
|----------|-----------|--------------|--------|
| [Name] | 100% | Available | None |
| [Name] | 50% | PTO [dates] | Coverage plan in place |

---

## Decisions Needed

1. **[Decision required]**
   - Context: [Why this decision is needed]
   - Options: [A, B, C]
   - Recommendation: [Your recommendation]
   - Needed by: [Date]

---

## Help Needed

- [Blocker requiring stakeholder intervention]
- [Resource or approval needed]

---

## Key Accomplishments

Highlights worth celebrating:

- [Achievement 1]
- [Achievement 2]

---

**Next Report:** [Date]
**Questions?** Contact [Name] at [email]
```

---

## Milestone Tracker

Track key project milestones.

```markdown
# Milestone Tracker

**Project:** [Project Name]
**Updated:** [Date]

| # | Milestone | Description | Target Date | Forecast Date | Actual Date | Status | % Complete | Owner | Dependencies | Notes |
|---|-----------|-------------|-------------|---------------|-------------|--------|------------|-------|--------------|-------|
| M1 | [Name] | [Brief desc] | 2024-03-15 | 2024-03-15 | 2024-03-14 | ✅ Complete | 100% | [Name] | None | Completed early |
| M2 | [Name] | [Brief desc] | 2024-04-01 | 2024-04-05 | - | 🟡 At Risk | 75% | [Name] | M1 | Delayed by [reason] |
| M3 | [Name] | [Brief desc] | 2024-04-15 | 2024-04-15 | - | 🟢 On Track | 25% | [Name] | M1, M2 | Started early |
| M4 | [Name] | [Brief desc] | 2024-05-01 | TBD | - | ⏸️ Not Started | 0% | [Name] | M3 | Awaiting M3 completion |

**Status Legend:**
- ✅ Complete
- 🟢 On Track
- 🟡 At Risk
- 🔴 Off Track
- ⏸️ Not Started
- 🚫 Blocked
```

---

## Gantt Chart Structure

Text-based timeline (for ASCII or markdown rendering).

```
Project Timeline (Weeks)
                    1   2   3   4   5   6   7   8   9   10  11  12
Phase 1: Planning   ████████
  - Requirements    ████
  - Design              ████

Phase 2: Build              ████████████████
  - Frontend Dev            ████████
  - Backend Dev                 ████████
  - Integration                         ████

Phase 3: Test                               ████████
  - QA Testing                              ████
  - UAT                                         ████

Phase 4: Deploy                                     ████
  - Deployment                                      ██
  - Training                                          ██

Milestones:
↓ M1: Requirements (Week 2)
        ↓ M2: Design Complete (Week 4)
                    ↓ M3: Code Complete (Week 8)
                                ↓ M4: Testing Done (Week 10)
                                        ↓ M5: Launch (Week 11)
```

### Dependency Diagram

```
[Requirements]
      ↓
[Design]
      ↓
[Development] → [Testing] → [Deployment]
      ↓              ↓
[Documentation]  [Training]
```

---

## Using These Templates

**Customization:**
- Copy the template that fits your need
- Replace bracketed placeholders with your content
- Add/remove sections as appropriate
- Adapt format (markdown, spreadsheet, docs) to your tools

**Best Practices:**
- Keep templates concise - don't add fields you won't use
- Review and update regularly - outdated plans are useless
- Share with stakeholders - plans are communication tools
- Version control - track changes over time
- Archive completed plans - reference for future projects

**Template Selection Guide:**

| Need | Use Template |
|------|-------------|
| Project kickoff | Project Charter |
| Task decomposition | WBS |
| Track risks | Risk Register |
| Define roles | RACI Matrix |
| Stakeholder summary | Project Plan (One-Pager) |
| Agile iteration | Sprint Planning |
| Continuous improvement | Retrospective |
| Regular updates | Status Report |
| Timeline overview | Gantt Chart / Milestone Tracker |

---
