# Diagnostic Frameworks

Comprehensive guides to structured diagnostic methodologies for root cause analysis and systematic troubleshooting.

## Table of Contents

- [5 Whys Method](#5-whys-method)
- [Fishbone/Ishikawa Diagram](#fishboneishikawa-diagram)
- [FMEA (Failure Mode and Effects Analysis)](#fmea-failure-mode-and-effects-analysis)
- [Binary Search Troubleshooting](#binary-search-troubleshooting)
- [Rubber Duck Debugging](#rubber-duck-debugging)
- [Comparative Analysis](#comparative-analysis)
- [Hypothesis-Driven Debugging](#hypothesis-driven-debugging)

---

## 5 Whys Method

### Overview

Simple but powerful technique for finding root causes by repeatedly asking "why" - typically 5 times, but may need more or fewer.

### How to Use

**Process:**
1. Start with the problem statement
2. Ask "Why did this happen?"
3. Answer with facts (not speculation)
4. Ask "Why?" about that answer
5. Repeat until you reach an actionable root cause

**Stop when:**
- You reach a cause you can fix
- You identify a process/system failure
- Further "whys" become speculative
- You've identified the true root (not just symptoms)

### Example 1: Software Bug

**Problem:** Users can't log in to the application

1. **Why?** Login API is returning 500 errors
2. **Why?** Database connection is failing
3. **Why?** Connection pool is exhausted
4. **Why?** Connections aren't being released after queries
5. **Why?** Developer forgot to call `connection.close()` in error handling paths
6. **Why?** (Optional) No code review checklist item for resource cleanup
7. **Root cause:** Missing resource cleanup in error paths + inadequate code review process

**Fixes:**
- Immediate: Add `connection.close()` in finally block
- Preventive: Add resource cleanup to code review checklist
- Systemic: Implement linting rule to detect unclosed resources

### Example 2: Process Failure

**Problem:** Customer shipments are consistently late

1. **Why?** Warehouse doesn't pack orders on time
2. **Why?** Warehouse runs out of packing materials mid-day
3. **Why?** Materials order quantity is insufficient
4. **Why?** Order quantity based on outdated sales projections
5. **Why?** Sales team doesn't share updated forecasts with warehouse
6. **Root cause:** Communication breakdown between sales and warehouse teams

**Fixes:**
- Immediate: Order more packing materials
- Preventive: Establish weekly forecast sharing meeting
- Systemic: Implement shared sales dashboard for warehouse visibility

### Guidelines

**Do:**
- Focus on facts, not blame
- Ask "why" about the system/process, not people
- Stop when you reach actionable root cause
- May need more than 5 whys for complex issues
- Document each step

**Don't:**
- Stop at symptoms ("database was slow")
- Blame individuals ("John made a mistake")
- Speculate without evidence
- Continue past actionable root cause

---

## Fishbone/Ishikawa Diagram

### Overview

Visual method to explore all possible causes of a problem, organized by category. Shaped like fish skeleton.

### Structure

```
                       ┌─── Category 1
                       │    ├─ Cause 1.1
                       │    └─ Cause 1.2
                       │
     Category 2 ───┐   │
     ├─ Cause 2.1  ├───┼──► PROBLEM
     └─ Cause 2.2  │   │
                       │
                       │    ├─ Cause 3.1
                       └─── Category 3
                            └─ Cause 3.2
```

### Common Categories

**6 M's (Manufacturing):**
- **M**an/People - Human factors
- **M**achine - Equipment, tools
- **M**ethod - Processes, procedures
- **M**aterial - Inputs, supplies
- **M**easurement - Metrics, monitoring
- **M**other Nature - Environment, conditions

**4 P's (Services):**
- **P**olicies - Rules, regulations
- **P**rocedures - Workflows, methods
- **P**eople - Skills, training
- **P**lant/Technology - Systems, tools

**Software Development:**
- **Code** - Logic, algorithms, syntax
- **Configuration** - Settings, environment
- **Dependencies** - Libraries, services
- **Data** - Inputs, state, database
- **Infrastructure** - Network, servers, resources
- **Process** - Workflow, communication

### How to Create

**Process:**

1. **Define problem clearly** - Write at the "head" of the fish

2. **Identify major categories** - Draw "bones" from spine (use 6 M's, 4 P's, or custom)

3. **Brainstorm causes** - For each category, list possible causes
   - What in this category could cause the problem?
   - Don't self-censor yet - list everything

4. **Add sub-causes** - Add "bones" off the causes
   - Why does this cause happen?
   - What are contributing factors?

5. **Analyze** - Look for:
   - Causes appearing in multiple categories (likely suspects)
   - Causes with many sub-causes (complex issues)
   - Causes you can verify quickly

6. **Prioritize** - Which causes are most likely? Most impactful?

7. **Test** - Design experiments to confirm/eliminate causes

### Example: Website Performance Problem

**Problem:** Website response time increased from 200ms to 2000ms

**Fishbone Analysis:**

**Code:**
- Inefficient database queries
  - Missing indexes
  - N+1 query problem
- Memory leaks
- Unoptimized algorithms

**Configuration:**
- Insufficient server resources
  - Too few worker processes
  - Low memory allocation
- Cache disabled
- Compression disabled

**Dependencies:**
- Database server slow
  - High query load
  - Disk I/O bottleneck
- Third-party API timeout increased
- CDN issues

**Data:**
- Dataset grew significantly
  - Tables not optimized for size
- User session data bloat

**Infrastructure:**
- Network latency
- Load balancer misconfigured
- Database connection pool exhausted

**Process:**
- No performance testing before deploy
- Monitoring alerts not configured

**Analysis:**
- "Database" appears in multiple categories → High priority
- "Recent deploy" timeline → Check what changed
- Test database queries first (likely culprit)

### Benefits

- Comprehensive (ensures you consider all angle)
- Collaborative (team can brainstorm together)
- Visual (easy to see relationships)
- Systematic (organized by category)

---

## FMEA (Failure Mode and Effects Analysis)

### Overview

Systematic method to identify potential failures, their causes, and effects. Prioritize by risk (likelihood × severity).

### FMEA Table Structure

| Failure Mode | Potential Causes | Effects | Severity (1-10) | Likelihood (1-10) | Detection (1-10) | RPN | Actions |
|--------------|-----------------|---------|-----------------|-------------------|------------------|-----|---------|
| What can fail | Why it fails | Impact if fails | How bad | How often | Can we catch it | Risk Priority Number | What to do |

**Severity (1-10):**
- 1-3: Minor (negligible impact)
- 4-6: Moderate (some impact, workarounds exist)
- 7-9: Major (significant impact, business disruption)
- 10: Catastrophic (system failure, data loss, safety risk)

**Likelihood (1-10):**
- 1-2: Remote (almost never happens)
- 3-4: Low (rarely happens)
- 5-6: Moderate (happens occasionally)
- 7-8: High (happens frequently)
- 9-10: Very High (happens almost always)

**Detection (1-10):**
- 1-2: Very high (almost certain to detect before impact)
- 3-4: High (good monitoring, likely to catch)
- 5-6: Moderate (may catch it, may not)
- 7-8: Low (difficult to detect)
- 9-10: Very low (almost impossible to detect)

**RPN (Risk Priority Number) = Severity × Likelihood × Detection**

Higher RPN = higher priority to address

### Example: Web Application FMEA

| Failure Mode | Causes | Effects | Sev | Like | Det | RPN | Actions |
|--------------|--------|---------|-----|------|-----|-----|---------|
| Database connection failure | Network issue, DB down, auth failure | Site down, users can't access | 9 | 3 | 2 | 54 | Connection pooling, retry logic, monitoring alerts |
| Memory leak | Unclosed resources, circular refs | Slow performance, eventual crash | 7 | 5 | 6 | 210 | Code review, profiling, automated leak detection |
| API rate limiting | Too many requests, no backoff | Features fail, user frustration | 6 | 7 | 4 | 168 | Implement rate limiting, request queue, exponential backoff |
| Cache invalidation bug | Race condition, stale data | Users see old data | 5 | 6 | 7 | 210 | Cache versioning, TTL tuning, cache key design review |
| Disk space exhaustion | Log growth, temp files | System crash, data loss | 9 | 4 | 3 | 108 | Log rotation, disk space monitoring, auto-cleanup |

**Priority order (by RPN):**
1. Memory leak (210) → Code review + profiling
2. Cache bug (210) → Cache design review
3. API rate limiting (168) → Implement backoff
4. Disk space (108) → Add monitoring
5. Database connection (54) → Already have good detection

### When to Use FMEA

- Proactive failure prevention
- Risk assessment for new features/systems
- Prioritizing monitoring and alerting
- Planning redundancy and failover
- Compliance requirements (safety-critical systems)

---

## Binary Search Troubleshooting

### Overview

Divide problem space in half repeatedly until issue is isolated. Efficient for large search spaces.

### When to Use

- Many possible causes
- Can test "halfway" point
- Problem space is divisible
- Each test gives clear yes/no answer

### Applications

**Version/Commit Bisect:**
- Problem: Bug introduced in last 100 commits
- Method: Test commit 50, determine which half has bug, repeat
- Tool: `git bisect`

**Configuration Narrowing:**
- Problem: One of 50 config settings causes crash
- Method: Disable half, see if crash persists, repeat

**Data Range:**
- Problem: Processing fails on subset of 1M records
- Method: Process first 500K, determine which half fails, repeat

**Time Range:**
- Problem: Service degraded sometime last week
- Method: Check midweek, determine before/after, repeat

**Code Section:**
- Problem: Function has bug but it's 500 lines
- Method: Comment out half, see if bug persists, repeat

### Process

1. **Define boundaries** - Working state vs broken state

2. **Find midpoint** - Halfway between boundaries

3. **Test midpoint** - Does problem occur at midpoint?

4. **Update boundaries**
   - If problem occurs: Midpoint becomes new "broken" boundary
   - If problem doesn't occur: Midpoint becomes new "working" boundary

5. **Repeat** - Until boundaries converge on single item

6. **Verify** - Confirm identified item is the cause

### Example: Git Bisect

```bash
# Scenario: Tests passed on commit A (100 commits ago), fail on commit B (current)

git bisect start
git bisect bad HEAD          # Current commit is bad
git bisect good A            # Commit A was good

# Git checks out commit 50 (midpoint)
run_tests                    # Tests fail

git bisect bad               # Mark as bad, Git moves to commit 25

run_tests                    # Tests pass

git bisect good              # Mark as good, Git moves to commit 37

# Continue until Git identifies the exact commit that introduced the bug

git bisect reset             # Return to original state
```

### Example: Configuration Narrowing

```python
# Problem: One of these 64 settings causes crash

settings = {
    'setting_1': 'value1',
    'setting_2': 'value2',
    # ... 62 more settings
    'setting_64': 'value64'
}

# Round 1: Test with first 32 settings
# Result: Still crashes → Problem is in first 32

# Round 2: Test with first 16 settings
# Result: Doesn't crash → Problem is in settings 17-32

# Round 3: Test with settings 17-24
# Result: Still crashes → Problem is in 17-24

# Round 4: Test with settings 17-20
# Result: Still crashes → Problem is in 17-20

# Round 5: Test with settings 17-18
# Result: Doesn't crash → Problem is setting 19 or 20

# Round 6: Test with setting 19 only
# Result: Crashes → FOUND IT: setting_19 causes crash

# Only 6 tests to find 1 issue among 64 settings!
```

### Advantages

- Efficient: O(log n) instead of O(n)
- Systematic: No guessing
- Exhaustive: Guaranteed to find issue
- 64 items → 6 tests
- 1024 items → 10 tests
- 1,000,000 items → 20 tests

### Requirements

- Divisible problem space
- Clear test that gives yes/no answer
- Monotonic property (problem exists from some point onward)

Use `scripts/binary_search_helper.py` for interactive guidance.

---

## Rubber Duck Debugging

### Overview

Explain the problem out loud (to a rubber duck, colleague, or even yourself) to gain new perspective.

### Why It Works

- **Forces clarity** - Must articulate problem precisely
- **Reveals assumptions** - Speaking exposes unexamined beliefs
- **Fresh perspective** - Hearing yourself triggers new connections
- **Step-by-step** - Explaining forces sequential thinking
- **No judgment** - Duck never makes you feel dumb

### How to Use

1. **Get a duck** - Actual rubber duck, stuffed animal, houseplant, or willing colleague

2. **Explain the problem** - Start from the beginning:
   - "What I'm trying to do is..."
   - "The code should..."
   - "But what actually happens is..."

3. **Walk through code/process** - Line by line, step by step:
   - "First, we..."
   - "Then we..."
   - "At this point, X should equal..."

4. **Speak assumptions** - Say what you're thinking:
   - "I assume this function returns..."
   - "This should always be true because..."
   - "Users would never..."

5. **Listen to yourself** - Often mid-explanation:
   - "Wait, actually..."
   - "Hmm, that doesn't make sense..."
   - "Oh! I never checked if..."

### Example

```
Developer: "Okay duck, my function should calculate the average of a list of numbers.
           Let me walk through the code.

           First, I initialize sum to 0.
           Then I loop through each number and add it to sum.
           Then I return sum divided by the length of the list.

           Wait... what if the list is empty? I'd be dividing by zero!
           That's the bug! I need to check if the list is empty first."
```

**Duck didn't say a word, but problem solved.**

### Benefits

- Free (except cost of duck)
- Always available
- Non-judgmental
- Surprisingly effective

### Variants

- **Rubber duck reviewing** - Explain proposed solution
- **Rubber duck designing** - Explain architecture plan
- **Rubber duck testing** - Explain test coverage
- **Rubber duck documenting** - Explain code to write docs

---

## Comparative Analysis

### Overview

Compare working vs non-working cases to identify differences. Differences are candidates for root cause.

### Method

**Process:**

1. **Identify comparison cases**
   - Working case (baseline)
   - Non-working case (problem)

2. **Document everything about each**
   - Environment
   - Configuration
   - Data/inputs
   - Dependencies/versions
   - Timing/sequence
   - System state

3. **Compare systematically**
   - What's different?
   - What's the same?
   - What's missing in one but not the other?

4. **Test differences**
   - Make working case more like non-working
   - Make non-working case more like working
   - See which change introduces/fixes problem

### Example Comparisons

**Works on Developer Machine, Fails in Production:**

| Aspect | Developer | Production | Potential Issue |
|--------|-----------|------------|----------------|
| OS | macOS | Linux | Case-sensitive paths |
| Environment | development | production | Different config |
| Database | SQLite | PostgreSQL | SQL dialect differences |
| Data volume | 100 records | 1M records | Performance/scalability |
| Network | localhost | multi-server | Network latency, firewalls |

**Test:** Run with production config on dev machine → Isolate config vs environment

**Works for User A, Fails for User B:**

| Aspect | User A | User B | Potential Issue |
|--------|--------|--------|----------------|
| Browser | Chrome 120 | Firefox 100 | Browser compatibility |
| Permissions | Admin | Standard | Permission issue |
| Data | Profile complete | Profile incomplete | Null handling |
| Location | US | EU | Timezone, locale |
| Session | Fresh login | Long session | Session expiry, stale data |

**Test:** User A tries with Firefox → Isolate browser issue

**Worked Yesterday, Fails Today:**

| Aspect | Yesterday | Today | Potential Issue |
|--------|-----------|-------|----------------|
| Code version | v1.2.3 | v1.2.4 | Recent deploy introduced bug |
| Dependencies | lib 2.0.0 | lib 2.1.0 | Dependency update broke something |
| Data | 10K records | 15K records | Crossed threshold |
| Traffic | 100 req/s | 200 req/s | Load-related |
| Certificate | Valid | Expired | SSL cert expired |

**Test:** Roll back to v1.2.3 → Isolate code changes

### Controlled Comparison

**Make working case progressively more like failing case:**

Start with working case, change one thing at a time toward failing case:

1. Working: Dev machine, dev config, small data
2. Test: Dev machine, dev config, **production data** → Still works
3. Test: Dev machine, **production config**, production data → Still works
4. Test: **Production machine**, production config, production data → **FAILS!**

**Conclusion:** Production machine environment is the issue (investigate OS, network, installed software)

### Benefits

- Isolates variables systematically
- Works when root cause is unknown
- Reveals environmental issues
- Provides baseline for testing fixes

---

## Hypothesis-Driven Debugging

### Overview

Generate testable hypotheses about root cause, design experiments to test each, systematically eliminate until cause found.

### Scientific Method for Debugging

**Process:**

1. **Observe** - Gather facts about the problem

2. **Hypothesize** - Generate possible explanations

3. **Predict** - "If hypothesis X is true, then Y should happen"

4. **Experiment** - Design test that could prove hypothesis wrong

5. **Analyze** - Did results match prediction?
   - Yes → Hypothesis supported (but not proven)
   - No → Hypothesis falsified (ruled out)

6. **Iterate** - Refine hypothesis or test next one

7. **Conclude** - When hypothesis consistently supported and fix works

### Hypothesis Quality

**Good Hypothesis:**
- Specific (not vague)
- Testable (can design experiment)
- Falsifiable (could be proven wrong)
- Based on evidence (not wild guess)

**Examples:**

❌ **Bad:** "Something is wrong with the database"
✅ **Good:** "Database connection pool is exhausted, causing new connections to time out"

❌ **Bad:** "The code has a bug"
✅ **Good:** "The null check is missing on line 47, causing NullPointerException when user.email is null"

❌ **Bad:** "Network is slow"
✅ **Good:** "DNS resolution is adding 2-second latency to each API call"

### Example: Debugging Login Failure

**Observation:** Users can't log in, getting "Invalid credentials" error

**Hypothesis 1:** Database is down
- **Prediction:** All database queries should fail
- **Test:** Query database directly → SUCCESS
- **Result:** Database is working fine → **Hypothesis falsified**

**Hypothesis 2:** Password hashing algorithm changed
- **Prediction:** Stored hashes won't match computed hashes
- **Test:** Compare hash format in database vs code → MATCH
- **Result:** Hashing algorithm unchanged → **Hypothesis falsified**

**Hypothesis 3:** User table was accidentally truncated
- **Prediction:** User table should be empty
- **Test:** `SELECT COUNT(*) FROM users` → Returns 10,000
- **Result:** Users exist in database → **Hypothesis falsified**

**Hypothesis 4:** Environment variable for database host is wrong
- **Prediction:** App is connecting to wrong database
- **Test:** Check `DB_HOST` env var → Points to old database!
- **Verify:** Check users in old database → Table WAS truncated
- **Result:** **ROOT CAUSE FOUND!** - Environment variable points to old database that was decommissioned

**Fix:** Update `DB_HOST` to correct database

**Verify:** Users can now log in → **Success!**

### Benefits

- Systematic (not random)
- Efficient (eliminates possibilities)
- Scientific (objective, evidence-based)
- Documentable (clear trail of reasoning)

Use `scripts/hypothesis_tracker.py` to track hypothesis testing process.

---

## Framework Selection Guide

Choose the right framework for your situation:

| Framework | Best For | When to Use |
|-----------|----------|-------------|
| **5 Whys** | Finding root cause | Simple problems, clear cause-effect chains |
| **Fishbone** | Complex problems | Many possible causes, need comprehensive view |
| **FMEA** | Risk assessment | Proactive failure prevention, prioritization |
| **Binary Search** | Large search space | Many candidates, can test midpoint |
| **Rubber Duck** | Stuck/confused | Can't see the problem clearly |
| **Comparative** | Environment issues | Works somewhere but not everywhere |
| **Hypothesis-Driven** | Scientific approach | Multiple theories, need systematic testing |

**Combine frameworks:** Use Fishbone to generate hypotheses, then hypothesis-driven testing to verify. Use 5 Whys after finding issue to dig to root cause.

---

**Remember:** All frameworks are tools. Use the one that fits the problem. Be systematic, document your process, and follow the evidence.
