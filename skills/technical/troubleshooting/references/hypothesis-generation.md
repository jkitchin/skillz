# Hypothesis Generation

Guide to forming, ranking, and testing hypotheses during troubleshooting.

## Forming Good Hypotheses

**Testable hypothesis structure:**
"If [cause] is true, then [observable consequence] should occur"

**Example:**
- ❌ Weak: "The database is slow"
- ✅ Strong: "If the database query is missing an index on the user_id column, then queries filtering by user_id should take >1 second while queries on indexed columns take <100ms"

## Hypothesis Ranking

**Rank by likelihood using:**
1. **Recent changes** - Most likely cause
2. **Common failures** - Frequent issues in this domain
3. **Evidence strength** - Symptoms match this cause
4. **Simplicity** - Simpler explanations more likely (Occam's Razor)

**Example ranking:**
1. HIGH: Config typo in recent deploy (recent change + common)
2. MEDIUM: Dependency version mismatch (somewhat common)
3. LOW: Hardware failure (rare, no evidence)

## Test Design

**Good test characteristics:**
- **Specific** - Tests one hypothesis, not multiple
- **Observable** - Clear pass/fail outcome
- **Quick** - Can execute rapidly
- **Safe** - Won't cause additional damage
- **Reversible** - Can undo if needed

## Common Failure Patterns

See domain-specific-patterns.md for comprehensive lists.

**Quick reference:**
- **Null/undefined**: Missing validation
- **Type errors**: Data format mismatch
- **Permissions**: Access control issues
- **Connection failures**: Network/service unavailable
- **Timeouts**: Performance bottlenecks
- **Resource exhaustion**: Memory/disk/connections depleted

## Avoiding False Leads

**Red herrings:**
- Correlation without causation
- Symptoms confused with root cause
- Distractions from obvious answer
- Confirmation bias (seeing what you expect)

**Stay objective:**
- Test hypotheses that could prove you wrong
- Don't commit to theory before evidence
- Follow data, not intuition alone
