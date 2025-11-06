# Information Gathering

Comprehensive guide to collecting diagnostic information systematically.

## Table of Contents

- [Essential Diagnostic Questions](#essential-diagnostic-questions)
- [Log Analysis Techniques](#log-analysis-techniques)
- [Error Message Interpretation](#error-message-interpretation)
- [Environment Documentation](#environment-documentation)
- [Reproduction Case Development](#reproduction-case-development)
- [Version and Dependency Checking](#version-and-dependency-checking)
- [Configuration Review](#configuration-review)

---

## Essential Diagnostic Questions

Ask these questions systematically to gather complete picture.

### Problem Description

**What's happening:**
- What are you trying to do? (intended action)
- What do you expect to happen? (expected behavior)
- What actually happens? (actual behavior)
- What error messages do you see? (exact text, screenshots)
- Can you show me? (demonstration, screen recording)

**Specifics:**
- Is there an error code or message?
- What's the exact wording of the error?
- Where does it occur? (which screen, function, step)
- What are the symptoms? (slow, crash, wrong output)

### Timeline Questions

**When:**
- When did you first notice this?
- Did it ever work correctly?
- If it worked before, when did it stop working?
- When was the last time it worked?
- Does it happen every time or sometimes?

**What Changed:**
- What changed since it last worked?
- Were there any updates? (code, dependencies, OS, config)
- Was anything deployed recently?
- Were there configuration changes?
- Did data volume/patterns change?
- Did user behavior change?

### Scope Questions

**Who:**
- Does this affect all users or specific users?
- Which users are affected? (role, location, device)
- Can you reproduce it with test account?
- Does it affect administrators too?

**Where:**
- Which environments are affected? (dev, staging, production)
- Does it work in some environments but not others?
- Which servers/regions/data centers?
- Which browsers/devices/platforms?

**When (timing):**
- Does it happen all the time or intermittently?
- Is there a pattern? (time of day, day of week)
- Does it correlate with traffic/load?
- Does it happen immediately or after some delay?

### Reproduction Questions

**Steps:**
- Can you reproduce it reliably?
- What are the exact steps to reproduce?
- Can you walk me through what you did?
- What data/input did you use?
- What was the starting state?

**Minimal Case:**
- What's the simplest scenario that fails?
- Can you remove any steps and still reproduce?
- What's the minimum data needed?
- Does it work with different data?

### Environment Questions

**System:**
- What OS and version?
- What browser and version?
- What application version?
- What database/runtime/framework versions?
- Are you behind a proxy/firewall/VPN?

**Configuration:**
- What settings are you using?
- Any non-default configuration?
- What environment variables are set?
- Which feature flags are enabled?

### Impact Questions

**Severity:**
- How many users are affected?
- What's the business impact?
- Is there a workaround?
- Is this blocking critical work?
- Is data at risk?

**Urgency:**
- How urgent is this?
- When do you need it fixed by?
- Is this impacting production?
- Are customers complaining?

---

## Log Analysis Techniques

### Where to Look

**Application Logs:**
- Application server logs
- Web server logs (nginx, Apache)
- Runtime/interpreter logs
- Framework logs

**System Logs:**
- Syslog/system event logs
- Kernel logs (dmesg)
- Authentication logs
- Cron/scheduled task logs

**Service Logs:**
- Database logs
- Message queue logs
- Cache server logs
- Load balancer logs

**Infrastructure Logs:**
- Container/orchestration logs (Docker, Kubernetes)
- Cloud provider logs (AWS CloudWatch, Azure Monitor)
- Network equipment logs
- Monitoring/APM logs

### What to Look For

**Error Indicators:**
- ERROR, FATAL, CRITICAL level messages
- Exception/stack trace entries
- WARN messages (may precede errors)
- 4xx, 5xx HTTP status codes
- Timeout messages
- Connection refused/failed
- Out of memory errors
- Permission denied

**Patterns:**
- Repeated errors (same message many times)
- Error sequences (A fails, then B, then C)
- Temporal patterns (errors cluster at certain times)
- Correlation (errors across multiple services at same time)

**Context Around Errors:**
- What was system doing before error?
- What was the last successful operation?
- What requests were being processed?
- What was resource usage?

### Log Analysis Process

**1. Identify failure time:**
```
When did problem occur?
Look for timestamp in logs around that time
```

**2. Search for errors:**
```bash
# Find errors in logs
grep -i "error" application.log
grep -E "ERROR|FATAL|EXCEPTION" application.log

# Find errors around specific time
grep "2024-01-15 14:3" application.log | grep -i error

# Count error types
grep -i error application.log | sort | uniq -c | sort -rn
```

**3. Get context:**
```bash
# Show 10 lines before and after error
grep -B 10 -A 10 "NullPointerException" application.log

# Follow log trail for specific request
grep "request_id=abc123" application.log
```

**4. Look for patterns:**
```bash
# Find all unique error messages
grep ERROR application.log | sed 's/[0-9]//g' | sort | uniq -c

# Errors by hour
grep ERROR application.log | awk '{print $1" "$2}' | cut -d: -f1 | uniq -c
```

### Common Log Patterns

**Stack Traces:**
```
ERROR: NullPointerException: Cannot invoke method on null object
    at com.example.UserService.getUser(UserService.java:47)
    at com.example.Controller.handleRequest(Controller.java:123)
    at javax.servlet.http.HttpServlet.service(HttpServlet.java:790)
```
**Read bottom-to-top:** Deepest call is at top, your code is usually near top

**Chained Exceptions:**
```
ERROR: RuntimeException: Failed to process request
    Caused by: SQLException: Connection timeout
    Caused by: SocketTimeoutException: Read timed out
```
**Root cause is deepest:** `SocketTimeoutException` is the real issue

**Request Tracing:**
```
INFO [req-abc123] Processing user login
DEBUG [req-abc123] Validating credentials
DEBUG [req-abc123] Querying database
ERROR [req-abc123] Database timeout after 30s
```
**Follow request ID** through all log entries to see flow

### Log Analysis Tools

**Command-line:**
- `grep` - Search for patterns
- `awk`/`sed` - Text processing
- `tail -f` - Follow log in real-time
- `less` - Browse large files
- `jq` - Parse JSON logs

**Log Management:**
- ELK Stack (Elasticsearch, Logstash, Kibana)
- Splunk
- Datadog
- CloudWatch Logs Insights
- Papertrail

---

## Error Message Interpretation

### Error Message Anatomy

```
ERROR: [ErrorType] ErrorMessage at LocationInfo
       ContextInfo
       StackTrace
```

**Example:**
```
ERROR: NullPointerException: Cannot read property 'email' of null
       at UserService.sendEmail (UserService.java:line 47)
       in thread "main"
```

**Components:**
- **ErrorType**: `NullPointerException` (what kind of error)
- **Message**: `Cannot read property 'email' of null` (what went wrong)
- **Location**: `UserService.java:line 47` (where it happened)
- **Context**: `in thread "main"` (under what conditions)

### Common Error Patterns

**Null/Undefined Reference:**
```
NullPointerException
TypeError: Cannot read property 'x' of null
AttributeError: 'NoneType' object has no attribute 'x'
```
**Meaning:** Trying to use an object that doesn't exist
**Fix:** Check if object exists before using

**Type Mismatch:**
```
TypeError: expected string, got int
ClassCastException: Cannot cast Integer to String
```
**Meaning:** Wrong data type provided
**Fix:** Convert to correct type or fix data source

**Not Found:**
```
FileNotFoundException: config.yaml not found
ModuleNotFoundError: No module named 'requests'
404 Not Found
```
**Meaning:** Resource doesn't exist
**Fix:** Check path, install dependency, verify URL

**Permission Denied:**
```
PermissionError: [Errno 13] Permission denied
403 Forbidden
Access denied for user
```
**Meaning:** Insufficient permissions
**Fix:** Check file permissions, user roles, API credentials

**Connection Failed:**
```
Connection refused
ECONNREFUSED
Timeout
Network unreachable
```
**Meaning:** Can't reach service/server
**Fix:** Check if service is running, network connectivity, firewall

**Resource Exhausted:**
```
OutOfMemoryError
Too many open files
Connection pool exhausted
Disk full
```
**Meaning:** Ran out of system resources
**Fix:** Increase limits, fix leaks, clean up resources

**Syntax Error:**
```
SyntaxError: unexpected token '}'
ParseError: missing semicolon
IndentationError
```
**Meaning:** Code has syntax mistake
**Fix:** Check syntax at indicated line

### Reading Stack Traces

**Stack trace order (varies by language):**

**Java/Python** (top = most recent):
```
Most recent call
    ↓
    ↓
    ↓
Original call (main/entry)
```

**JavaScript** (bottom = most recent):
```
Original call (entry)
    ↓
    ↓
    ↓
Most recent call
```

**What to look for:**
1. **Your code** - Look for your package/module names
2. **Error location** - First mention of your code
3. **Entry point** - How did we get there?
4. **Framework code** - Usually not the issue (unless framework bug)

**Example (Java):**
```
Exception in thread "main" java.lang.NullPointerException
    at com.myapp.UserService.sendEmail(UserService.java:47)    ← YOUR CODE (likely issue)
    at com.myapp.Controller.processRequest(Controller.java:123) ← YOUR CODE
    at javax.servlet.http.HttpServlet.service(HttpServlet.java:790) ← Framework
    at org.apache.catalina.core.StandardWrapper.service(StandardWrapper.java:1280) ← Framework
```

**Focus on:** `UserService.java:47` - that's where null happened in your code

---

## Environment Documentation

### Complete Environment Checklist

**Software Versions:**
```
Application version: _______
Runtime/Platform: _______ (e.g., Node 18.2.0, Python 3.11, Java 17)
OS: _______ (e.g., Ubuntu 22.04, macOS 14.1, Windows Server 2019)
Database: _______ (e.g., PostgreSQL 15.2, MySQL 8.0)
Dependencies: Check package.json, requirements.txt, pom.xml
```

**Configuration:**
```
Config files: List all config files and their locations
Environment variables: List all env vars (redact secrets)
Feature flags: Which features are enabled/disabled
Settings: Any non-default settings
```

**Infrastructure:**
```
Server/Container: Physical, VM, container?
Network: Topology, firewalls, load balancers
Storage: Disk type, mount points, capacity
Resources: CPU, RAM, disk allocated
```

**Dependencies:**
```
External services: APIs, databases, message queues
Service versions: Check versions of all external services
Network endpoints: URLs, IPs, ports
Authentication: How services authenticate
```

### Environment Comparison Matrix

Compare environments to identify differences:

| Component | Development | Staging | Production |
|-----------|-------------|---------|------------|
| App Version | v1.2.4-dev | v1.2.3 | v1.2.3 |
| Node Version | 18.2.0 | 18.2.0 | 16.14.0 ← DIFFERENCE |
| Database | PostgreSQL 15 | PostgreSQL 15 | PostgreSQL 14 ← DIFFERENCE |
| Config | dev.yaml | staging.yaml | prod.yaml |
| Data Volume | 1K records | 10K records | 1M records ← DIFFERENCE |
| Traffic | 1 req/min | 10 req/min | 1000 req/min ← DIFFERENCE |

**Investigate differences** - these are likely culprits for "works in dev, fails in prod"

---

## Reproduction Case Development

### Goals of Reproduction Case

**A good reproduction case is:**
1. **Reliable** - Reproduces every time (or predictably)
2. **Minimal** - Fewest steps, smallest data
3. **Isolated** - Removes non-essential dependencies
4. **Documented** - Clear steps anyone can follow
5. **Self-contained** - Includes all necessary data/config

### Building Minimal Reproducible Example (MRE)

**Process:**

**1. Start with full reproduction:**
```
All steps user took
All data they used
Their full environment
```

**2. Remove non-essential steps:**
```
Can you skip any steps and still reproduce?
What's the shortest path to the error?
```

**3. Minimize data:**
```
Does it work with less data?
What's the simplest input that fails?
Can you use dummy/synthetic data?
```

**4. Remove dependencies:**
```
Can you reproduce without database?
Can you mock external APIs?
Can you use local file instead of remote?
```

**5. Simplify environment:**
```
Does it need production config?
Can you reproduce with default settings?
Can you reproduce in fresh environment?
```

**6. Document and verify:**
```
Write step-by-step instructions
Test instructions on fresh setup
Ensure anyone can reproduce
```

### MRE Template

```markdown
## Minimal Reproducible Example

**Environment:**
- OS: Ubuntu 22.04
- Python: 3.11.2
- Required packages: requests==2.28.0

**Setup:**
1. Install dependencies: `pip install requests==2.28.0`
2. Create test file: `test.py` (code below)

**Code:**
[Paste minimal code that demonstrates issue]

**Data:**
[Paste minimal data needed]

**Steps to Reproduce:**
1. Run: `python test.py`
2. Observe error (actual output below)

**Expected Behavior:**
[What should happen]

**Actual Behavior:**
[What actually happens]

**Error Output:**
[Paste exact error message]
```

### Example: From Complex to Minimal

**Initial Report:**
"My web app crashes when I try to export the user report from the admin dashboard after filtering by date range and selecting multiple user types while the database backup is running."

**Minimal Reproduction:**
```python
# Reproduces crash without web app, dashboard, filters, or backup
import database

users = database.query("SELECT * FROM users WHERE type IN ('admin', 'customer')")
print(len(users))  # Crashes here

# Error: TypeError: object of type 'NoneType' has no len()
# Root cause: query returns None when multiple types specified
```

**That's an MRE:** Simple code, minimal data, reproduces issue reliably.

---

## Version and Dependency Checking

### Capturing Versions

**Application version:**
```bash
# Git commit/tag
git describe --always --tags

# Package version
npm version
python --version
java -version
```

**Dependencies:**
```bash
# Node.js
npm list --depth=0

# Python
pip freeze

# Ruby
bundle list

# Java
mvn dependency:tree
```

**System:**
```bash
# OS version
uname -a                    # Linux/Mac
lsb_release -a             # Linux
sw_vers                     # macOS
systeminfo                  # Windows

# Kernel version
uname -r
```

### Dependency Conflicts

**Common scenarios:**

**Version mismatch:**
```
App depends on libA 2.0
  libA 2.0 depends on libB 1.5

App also depends on libC 3.0
  libC 3.0 depends on libB 2.0

CONFLICT: libB needs to be both 1.5 and 2.0
```

**Check for conflicts:**
```bash
# Node.js
npm ls <package>      # Show all versions of package

# Python
pip list | grep <package>
```

**Transitive dependencies:**
```
You only specified dependency A
But A depends on B, C, D
And B depends on E, F
And C depends on E (different version!)
```

**Solution:** Lock file captures entire tree
- `package-lock.json` (Node.js)
- `Pipfile.lock` (Python)
- `Gemfile.lock` (Ruby)

### Comparing Environments

**Compare lock files between environments:**
```bash
diff production/package-lock.json staging/package-lock.json
```

Differences = potential cause of environment-specific bugs

---

## Configuration Review

### Configuration Sources

**Files:**
- Application config files (.yaml, .json, .ini, .conf)
- Framework config files
- Server config (nginx.conf, apache.conf)
- System config (/etc/)

**Environment Variables:**
```bash
# List all environment variables
printenv

# Check specific variable
echo $DATABASE_URL
echo %PATH%  # Windows
```

**Command-line Arguments:**
- How was application started?
- What flags/options were passed?

**Database/Remote Config:**
- Feature flags
- Remote configuration services
- Database-stored settings

### Config Review Checklist

**For each config setting:**
- What's the current value?
- What's the default value?
- What's the expected value?
- Was it recently changed?
- Is it different between environments?

**Common config issues:**

**Typos:**
```yaml
databse_url: postgres://...   # Typo: databse instead of database
```

**Wrong environment:**
```yaml
database_url: localhost:5432   # Production using localhost!
```

**Missing values:**
```yaml
api_key:    # Empty!
```

**Wrong type:**
```yaml
port: "8080"   # String instead of number
max_connections: 100.5   # Float instead of int
```

**Path issues:**
```yaml
log_file: logs/app.log   # Relative path, but where's working directory?
```

**Sensitive data hardcoded:**
```yaml
password: supersecret123   # In version control!
```

### Configuration Diff

Compare configs between environments:

```bash
# Show differences
diff dev.yaml prod.yaml

# Side-by-side comparison
diff -y dev.yaml prod.yaml
```

Any differences are suspects when troubleshooting environment-specific issues.

---

## Information Gathering Workflow

**When debugging, systematically gather:**

1. ✅ **Problem description** - Expected vs actual behavior
2. ✅ **Error messages** - Exact text, full stack trace
3. ✅ **Logs** - Application, system, service logs around failure time
4. ✅ **Environment** - Versions, config, infrastructure
5. ✅ **Reproduction steps** - Minimal reliable reproduction
6. ✅ **Timeline** - When it started, what changed
7. ✅ **Scope** - Who/where/when affected
8. ✅ **Impact** - How many users, business impact

**Document everything** as you gather - you'll need it when hypothesizing causes.

---

**Remember:** Time spent gathering information is never wasted. Complete information makes diagnosis faster and more accurate.
