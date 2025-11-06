# Domain-Specific Troubleshooting Patterns

Common failure modes and diagnostic approaches by domain.

## Software/Code Debugging

### Syntax Errors
**Symptoms:** Code won't compile/run
**Diagnostic:** Compiler/interpreter error message points to line
**Fix:** Correct syntax at indicated location

### Null/Undefined References
**Symptoms:** NullPointerException, TypeError, AttributeError
**Diagnostic:** Stack trace shows where null was accessed
**Fix:** Add null checks, ensure object initialization

### Logic Errors
**Symptoms:** Wrong output, incorrect calculations
**Diagnostic:** Debug with breakpoints, trace execution
**Fix:** Review algorithm, check conditionals

### Race Conditions
**Symptoms:** Intermittent failures, non-deterministic behavior
**Diagnostic:** Happens under load, timing-dependent
**Fix:** Proper synchronization, locks, atomic operations

### Memory Leaks
**Symptoms:** Gradual slowdown, eventual crash
**Diagnostic:** Memory usage grows over time
**Fix:** Profile memory, find unreleased resources

## System/Infrastructure

### Network Connectivity
**Symptoms:** Connection refused, timeouts
**Diagnostic:** ping, telnet, curl to test connectivity
**Fix:** Check firewall, routing, service status

### Permission Denied
**Symptoms:** Access denied, 403 errors
**Diagnostic:** Check file permissions, user roles
**Fix:** Grant appropriate permissions

### Port Conflicts
**Symptoms:** "Address already in use"
**Diagnostic:** `netstat -tulpn | grep <port>`
**Fix:** Stop conflicting service or use different port

### Resource Exhaustion
**Symptoms:** Out of memory, disk full, too many connections
**Diagnostic:** Monitor resources (top, df, netstat)
**Fix:** Increase limits, fix leaks, clean up

### Service Not Running
**Symptoms:** Connection refused to known service
**Diagnostic:** `systemctl status <service>` or `ps aux | grep <service>`
**Fix:** Start service, check why it stopped

## Configuration

### Typos
**Symptoms:** Setting not taking effect, errors about unknown keys
**Diagnostic:** Compare against documentation, check for misspellings
**Fix:** Correct spelling

### Wrong Environment
**Symptoms:** Works in dev, fails in prod
**Diagnostic:** Compare config files between environments
**Fix:** Use correct config for environment

### Path Issues
**Symptoms:** File not found despite file existing
**Diagnostic:** Absolute vs relative paths, working directory
**Fix:** Use absolute paths or verify working directory

### Type Mismatches
**Symptoms:** Config parser errors
**Diagnostic:** Config expects number but got string
**Fix:** Correct value type in config

## Data Issues

### Invalid Format
**Symptoms:** Parse errors, validation failures
**Diagnostic:** Examine actual data format vs expected
**Fix:** Convert data to correct format

### Null/Empty Values
**Symptoms:** Errors when processing data
**Diagnostic:** Missing required fields
**Fix:** Validate inputs, handle nulls gracefully

### Character Encoding
**Symptoms:** Garbled text, special characters wrong
**Diagnostic:** Check encoding (UTF-8, Latin-1, etc.)
**Fix:** Convert to correct encoding

### Schema Mismatch
**Symptoms:** Data doesn't match expectations
**Diagnostic:** Compare data structure to schema
**Fix:** Migrate data or update schema

## Integration Issues

### API Changes
**Symptoms:** Requests fail unexpectedly
**Diagnostic:** API endpoint or format changed
**Fix:** Update to new API version

### Authentication Failures
**Symptoms:** 401 Unauthorized
**Diagnostic:** Check credentials, tokens, expiration
**Fix:** Refresh credentials, update config

### Rate Limiting
**Symptoms:** 429 Too Many Requests
**Diagnostic:** Exceeded API rate limits
**Fix:** Implement backoff, reduce request rate

### Timeouts
**Symptoms:** Requests time out
**Diagnostic:** Network slow, service overloaded
**Fix:** Increase timeout, optimize requests, add retry logic

## Hardware Issues

### Power Problems
**Symptoms:** Unexpected shutdowns, reboots
**Diagnostic:** Check power supply, cables
**Fix:** Replace faulty power components

### Connection Failures
**Symptoms:** Device not detected
**Diagnostic:** Check physical connections
**Fix:** Reseat cables, replace if damaged

### Overheating
**Symptoms:** Performance throttling, shutdowns
**Diagnostic:** Check temperatures, fans
**Fix:** Improve cooling, clean dust

### Disk Failures
**Symptoms:** I/O errors, file corruption
**Diagnostic:** SMART status, disk diagnostics
**Fix:** Replace failing disk, restore from backup

## Process/Workflow Issues

### Missing Steps
**Symptoms:** Process incomplete, errors downstream
**Diagnostic:** Map actual vs expected process flow
**Fix:** Document and enforce all required steps

### Communication Gaps
**Symptoms:** Misunderstandings, delays
**Diagnostic:** Identify handoff points where info lost
**Fix:** Clarify communication protocols

### Unmet Dependencies
**Symptoms:** Can't proceed, waiting on prerequisites
**Diagnostic:** Check dependency status
**Fix:** Ensure dependencies complete before starting

### Timing Problems
**Symptoms:** Steps out of sequence
**Diagnostic:** Review order of operations
**Fix:** Enforce proper sequencing

## Quick Diagnostic Reference

| Symptom | Likely Cause | First Check |
|---------|--------------|-------------|
| Error message | See exact text | Logs, stack trace |
| Slow performance | Resource issue | CPU, memory, disk I/O |
| Intermittent failure | Race condition or load | Reproduce under load |
| Works locally, not prod | Environment difference | Compare configs, versions |
| Sudden failure | Recent change | What was deployed/changed |
| Gradual degradation | Resource leak | Memory/disk usage over time |
| Connection refused | Service down or firewall | Check service status, network |
| Permission denied | Access control | Check permissions, auth |
| Data corruption | Encoding or schema | Validate data format |
| Timeout | Network or performance | Check latency, load |
