# Emacs Lisp TDD Guide

Comprehensive guide to Test-Driven Development in Emacs Lisp using ERT.

## ERT (Emacs Lisp Regression Testing)

Built-in testing framework for Emacs Lisp (Emacs 24+).

### Basic Test Structure

```elisp
(require 'ert)

(ert-deftest test-addition ()
  "Test that addition works correctly."
  (should (= (+ 2 3) 5))
  (should (= (+ 0 0) 0))
  (should (= (+ -1 1) 0)))
```

### Running Tests

**Interactive:**
```elisp
M-x ert RET t RET          ; Run all tests
M-x ert RET test-name RET  ; Run specific test
M-x ert RET "prefix-*" RET ; Run tests matching pattern
```

**Batch Mode (CI):**
```bash
emacs -batch -l ert -l my-package.el -l test-my-package.el -f ert-run-tests-batch-and-exit
```

## ERT Assertions

### should
Basic assertion - expression should be true:

```elisp
(should (= 5 (+ 2 3)))
(should (string= "hello" (upcase "HELLO")))
(should (member 'x '(x y z)))
```

### should-not
Expression should be false/nil:

```elisp
(should-not (= 5 6))
(should-not (string= "hello" "world"))
(should-not (null '(1 2 3)))
```

### should-error
Should signal an error:

```elisp
(should-error (/ 1 0))
(should-error (error "test error"))

;; Match specific error type
(should-error (/ 1 0) :type 'arith-error)

;; Match error message
(should-error (error "invalid input") :type 'error)
```

## Test Organization

### File Structure

```
my-package/
├── my-package.el          ; Package code
└── test/
    └── test-my-package.el ; Tests
```

### Test File Template

```elisp
;;; test-my-package.el --- Tests for my-package

(require 'ert)
(require 'my-package)

(ert-deftest test-my-function ()
  "Test my-function with various inputs."
  (should (= (my-function 5) 10))
  (should (= (my-function 0) 0)))

(provide 'test-my-package)
;;; test-my-package.el ends here
```

## Testing Patterns

### Testing Simple Functions

```elisp
(ert-deftest test-double ()
  "Test doubling a number."
  (should (= (my-double 5) 10))
  (should (= (my-double 0) 0))
  (should (= (my-double -3) -6)))
```

### Testing Buffer Manipulation

```elisp
(ert-deftest test-insert-greeting ()
  "Test inserting greeting in buffer."
  (with-temp-buffer
    (my-insert-greeting "Alice")
    (should (string= (buffer-string) "Hello, Alice!"))))

(ert-deftest test-buffer-uppercase ()
  "Test uppercasing buffer content."
  (with-temp-buffer
    (insert "hello world")
    (my-uppercase-buffer)
    (should (string= (buffer-string) "HELLO WORLD"))))
```

### Testing Interactive Functions

```elisp
(ert-deftest test-interactive-command ()
  "Test interactive command."
  (with-temp-buffer
    (insert "test")
    (goto-char (point-min))
    (call-interactively 'my-command)
    (should (string= (buffer-string) "TEST"))))
```

### Mocking User Input

```elisp
(ert-deftest test-read-input ()
  "Test function that reads user input."
  (cl-letf (((symbol-function 'read-string)
             (lambda (prompt) "mocked-input")))
    (should (string= (my-function-that-reads) "mocked-input"))))

(ert-deftest test-y-or-n-p ()
  "Test function that asks yes/no question."
  (cl-letf (((symbol-function 'y-or-n-p)
             (lambda (prompt) t)))
    (should (my-function-that-asks))))
```

### Testing with let-binding

```elisp
(ert-deftest test-with-custom-variable ()
  "Test function with custom variable value."
  (let ((my-package-setting 'custom-value))
    (should (eq (my-get-setting) 'custom-value))))
```

### Testing File Operations

```elisp
(ert-deftest test-file-read ()
  "Test reading file."
  (let ((temp-file (make-temp-file "test")))
    (unwind-protect
        (progn
          (with-temp-file temp-file
            (insert "test content"))
          (should (string= (my-read-file temp-file) "test content")))
      (delete-file temp-file))))
```

### Testing Hooks

```elisp
(ert-deftest test-hook-function ()
  "Test that hook function is called."
  (let ((called nil))
    (add-hook 'my-hook (lambda () (setq called t)))
    (run-hooks 'my-hook)
    (should called)))
```

### Testing Advices

```elisp
(ert-deftest test-advice ()
  "Test function advice."
  (let ((called-with nil))
    (advice-add 'my-function :before
                (lambda (&rest args) (setq called-with args)))
    (unwind-protect
        (progn
          (my-function 1 2 3)
          (should (equal called-with '(1 2 3))))
      (advice-remove 'my-function
                     (lambda (&rest args) (setq called-with args))))))
```

## TDD Workflow Example

### RED Phase - Write Failing Test

```elisp
;; test/test-string-utils.el
(ert-deftest test-reverse-string ()
  "Test reversing a string."
  (should (string= (reverse-string "hello") "olleh")))
```

Run: `M-x ert RET test-reverse-string RET`
Result: ❌ `void-function reverse-string`

### GREEN Phase - Minimal Implementation

```elisp
;; string-utils.el
(defun reverse-string (s)
  "Reverse string S."
  (concat (reverse (string-to-list s))))
```

Run: `M-x ert RET test-reverse-string RET`
Result: ✅ Test passes

### REFACTOR Phase

```elisp
;; string-utils.el
(defun reverse-string (s)
  "Reverse string S.

S must be a string. Returns reversed string."
  (if (not (stringp s))
      (error "Argument must be a string")
    (concat (reverse (string-to-list s)))))
```

Add test for error case:
```elisp
(ert-deftest test-reverse-string-error ()
  "Test reverse-string with non-string."
  (should-error (reverse-string 123) :type 'error))
```

Run all tests: ✅ All pass

## Test Fixtures (Setup/Teardown)

### Using unwind-protect

```elisp
(ert-deftest test-with-setup-teardown ()
  "Test with setup and teardown."
  (let ((original-value my-var))
    (unwind-protect
        (progn
          ;; Setup
          (setq my-var 'test-value)
          ;; Test
          (should (eq my-var 'test-value)))
      ;; Teardown (always runs)
      (setq my-var original-value))))
```

### Fixture Macro

```elisp
(defmacro with-test-env (&rest body)
  "Execute BODY in test environment."
  `(let ((original-value my-var))
     (unwind-protect
         (progn
           (setq my-var 'test-value)
           ,@body)
       (setq my-var original-value))))

(ert-deftest test-using-fixture ()
  "Test using custom fixture."
  (with-test-env
    (should (eq my-var 'test-value))))
```

## ERT Best Practices

### Test Naming

```elisp
;; Good: Descriptive
(ert-deftest test-add-positive-numbers ()
  "Test adding two positive numbers.")

(ert-deftest test-divide-by-zero-signals-error ()
  "Test that dividing by zero signals an error.")

;; Bad: Vague
(ert-deftest test1 ()
  "Test.")
```

### Docstrings

```elisp
;; Good: Clear description
(ert-deftest test-parse-date ()
  "Test parsing date string in ISO format.
Should handle YYYY-MM-DD format and return proper date object."
  ...)

;; Bad: No docstring or useless docstring
(ert-deftest test-parse-date ()
  ...)
```

### Test Independence

```elisp
;; Good: Independent tests
(ert-deftest test-1 ()
  (with-temp-buffer
    (my-function)
    (should ...)))

(ert-deftest test-2 ()
  (with-temp-buffer  ; Fresh buffer
    (my-other-function)
    (should ...)))

;; Bad: Tests share state
(defvar test-buffer nil)

(ert-deftest test-1 ()
  (setq test-buffer (get-buffer-create "*test*"))
  ...)

(ert-deftest test-2 ()
  (with-current-buffer test-buffer  ; Depends on test-1!
    ...))
```

## Buttercup (BDD Alternative)

Behavior-driven development style testing for Elisp.

### Installation

```bash
# Install from MELPA
M-x package-install RET buttercup RET
```

### Example

```elisp
(require 'buttercup)

(describe "String reversal"
  (it "reverses a simple string"
    (expect (reverse-string "hello") :to-equal "olleh"))

  (it "handles empty string"
    (expect (reverse-string "") :to-equal ""))

  (it "signals error for non-string"
    (expect (reverse-string 123) :to-throw 'error)))
```

### Running Buttercup

```bash
buttercup -L . test/
```

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        emacs_version: ['27.1', '28.1', '29.1']
    steps:
      - uses: actions/checkout@v2
      - uses: purcell/setup-emacs@master
        with:
          version: ${{ matrix.emacs_version }}
      - name: Run tests
        run: |
          emacs -batch -l ert -l my-package.el -l test/test-my-package.el \
          -f ert-run-tests-batch-and-exit
```

## Package Test Requirements

For MELPA packages:

```elisp
;; my-package.el
;;; Code:

(defun my-function (x)
  "Double X."
  (* x 2))

(provide 'my-package)
;;; my-package.el ends here
```

```elisp
;; test/test-my-package.el
;;; Code:

(require 'ert)
(require 'my-package)

(ert-deftest test-my-function ()
  (should (= (my-function 5) 10)))

(provide 'test-my-package)
;;; test-my-package.el ends here
```

## Resources

**ERT manual:** https://www.gnu.org/software/emacs/manual/html_node/ert/
**Buttercup:** https://github.com/jorgenschaefer/emacs-buttercup
**Package testing:** https://github.com/melpa/melpa/blob/master/CONTRIBUTING.md
