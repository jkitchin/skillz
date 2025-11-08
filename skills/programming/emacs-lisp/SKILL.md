---
name: emacs-lisp
description: |
  Expert guidance for writing professional Emacs Lisp code following the
  community-driven Emacs Lisp Style Guide by Bozhidar Batsov. Covers layout,
  naming conventions, syntax preferences, macros, documentation, and best practices.

  Use this skill when writing Emacs configuration, creating Emacs packages,
  writing interactive commands, developing major/minor modes, or refactoring
  Emacs Lisp code.

  Emphasizes: lexical binding, proper naming, documentation standards, autoloading,
  and community conventions.
allowed-tools:
  - Read
  - Edit
  - Write
  - Bash
  - Grep
  - Glob
---

# Emacs Lisp Style Guide Skill

This skill provides expert guidance for writing professional, maintainable Emacs Lisp code following the community-driven Emacs Lisp Style Guide.

## When to Use This Skill

Use this skill when:
- Writing Emacs configuration files (~/.emacs.d/init.el)
- Creating Emacs packages or libraries
- Developing major or minor modes
- Writing interactive commands and functions
- Creating keybindings and hooks
- Refactoring existing Emacs Lisp code
- Setting up package metadata and autoloads
- Writing Emacs Lisp macros
- Creating custom variables and faces

## Core Principles

### 1. Source Code Layout

#### Indentation and Spacing

**Use Spaces, Never Tabs:**
```elisp
;; Good - spaces for indentation
(defun my-function (arg1 arg2)
  (let ((result (+ arg1 arg2)))
    (message "Result: %s" result)))

;; Bad - mixing tabs and spaces (never do this)
```

**Align Function Arguments Vertically:**
```elisp
;; Good - arguments aligned
(defun long-function-name (arg1
                           arg2
                           arg3)
  (body))

;; Good - first arg on new line, aligned with function name
(defun long-function-name
    (arg1 arg2 arg3)
  (body))

;; Bad - inconsistent alignment
(defun long-function-name (arg1
    arg2
  arg3)
  (body))
```

**Special Form Indentation:**
```elisp
;; Special forms use 4-space indent for special arguments,
;; 2-space for body

;; Good - if with proper indentation
(if some-condition
    (do-something)  ; 4 spaces for then-clause
  (do-else))        ; 2 spaces for else-clause

;; Good - let bindings
(let ((var1 value1)
      (var2 value2))  ; 4 spaces for bindings
  (body-form-1)       ; 2 spaces for body
  (body-form-2))

;; Good - when (no special argument)
(when condition
  (do-something)
  (do-more))
```

**Maximum Line Length:**
```elisp
;; Aim for 80 characters per line where reasonable
;; Break long lines logically

;; Good
(defun my-function ()
  "Short description."
  (let ((long-variable-name
         (some-function-that-returns-value)))
    long-variable-name))

;; Acceptable for long strings
(message "This is a very long message that exceeds 80 characters but is kept on one line for readability")
```

**Parentheses:**
```elisp
;; Good - space between text and opening paren
(defun my-function (arg)
  (list arg))

;; Bad - no space after function name
(defun my-function(arg)
  (list arg))

;; Good - all closing parens on same line
(defun my-function ()
  (when condition
    (let ((x 1))
      (+ x 2))))

;; Bad - closing parens on separate lines
(defun my-function ()
  (when condition
    (let ((x 1))
      (+ x 2)
    )
  )
)
```

**Vertical Whitespace:**
```elisp
;; Good - empty line between top-level forms
(defvar my-var1 "value1")

(defvar my-var2 "value2")

(defun my-function ()
  "Do something.")

;; Exception: related definitions can be grouped
(defvar my-var1 "value1")
(defvar my-var2 "value2")
(defvar my-var3 "value3")
```

### 2. Syntax and Idioms

#### Prefer Higher-Level Constructs

**Use `when` Instead of Single-Branch `if`:**
```elisp
;; Good
(when condition
  (do-something)
  (do-more))

;; Bad
(if condition
    (progn
      (do-something)
      (do-more)))
```

**Use `unless` for Negated Conditions:**
```elisp
;; Good
(unless condition
  (do-something))

;; Bad
(when (not condition)
  (do-something))

;; Bad
(if (not condition)
    (do-something))
```

**Use `not` Instead of `null`:**
```elisp
;; Good - checking for nil/false
(not value)

;; Good - checking for empty list specifically
(null my-list)

;; Bad - using null for general nil check
(null value)
```

**Leverage Comparison Functions:**
```elisp
;; Good - direct comparison
(< 5 x 10)

;; Bad - nested conditions
(and (> x 5) (< x 10))

;; Good - multiple comparisons
(= x y z)

;; Bad - nested comparisons
(and (= x y) (= y z))
```

**Use `cond` with `t` as Catch-All:**
```elisp
;; Good
(cond
 ((eq x 'foo) (handle-foo))
 ((eq x 'bar) (handle-bar))
 (t (handle-default)))

;; Bad - else at the end
(cond
 ((eq x 'foo) (handle-foo))
 ((eq x 'bar) (handle-bar))
 (else (handle-default)))
```

**Increment and Decrement:**
```elisp
;; Good
(1+ x)
(1- x)

;; Bad
(+ x 1)
(- x 1)
```

**Prefer `with-eval-after-load`:**
```elisp
;; Good - modern approach
(with-eval-after-load 'company
  (setq company-idle-delay 0.2))

;; Acceptable but older
(eval-after-load 'company
  '(setq company-idle-delay 0.2))
```

### 3. Naming Conventions

#### General Naming Rules

**Use Kebab-Case (Lisp-Case):**
```elisp
;; Good
(defun my-function-name ()
  (let ((my-variable 10))
    my-variable))

;; Bad - using snake_case
(defun my_function_name ()
  (let ((my_variable 10))
    my_variable))

;; Bad - using camelCase
(defun myFunctionName ()
  (let ((myVariable 10))
    myVariable))
```

**Library Name Prefixes:**
```elisp
;; For a library named 'my-package', prefix all public symbols
(defun my-package-do-something ()
  "Public function.")

(defvar my-package-option t
  "Public variable.")

;; Private functions use double-hyphen
(defun my-package--internal-helper ()
  "Private function.")

;; Example from projectile package
(defun projectile-find-file ()
  "Public command.")

(defun projectile--get-project-root ()
  "Private helper.")
```

**Unused Variables:**
```elisp
;; Good - prefix with underscore
(lambda (x _y _z)
  (* x 2))

;; Good - in destructuring
(pcase-let ((`(,first . ,_rest) my-list))
  first)
```

#### Predicates

**Single-Word Predicates:**
```elisp
;; Good - ends with 'p'
(defun evenp (n)
  (zerop (mod n 2)))

(defun emptyp (list)
  (null list))
```

**Multi-Word Predicates:**
```elisp
;; Good - ends with '-p'
(defun buffer-live-p (buffer)
  (and buffer (buffer-name buffer)))

(defun my-package-valid-state-p (state)
  (member state '(ready active complete)))
```

#### Special Naming

**Face Names:**
```elisp
;; Good - no -face suffix
(defface my-package-error
  '((t :foreground "red"))
  "Face for error messages.")

;; Bad - redundant -face suffix
(defface my-package-error-face
  '((t :foreground "red"))
  "Face for error messages.")
```

### 4. Functions

#### Lambda Functions

**Usage Guidelines:**
```elisp
;; Good - lambda for local use
(mapcar (lambda (x) (* x 2)) '(1 2 3))

;; Bad - lambda in hooks (use named functions)
(add-hook 'emacs-lisp-mode-hook
          (lambda () (eldoc-mode 1)))

;; Good - named function in hooks
(defun my-elisp-mode-setup ()
  "Setup for Emacs Lisp mode."
  (eldoc-mode 1))

(add-hook 'emacs-lisp-mode-hook #'my-elisp-mode-setup)
```

**Never Hard-Quote Lambdas:**
```elisp
;; Bad - hard-quoted lambda
(mapcar '(lambda (x) (* x 2)) list)

;; Good - unquoted lambda
(mapcar (lambda (x) (* x 2)) list)

;; Good - function quote (sharp-quote)
(mapcar #'(lambda (x) (* x 2)) list)
```

**Don't Wrap Existing Functions:**
```elisp
;; Bad - unnecessary wrapper
(add-hook 'prog-mode-hook
          (lambda () (linum-mode 1)))

;; Good - direct function reference
(add-hook 'prog-mode-hook #'linum-mode)
```

#### Function Quotes

**Use Sharp-Quote for Function Names:**
```elisp
;; Good - function quote aids byte-compilation
(mapcar #'upcase '("foo" "bar"))
(add-hook 'text-mode-hook #'turn-on-auto-fill)
(global-set-key (kbd "C-c f") #'find-file)

;; Acceptable but less optimal
(mapcar 'upcase '("foo" "bar"))

;; Bad - unquoted (only works for interactive commands)
(global-set-key (kbd "C-c f") 'find-file)
```

#### Parameter Limits

**Avoid Too Many Parameters:**
```elisp
;; Bad - too many positional parameters
(defun create-user (name email age address city state zip country)
  ...)

;; Good - use plist or alist for many parameters
(defun create-user (name email &rest properties)
  "Create user with NAME and EMAIL.
PROPERTIES is a plist that may include :age, :address, :city, etc."
  (let ((age (plist-get properties :age))
        (address (plist-get properties :address)))
    ...))

;; Good - use cl-defun with keyword arguments
(cl-defun create-user (name email &key age address city state)
  "Create user with NAME and EMAIL and optional properties."
  ...)
```

### 5. Macros

#### When to Use Macros

**Only Use Macros When Functions Won't Work:**
```elisp
;; Good - macro needed for special evaluation
(defmacro with-timing (label &rest body)
  "Execute BODY and print execution time with LABEL."
  (declare (indent 1))
  (let ((start (make-symbol "start")))
    `(let ((,start (current-time)))
       (prog1 (progn ,@body)
         (message "%s: %.2fs" ,label
                  (float-time (time-since ,start)))))))

;; Bad - function would work fine here
(defmacro double (x)
  `(* 2 ,x))

;; Good - use function instead
(defun double (x)
  (* 2 x))
```

#### Macro Design

**Extract Examples First:**
```elisp
;; Before writing macro, write example usage:
;;
;; (with-temporary-buffer
;;   (insert "content")
;;   (buffer-string))
;;
;; Then implement to match desired usage
```

**Use Declare for Metadata:**
```elisp
;; Good - complete macro with declare
(defmacro with-my-mode (buffer &rest body)
  "Execute BODY with my-mode enabled in BUFFER."
  (declare (indent 1)
           (debug (sexp body)))
  `(with-current-buffer ,buffer
     (my-mode 1)
     (unwind-protect
         (progn ,@body)
       (my-mode -1))))
```

**Prefer Syntax Quoting (Backquote):**
```elisp
;; Good - backquote with unquote
(defmacro my-when (condition &rest body)
  `(if ,condition
       (progn ,@body)))

;; Bad - manual list construction
(defmacro my-when (condition &rest body)
  (list 'if condition
        (cons 'progn body)))
```

### 6. Documentation

#### Docstrings

**Function Docstrings:**
```elisp
;; Good - complete docstring
(defun my-package-process-file (filename)
  "Process the file at FILENAME.

The file is read, processed, and the result is returned.
If FILENAME doesn't exist, signal an error.

Returns the processed content as a string."
  (with-temp-buffer
    (insert-file-contents filename)
    (buffer-string)))

;; Good - document all parameters
(defun my-package-create-entry (name email &optional age)
  "Create entry with NAME and EMAIL.

NAME should be a non-empty string.
EMAIL should be a valid email address.
Optional AGE should be a positive integer.

Returns the created entry as a plist."
  ...)
```

**Variable Docstrings:**
```elisp
;; Good
(defvar my-package-timeout 30
  "Timeout in seconds for network operations.

This value is used by all network functions in my-package.
Set to nil to disable timeout.")

;; Good - custom variable with type
(defcustom my-package-auto-save t
  "Whether to auto-save files.

When non-nil, files are automatically saved after modifications."
  :type 'boolean
  :group 'my-package)
```

**First Line Summary:**
```elisp
;; Good - first line is complete sentence
(defun my-function (arg)
  "Process ARG and return result.

Additional details about the function behavior.
More explanation if needed."
  ...)

;; Bad - incomplete first line
(defun my-function (arg)
  "Process ARG and
return result."
  ...)
```

#### Comments

**Comment Levels:**
```elisp
;;; Section Heading
;;; Used for major sections of the file

;;; Commentary:
;; Package description and usage examples

;;; Code:

;; Top-level comment
;; Used for explanations before functions

(defun my-function ()
  ;; Comment within function
  ;; Used for explaining code blocks
  (let ((x 10))  ; Margin comment for inline explanation
    x))
```

**Comment Best Practices:**
```elisp
;; Good - explains why, not what
(defun my-function (items)
  ;; Sort items first to optimize lookup performance
  (let ((sorted-items (sort items #'string<)))
    (process-items sorted-items)))

;; Bad - redundant comment
(defun my-function (items)
  ;; Sort the items
  (sort items #'string<))

;; Good - self-documenting code instead of comments
(defun process-active-buffers ()
  (let ((active-buffers (get-active-buffers)))
    (mapc #'process-buffer active-buffers)))

;; Bad - needs comments to explain
(defun process ()
  ;; Get buffers that are active
  (let ((bufs (get-bufs)))
    ;; Process each one
    (mapc #'proc bufs)))
```

**Annotation Keywords:**
```elisp
;; TODO: Add support for remote files
;; FIXME: This breaks with Unicode characters
;; OPTIMIZE: Cache results for better performance
;; HACK: Workaround for bug in package.el
;; REVIEW: Is this the best approach?

;; Good - with initials and date
;; TODO(user 2025-01-15): Implement async version
;; FIXME(user 2025-01-15): Handle edge case with empty input
```

### 7. Loading and Autoloading

#### Module Structure

**Proper File Footer:**
```elisp
;;; my-package.el --- Short description  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <your.email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: convenience, tools
;; URL: https://github.com/user/my-package

;;; Commentary:

;; Longer description of what the package does.
;; Usage examples.

;;; Code:

(defun my-package-do-something ()
  "Do something useful.")

(provide 'my-package)
;;; my-package.el ends here
```

#### Require vs Load

**Use `require` for Dependencies:**
```elisp
;; Good - idempotent, safe
(require 'subr-x)
(require 'cl-lib)

;; Bad - not idempotent
(load "subr-x")
(load-library "cl-lib")
```

#### Autoload Cookies

**Use Autoload for User-Facing Commands:**
```elisp
;; Good - autoload interactive command
;;;###autoload
(defun my-package-start ()
  "Start my-package."
  (interactive)
  (my-package-mode 1))

;; Good - autoload mode definition
;;;###autoload
(define-minor-mode my-package-mode
  "Toggle my-package mode."
  :lighter " MyPkg"
  (if my-package-mode
      (my-package--enable)
    (my-package--disable)))

;; Bad - autoloading internal function
;;;###autoload
(defun my-package--internal-helper ()
  "Internal helper.")
```

### 8. Lexical Binding

**Always Enable Lexical Binding:**
```elisp
;;; my-package.el --- Description  -*- lexical-binding: t; -*-

;; Lexical binding is faster and more intuitive
;; Always use it for new code
```

**Benefits of Lexical Binding:**
```elisp
;; With lexical binding, closures work correctly
(defun make-counter ()
  "Create a counter function."
  (let ((count 0))
    (lambda ()
      (setq count (1+ count))
      count)))

;; Usage:
;; (setq counter1 (make-counter))
;; (funcall counter1) => 1
;; (funcall counter1) => 2
```

### 9. List Operations

**Use `dolist` for Iteration:**
```elisp
;; Good - clear and idiomatic
(dolist (item items)
  (process-item item))

;; Bad - verbose
(mapc (lambda (item)
        (process-item item))
      items)
```

**Don't Use `mapcar` for Side Effects:**
```elisp
;; Bad - mapcar creates unused list
(mapcar #'process-item items)

;; Good - dolist for side effects
(dolist (item items)
  (process-item item))

;; Good - seq-do for functional style
(seq-do #'process-item items)
```

### 10. Interactive Commands

**Proper Interactive Specifications:**
```elisp
;; Good - with prompt
(defun my-package-open-file (filename)
  "Open FILENAME in my-package."
  (interactive "FOpen file: ")
  (find-file filename))

;; Good - multiple arguments
(defun my-package-replace (from to)
  "Replace FROM with TO in current buffer."
  (interactive
   (list (read-string "Replace: ")
         (read-string "With: ")))
  (save-excursion
    (goto-char (point-min))
    (while (search-forward from nil t)
      (replace-match to))))

;; Good - using prefix argument
(defun my-package-insert-n (n)
  "Insert N copies of something."
  (interactive "p")
  (dotimes (_ n)
    (insert "text")))
```

### 11. Error Handling

**Signal Errors Appropriately:**
```elisp
;; Good - descriptive error
(defun my-package-get-user (id)
  "Get user by ID."
  (or (gethash id my-package-users)
      (error "User not found: %s" id)))

;; Good - user-error for user mistakes
(defun my-package-save ()
  "Save current data."
  (interactive)
  (unless (my-package--has-changes-p)
    (user-error "No changes to save")))

;; Good - with condition-case
(defun my-package-load-file (filename)
  "Load FILENAME safely."
  (condition-case err
      (with-temp-buffer
        (insert-file-contents filename)
        (read (current-buffer)))
    (file-error
     (message "Cannot read file: %s" (error-message-string err))
     nil)))
```

## Best Practices Summary

When writing Emacs Lisp code:

1. **Always enable lexical binding** in file header
2. **Prefix all public symbols** with package name
3. **Use `defcustom`** for user-configurable variables
4. **Add autoload cookies** to interactive commands
5. **Write comprehensive docstrings** for all public functions
6. **Prefer built-in functions** over reinventing functionality
7. **Use `declare`** in macros for indentation and debugging
8. **Handle errors gracefully** with informative messages
9. **Follow naming conventions** consistently
10. **Keep functions focused** on single responsibilities

## Common Patterns

### Package Template

```elisp
;;; my-package.el --- Brief description  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: convenience
;; URL: https://github.com/user/my-package

;;; Commentary:

;; Longer description and usage examples.

;;; Code:

(require 'cl-lib)

(defgroup my-package nil
  "My package customization group."
  :group 'convenience
  :prefix "my-package-")

(defcustom my-package-option t
  "Enable some option."
  :type 'boolean
  :group 'my-package)

;;;###autoload
(defun my-package-command ()
  "Main command."
  (interactive)
  (message "Hello from my-package!"))

(provide 'my-package)
;;; my-package.el ends here
```

## Resources

- [Emacs Lisp Style Guide](https://github.com/bbatsov/emacs-lisp-style-guide)
- [Emacs Lisp Reference Manual](https://www.gnu.org/software/emacs/manual/html_node/elisp/)
- [Emacs Lisp Packaging Guide](https://www.gnu.org/software/emacs/manual/html_node/elisp/Packaging.html)
