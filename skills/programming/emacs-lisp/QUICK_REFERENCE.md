# Emacs Lisp Quick Reference

A concise reference for Emacs Lisp patterns, naming conventions, and common idioms based on community best practices.

## File Header Template

```elisp
;;; package-name.el --- Brief description  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <email@example.com>
;; Maintainer: Your Name <email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: convenience, tools
;; URL: https://github.com/user/package-name

;;; Commentary:

;; Longer description of what the package does.
;; Usage examples and configuration tips.

;;; Code:

;; Your code here

(provide 'package-name)
;;; package-name.el ends here
```

## Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Function | `package-function-name` | `projectile-find-file` |
| Variable | `package-variable-name` | `projectile-project-root` |
| Constant | `package-constant-name` | `package-version` |
| Private function | `package--private-function` | `projectile--get-root` |
| Unused variable | `_variable` | `(lambda (x _y) x)` |
| Single-word predicate | `wordp` | `evenp`, `emptyp` |
| Multi-word predicate | `phrase-p` | `buffer-live-p` |
| Face | `package-face-name` | `font-lock-keyword-face` |
| Group | `package` | `projectile` |

## Indentation Rules

```elisp
;; Standard function call - 2 spaces
(function-name arg1
  arg2
  arg3)

;; Special forms - 4 spaces for special args, 2 for body
(if condition
    then-clause    ; 4 spaces
  else-clause)     ; 2 spaces

(let ((var1 val1)
      (var2 val2)) ; 4 spaces for bindings
  body-form-1)     ; 2 spaces for body

(defun name (args)
  "Docstring."
  body)            ; 2 spaces
```

## Common Patterns

### Defining Variables

```elisp
;; Simple variable
(defvar my-package-variable "default"
  "Description of variable.")

;; Customizable variable
(defcustom my-package-option t
  "Description of option."
  :type 'boolean
  :group 'my-package)

;; Constant (by convention)
(defconst my-package-version "1.0.0"
  "Package version.")

;; Buffer-local variable
(defvar-local my-package-local-var nil
  "Buffer-local variable.")
```

### Defining Functions

```elisp
;; Basic function
(defun my-package-function (arg1 arg2)
  "Process ARG1 and ARG2.

Returns the processed result."
  (+ arg1 arg2))

;; Function with optional arguments
(defun my-package-function (required &optional opt1 opt2)
  "Process REQUIRED with optional OPT1 and OPT2."
  (list required opt1 opt2))

;; Function with rest arguments
(defun my-package-function (first &rest rest)
  "Process FIRST and REST arguments."
  (cons first rest))

;; Function with keyword arguments (requires cl-lib)
(cl-defun my-package-function (required &key opt1 opt2)
  "Process REQUIRED with optional :opt1 and :opt2."
  (list required opt1 opt2))
```

### Interactive Commands

```elisp
;; Simple interactive
(defun my-package-command ()
  "Do something interactively."
  (interactive)
  (message "Command executed"))

;; With arguments
(defun my-package-insert (text)
  "Insert TEXT at point."
  (interactive "sText to insert: ")
  (insert text))

;; With region
(defun my-package-process-region (start end)
  "Process region from START to END."
  (interactive "r")
  (save-excursion
    (goto-char start)
    ;; process...
    ))

;; With prefix argument
(defun my-package-insert-n (n)
  "Insert N copies."
  (interactive "p")
  (dotimes (_ n)
    (insert "text")))

;; Complex interactive form
(defun my-package-replace (from to)
  "Replace FROM with TO."
  (interactive
   (list (read-string "From: ")
         (read-string "To: ")))
  ;; implementation
  )
```

## Control Flow

### Conditionals

```elisp
;; if - two branches
(if condition
    then-form
  else-form)

;; when - single branch (if with implicit progn)
(when condition
  (do-something)
  (do-more))

;; unless - negated condition
(unless condition
  (do-something))

;; cond - multiple conditions
(cond
 ((= x 1) (handle-one))
 ((= x 2) (handle-two))
 (t (handle-default)))

;; pcase - pattern matching
(pcase value
  ('symbol (handle-symbol))
  ((pred numberp) (handle-number))
  (`(,car . ,cdr) (handle-cons))
  (_ (handle-default)))
```

### Loops

```elisp
;; dolist - iterate over list
(dolist (item items)
  (process item))

;; dotimes - numeric iteration
(dotimes (i 10)
  (insert (format "%d " i)))

;; while - general loop
(while condition
  (do-something))

;; cl-loop - powerful iteration
(cl-loop for x in items
         when (> x 5)
         collect x)
```

## Working with Lists

```elisp
;; Create list
(list 1 2 3)
'(1 2 3)

;; Cons cell
(cons 1 2)           ; (1 . 2)
(cons 1 '(2 3))      ; (1 2 3)

;; First/rest
(car '(1 2 3))       ; 1
(cdr '(1 2 3))       ; (2 3)

;; Append lists
(append '(1 2) '(3 4))  ; (1 2 3 4)

;; Map over list
(mapcar #'1+ '(1 2 3))  ; (2 3 4)

;; Filter list
(seq-filter #'evenp '(1 2 3 4))  ; (2 4)

;; Reduce list
(seq-reduce #'+ '(1 2 3 4) 0)  ; 10

;; Association list (alist)
(setq alist '((foo . 1) (bar . 2)))
(alist-get 'foo alist)  ; 1

;; Property list (plist)
(setq plist '(:foo 1 :bar 2))
(plist-get plist :foo)  ; 1
```

## Working with Strings

```elisp
;; Concatenate
(concat "hello" " " "world")

;; Format string
(format "Value: %s, Number: %d" "foo" 42)

;; String comparison
(string= "foo" "foo")
(string-match "regexp" "string")

;; Split/join
(split-string "foo,bar,baz" ",")
(string-join '("foo" "bar") ",")

;; Case conversion
(upcase "foo")
(downcase "BAR")
(capitalize "foo")
```

## Working with Buffers

```elisp
;; Current buffer
(current-buffer)

;; Switch buffer temporarily
(with-current-buffer buffer-or-name
  ;; operations in that buffer
  )

;; Create temporary buffer
(with-temp-buffer
  (insert "content")
  (buffer-string))

;; Save excursion (restore point and mark)
(save-excursion
  (goto-char (point-min))
  ;; do something
  )

;; Save restriction (restore narrowing)
(save-restriction
  (narrow-to-region start end)
  ;; do something
  )
```

## Working with Files

```elisp
;; Read file
(with-temp-buffer
  (insert-file-contents filename)
  (buffer-string))

;; Write file
(with-temp-file filename
  (insert content))

;; Check file exists
(file-exists-p filename)

;; Get directory name
(file-name-directory "/path/to/file.el")  ; "/path/to/"

;; Get file name
(file-name-nondirectory "/path/to/file.el")  ; "file.el"

;; Expand file name
(expand-file-name "~/.emacs.d/init.el")
```

## Hooks

```elisp
;; Add to hook (use function quote)
(add-hook 'text-mode-hook #'turn-on-auto-fill)

;; Remove from hook
(remove-hook 'text-mode-hook #'turn-on-auto-fill)

;; Local hook
(add-hook 'text-mode-hook #'my-function nil t)

;; Run hooks
(run-hooks 'my-mode-hook)
```

## Keybindings

```elisp
;; Global keybinding
(global-set-key (kbd "C-c f") #'find-file)

;; Mode-specific keybinding
(define-key my-mode-map (kbd "C-c C-c") #'my-command)

;; Key sequence
(global-set-key (kbd "C-c C-x C-f") #'my-command)

;; Unbind key
(global-unset-key (kbd "C-x C-c"))
```

## Defining Modes

### Minor Mode

```elisp
;;;###autoload
(define-minor-mode my-mode
  "Toggle my-mode."
  :lighter " MyMode"
  :keymap (let ((map (make-sparse-keymap)))
            (define-key map (kbd "C-c m") #'my-command)
            map)
  (if my-mode
      (my-mode--enable)
    (my-mode--disable)))
```

### Major Mode

```elisp
(defvar my-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "C-c C-c") #'my-compile)
    map)
  "Keymap for `my-mode'.")

;;;###autoload
(define-derived-mode my-mode prog-mode "My"
  "Major mode for editing my files."
  (setq-local comment-start "#")
  (setq-local comment-start-skip "#+\\s-*"))
```

## Macros

```elisp
;; Basic macro
(defmacro with-my-setup (&rest body)
  "Execute BODY with my setup."
  (declare (indent 0))
  `(let ((saved-value my-var))
     (unwind-protect
         (progn ,@body)
       (setq my-var saved-value))))

;; Macro with debug spec
(defmacro my-when (condition &rest body)
  "Like `when' but with special handling."
  (declare (indent 1)
           (debug (form body)))
  `(if ,condition
       (progn ,@body)))
```

## Error Handling

```elisp
;; Signal error
(error "Error message: %s" value)

;; User error (doesn't show backtrace)
(user-error "Invalid input")

;; Catch errors
(condition-case err
    (risky-operation)
  (file-error
   (message "File error: %s" (error-message-string err)))
  (error
   (message "General error: %s" (error-message-string err))))

;; Ignore errors
(ignore-errors
  (risky-operation))
```

## Advice

```elisp
;; Add advice before function
(advice-add 'function-name :before #'my-advice-function)

;; Add advice after function
(advice-add 'function-name :after #'my-advice-function)

;; Wrap function (around advice)
(defun my-advice-function (orig-fun &rest args)
  "Advice for function."
  (message "Before")
  (let ((result (apply orig-fun args)))
    (message "After")
    result))

(advice-add 'function-name :around #'my-advice-function)

;; Remove advice
(advice-remove 'function-name #'my-advice-function)
```

## Common Idioms

### Check if feature loaded

```elisp
(when (featurep 'package-name)
  (do-something))
```

### Load conditionally

```elisp
(with-eval-after-load 'package-name
  (setq package-name-option t))
```

### Optional require

```elisp
(when (require 'optional-package nil t)
  (use-optional-package))
```

### Get user input

```elisp
(read-string "Prompt: ")
(read-number "Number: ")
(yes-or-no-p "Really?")
(y-or-n-p "Quick question?")
(completing-read "Choose: " '("option1" "option2"))
```

### Message vs. print

```elisp
;; For user messages (shows in echo area)
(message "Status: %s" value)

;; For debugging (returns value)
(prin1 value)

;; Insert in buffer
(insert (format "Text: %s" value))
```

## Interactive Codes Quick Reference

| Code | Meaning |
|------|---------|
| `p` | Prefix argument as number |
| `P` | Raw prefix argument |
| `r` | Region (start and end) |
| `s` | String (read from minibuffer) |
| `f` | Existing file name |
| `F` | File name (may not exist) |
| `d` | Directory name |
| `b` | Buffer name |
| `B` | Buffer name (may not exist) |
| `n` | Number |
| `c` | Character |
| `k` | Key sequence |
| `x` | Lisp expression |

## Best Practices Checklist

- [ ] File header includes lexical-binding declaration
- [ ] All public symbols prefixed with package name
- [ ] Private functions use double-hyphen prefix
- [ ] Interactive commands have autoload cookies
- [ ] All functions have docstrings
- [ ] First line of docstring is complete sentence
- [ ] Custom variables use `defcustom` not `defvar`
- [ ] Hooks use function quotes (`#'function`)
- [ ] Predicates end in `p` or `-p`
- [ ] File ends with `(provide 'package-name)` and footer
- [ ] No hard tabs, only spaces
- [ ] Line length under 80 characters where reasonable
- [ ] All closing parens on same line
