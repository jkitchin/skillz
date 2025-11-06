;;; my-mode.el --- Major mode for My Language  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <your.email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: languages
;; URL: https://github.com/user/my-mode

;;; Commentary:

;; A major mode for editing My Language files.
;;
;; Features:
;; - Syntax highlighting
;; - Indentation
;; - Comment handling
;; - Imenu support
;; - Eldoc integration
;;
;; To use, add to your init file:
;;
;;   (require 'my-mode)
;;   (add-to-list 'auto-mode-alist '("\\.my\\'" . my-mode))

;;; Code:

(require 'prog-mode)

;;; Customization

(defgroup my-mode nil
  "Major mode for My Language."
  :prefix "my-mode-"
  :group 'languages
  :link '(url-link :tag "GitHub" "https://github.com/user/my-mode"))

(defcustom my-mode-indent-offset 2
  "Number of spaces for each indentation level."
  :type 'integer
  :safe #'integerp
  :group 'my-mode)

(defcustom my-mode-tab-width 2
  "Width of a tab character."
  :type 'integer
  :safe #'integerp
  :group 'my-mode)

;;; Syntax Highlighting

(defconst my-mode-keywords
  '("if" "then" "else" "while" "for" "def" "class" "return" "import")
  "Keywords for My Language.")

(defconst my-mode-builtins
  '("print" "input" "len" "range" "map")
  "Built-in functions for My Language.")

(defconst my-mode-constants
  '("true" "false" "nil")
  "Constants in My Language.")

(defvar my-mode-font-lock-keywords
  `(
    ;; Comments: # to end of line
    ("#.*$" . font-lock-comment-face)

    ;; Strings: "..." and '...'
    ("\"\\([^\"\\]\\|\\\\.\\)*\"" . font-lock-string-face)
    ("'\\([^'\\]\\|\\\\.\\)*'" . font-lock-string-face)

    ;; Keywords
    (,(regexp-opt my-mode-keywords 'words) . font-lock-keyword-face)

    ;; Built-in functions
    (,(regexp-opt my-mode-builtins 'words) . font-lock-builtin-face)

    ;; Constants
    (,(regexp-opt my-mode-constants 'words) . font-lock-constant-face)

    ;; Function definitions: def function_name
    ("\\<def\\>[ \t]+\\([a-zA-Z_][a-zA-Z0-9_]*\\)"
     (1 font-lock-function-name-face))

    ;; Class definitions: class ClassName
    ("\\<class\\>[ \t]+\\([A-Z][a-zA-Z0-9_]*\\)"
     (1 font-lock-type-face))

    ;; Variables: assignments
    ("\\<\\([a-zA-Z_][a-zA-Z0-9_]*\\)[ \t]*="
     (1 font-lock-variable-name-face))

    ;; Numbers
    ("\\<[0-9]+\\(\\.[0-9]+\\)?\\>" . font-lock-constant-face)
    )
  "Font lock keywords for My Mode.")

;;; Syntax Table

(defvar my-mode-syntax-table
  (let ((table (make-syntax-table)))
    ;; Comments: # starts comment, newline ends it
    (modify-syntax-entry ?# "<" table)
    (modify-syntax-entry ?\n ">" table)

    ;; Strings
    (modify-syntax-entry ?\" "\"" table)
    (modify-syntax-entry ?\' "\"" table)

    ;; Operators
    (modify-syntax-entry ?+ "." table)
    (modify-syntax-entry ?- "." table)
    (modify-syntax-entry ?* "." table)
    (modify-syntax-entry ?/ "." table)
    (modify-syntax-entry ?% "." table)
    (modify-syntax-entry ?< "." table)
    (modify-syntax-entry ?> "." table)
    (modify-syntax-entry ?= "." table)

    ;; Parentheses and brackets
    (modify-syntax-entry ?\( "()" table)
    (modify-syntax-entry ?\) ")(" table)
    (modify-syntax-entry ?\[ "(]" table)
    (modify-syntax-entry ?\] ")[" table)
    (modify-syntax-entry ?\{ "(}" table)
    (modify-syntax-entry ?\} "){" table)

    ;; Underscore is part of word
    (modify-syntax-entry ?_ "w" table)

    table)
  "Syntax table for `my-mode'.")

;;; Indentation

(defun my-mode-indent-line ()
  "Indent current line for My Language."
  (interactive)
  (let ((indent-col 0)
        (cur-indent (current-indentation)))
    (save-excursion
      (beginning-of-line)
      (if (bobp)
          (setq indent-col 0)
        ;; Not at beginning of buffer
        (let ((not-indented t))
          ;; If current line closes a block
          (when (looking-at "^[ \t]*\\(end\\|else\\|elif\\)\\>")
            (save-excursion
              (forward-line -1)
              (setq indent-col (max 0 (- (current-indentation)
                                         my-mode-indent-offset)))))

          ;; Find previous non-empty line
          (when not-indented
            (save-excursion
              (while (and not-indented (not (bobp)))
                (forward-line -1)
                (unless (looking-at "^[ \t]*$")
                  ;; Previous line opens a block
                  (if (looking-at "^.*\\<\\(if\\|while\\|for\\|def\\|class\\)\\>.*:")
                      (progn
                        (setq indent-col (+ (current-indentation)
                                           my-mode-indent-offset))
                        (setq not-indented nil))
                    ;; Use same indentation as previous line
                    (setq indent-col (current-indentation))
                    (setq not-indented nil)))))))))

    ;; Apply indentation
    (if (<= (current-column) (current-indentation))
        (indent-line-to indent-col)
      (save-excursion
        (indent-line-to indent-col)))))

;;; Imenu Support

(defvar my-mode-imenu-generic-expression
  `(("Functions" "^[ \t]*def[ \t]+\\([a-zA-Z_][a-zA-Z0-9_]*\\)" 1)
    ("Classes" "^[ \t]*class[ \t]+\\([A-Z][a-zA-Z0-9_]*\\)" 1))
  "Imenu expression for My Mode.")

;;; Eldoc Support

(defun my-mode-eldoc-function ()
  "Return documentation for symbol at point."
  (let ((sym (thing-at-point 'symbol)))
    (when sym
      (cond
       ((string= sym "print")
        "print(value) - Print value to output")
       ((string= sym "len")
        "len(sequence) - Return length of sequence")
       ;; Add more built-in function docs here
       ))))

;;; Navigation

(defun my-mode-beginning-of-defun ()
  "Move to the beginning of the current function definition."
  (interactive)
  (re-search-backward "^[ \t]*def[ \t]+" nil 'move))

(defun my-mode-end-of-defun ()
  "Move to the end of the current function definition."
  (interactive)
  (re-search-forward "^[ \t]*def[ \t]+" nil 'move)
  (my-mode-beginning-of-defun)
  (forward-line 1)
  (let ((indent (current-indentation)))
    (while (and (not (eobp))
                (or (looking-at "^[ \t]*$")
                    (> (current-indentation) indent)))
      (forward-line 1))))

;;; Interactive Commands

;;;###autoload
(defun my-mode-run-file ()
  "Run the current My Language file."
  (interactive)
  (save-buffer)
  (let ((filename (buffer-file-name)))
    (if filename
        (compile (format "my-interpreter %s" (shell-quote-argument filename)))
      (user-error "Buffer is not visiting a file"))))

;;;###autoload
(defun my-mode-format-buffer ()
  "Format the current buffer according to My Language style."
  (interactive)
  (save-excursion
    (goto-char (point-min))
    (while (not (eobp))
      (unless (looking-at "^[ \t]*$")
        (my-mode-indent-line))
      (forward-line 1)))
  (message "Buffer formatted"))

;;; Keymap

(defvar my-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "C-c C-c") #'my-mode-run-file)
    (define-key map (kbd "C-c C-f") #'my-mode-format-buffer)
    (define-key map (kbd "C-M-a") #'my-mode-beginning-of-defun)
    (define-key map (kbd "C-M-e") #'my-mode-end-of-defun)
    map)
  "Keymap for `my-mode'.")

;;; Mode Definition

;;;###autoload
(define-derived-mode my-mode prog-mode "My"
  "Major mode for editing My Language files.

\\{my-mode-map}"
  :syntax-table my-mode-syntax-table

  ;; Comments
  (setq-local comment-start "# ")
  (setq-local comment-start-skip "#+\\s-*")
  (setq-local comment-end "")

  ;; Font lock
  (setq-local font-lock-defaults '(my-mode-font-lock-keywords))

  ;; Indentation
  (setq-local indent-line-function #'my-mode-indent-line)
  (setq-local tab-width my-mode-tab-width)

  ;; Imenu
  (setq-local imenu-generic-expression my-mode-imenu-generic-expression)

  ;; Eldoc
  (when (fboundp 'eldoc-mode)
    (setq-local eldoc-documentation-function #'my-mode-eldoc-function)
    (eldoc-mode 1))

  ;; Navigation
  (setq-local beginning-of-defun-function #'my-mode-beginning-of-defun)
  (setq-local end-of-defun-function #'my-mode-end-of-defun))

;;; Autoload

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.my\\'" . my-mode))

(provide 'my-mode)
;;; my-mode.el ends here
