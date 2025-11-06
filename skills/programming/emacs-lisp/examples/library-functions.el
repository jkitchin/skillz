;;; library-functions.el --- Example library functions  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <your.email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: examples

;;; Commentary:

;; This file demonstrates various common patterns and idioms in Emacs Lisp:
;;
;; - String manipulation
;; - List processing
;; - File operations
;; - Buffer manipulation
;; - Error handling
;; - Macro definitions
;; - Hook usage

;;; Code:

;;; String Functions

(defun library-capitalize-words (str)
  "Capitalize each word in STR.

Example:
  (library-capitalize-words \"hello world\")
  => \"Hello World\""
  (mapconcat #'capitalize (split-string str) " "))

(defun library-reverse-string (str)
  "Reverse the characters in STR."
  (concat (reverse (string-to-list str))))

(defun library-string-contains-p (str substring)
  "Return non-nil if STR contains SUBSTRING."
  (string-match-p (regexp-quote substring) str))

(defun library-trim-and-squeeze (str)
  "Trim whitespace and squeeze multiple spaces in STR."
  (let ((trimmed (string-trim str)))
    (replace-regexp-in-string "\\s-+" " " trimmed)))

;;; List Functions

(defun library-flatten (list)
  "Flatten nested LIST structure.

Example:
  (library-flatten '(1 (2 3) ((4))))
  => (1 2 3 4)"
  (cond
   ((null list) nil)
   ((listp (car list))
    (append (library-flatten (car list))
            (library-flatten (cdr list))))
   (t (cons (car list) (library-flatten (cdr list))))))

(defun library-partition (pred list)
  "Partition LIST into two lists based on PRED.

Returns cons cell (MATCHING . NON-MATCHING).

Example:
  (library-partition #'evenp '(1 2 3 4 5))
  => ((2 4) . (1 3 5))"
  (let (matching non-matching)
    (dolist (item list)
      (if (funcall pred item)
          (push item matching)
        (push item non-matching)))
    (cons (nreverse matching) (nreverse non-matching))))

(defun library-take (n list)
  "Return first N elements from LIST."
  (butlast list (- (length list) n)))

(defun library-group-by (key-fn list)
  "Group LIST elements by KEY-FN.

Returns alist where keys are results of KEY-FN and
values are lists of items with that key."
  (let (result)
    (dolist (item list)
      (let* ((key (funcall key-fn item))
             (existing (assoc key result)))
        (if existing
            (setcdr existing (cons item (cdr existing)))
          (push (list key item) result))))
    (mapcar (lambda (group)
              (cons (car group) (nreverse (cdr group))))
            (nreverse result))))

;;; File Functions

(defun library-read-lines (filename)
  "Read FILENAME and return list of lines."
  (with-temp-buffer
    (insert-file-contents filename)
    (split-string (buffer-string) "\n" t)))

(defun library-write-lines (filename lines)
  "Write LINES to FILENAME, one per line."
  (with-temp-file filename
    (insert (string-join lines "\n"))))

(defun library-file-size-human (filename)
  "Return human-readable size of FILENAME."
  (let* ((size (float (file-attribute-size (file-attributes filename))))
         (units '("B" "KB" "MB" "GB" "TB")))
    (while (and (> size 1024) (cdr units))
      (setq size (/ size 1024.0)
            units (cdr units)))
    (format "%.2f %s" size (car units))))

(defun library-directory-files-recursive (directory &optional match)
  "List all files in DIRECTORY recursively.

If MATCH is provided, only return files matching that regexp."
  (let (result)
    (dolist (file (directory-files directory t))
      (unless (member (file-name-nondirectory file) '("." ".."))
        (if (file-directory-p file)
            (setq result (append result
                                (library-directory-files-recursive file match)))
          (when (or (null match)
                    (string-match-p match file))
            (push file result)))))
    (nreverse result)))

;;; Buffer Functions

(defun library-buffer-string-no-properties ()
  "Return buffer contents without text properties."
  (buffer-substring-no-properties (point-min) (point-max)))

(defun library-count-matches-in-buffer (regexp)
  "Count matches of REGEXP in current buffer."
  (save-excursion
    (goto-char (point-min))
    (let ((count 0))
      (while (re-search-forward regexp nil t)
        (setq count (1+ count)))
      count)))

(defun library-replace-in-buffer (from to)
  "Replace all occurrences of FROM with TO in current buffer."
  (save-excursion
    (goto-char (point-min))
    (let ((count 0))
      (while (search-forward from nil t)
        (replace-match to nil t)
        (setq count (1+ count)))
      count)))

(defun library-sort-lines-in-region (beg end)
  "Sort lines in region from BEG to END."
  (interactive "r")
  (save-excursion
    (save-restriction
      (narrow-to-region beg end)
      (goto-char (point-min))
      (sort-subr nil #'forward-line #'end-of-line))))

;;; Functional Programming

(defun library-compose (&rest functions)
  "Return composition of FUNCTIONS.

Example:
  (funcall (library-compose #'1+ #'1+) 5)
  => 7"
  (lambda (x)
    (seq-reduce (lambda (acc fn) (funcall fn acc))
                (reverse functions)
                x)))

(defun library-partial (fn &rest args)
  "Return partial application of FN with ARGS."
  (lambda (&rest more-args)
    (apply fn (append args more-args))))

(defun library-memoize (fn)
  "Return memoized version of FN."
  (let ((cache (make-hash-table :test 'equal)))
    (lambda (&rest args)
      (or (gethash args cache)
          (puthash args (apply fn args) cache)))))

;;; Macros

(defmacro library-with-gensyms (syms &rest body)
  "Execute BODY with SYMS bound to unique symbols.

This is useful for writing hygienic macros."
  (declare (indent 1))
  `(let ,(mapcar (lambda (sym)
                   `(,sym (make-symbol ,(symbol-name sym))))
                 syms)
     ,@body))

(defmacro library-measure-time (&rest body)
  "Measure and return time taken to execute BODY."
  (declare (indent 0))
  (let ((start (make-symbol "start")))
    `(let ((,start (current-time)))
       (prog1 (progn ,@body)
         (message "Elapsed: %.3f seconds"
                  (float-time (time-since ,start)))))))

(defmacro library-with-temp-directory (var &rest body)
  "Execute BODY with VAR bound to a temporary directory.

The directory is automatically deleted after BODY completes."
  (declare (indent 1))
  (library-with-gensyms (dir)
    `(let ((,dir (make-temp-file "elisp-" t)))
       (unwind-protect
           (let ((,var ,dir))
             ,@body)
         (when (file-exists-p ,dir)
           (delete-directory ,dir t))))))

;;; Error Handling

(defun library-safe-divide (a b)
  "Safely divide A by B.

Returns result or nil if division fails."
  (condition-case nil
      (/ a b)
    (arith-error nil)))

(defun library-retry (fn max-attempts delay)
  "Retry FN up to MAX-ATTEMPTS times with DELAY seconds between attempts.

Returns result of FN or signals the last error."
  (let ((attempt 0)
        result)
    (while (< attempt max-attempts)
      (condition-case err
          (progn
            (setq result (funcall fn))
            (setq attempt max-attempts))
        (error
         (setq attempt (1+ attempt))
         (when (< attempt max-attempts)
           (sleep-for delay))
         (when (= attempt max-attempts)
           (signal (car err) (cdr err))))))
    result))

;;; Data Structures

(defun library-make-stack ()
  "Create a new stack."
  (list 'stack nil))

(defun library-stack-push (stack item)
  "Push ITEM onto STACK."
  (setcdr stack (cons item (cdr stack))))

(defun library-stack-pop (stack)
  "Pop and return top item from STACK."
  (when (cdr stack)
    (prog1 (cadr stack)
      (setcdr stack (cddr stack)))))

(defun library-stack-peek (stack)
  "Return top item from STACK without removing it."
  (cadr stack))

(defun library-stack-empty-p (stack)
  "Return non-nil if STACK is empty."
  (null (cdr stack)))

;;; Predicates

(defun library-positive-integer-p (n)
  "Return non-nil if N is a positive integer."
  (and (integerp n) (> n 0)))

(defun library-non-empty-string-p (s)
  "Return non-nil if S is a non-empty string."
  (and (stringp s) (not (string-empty-p s))))

(defun library-valid-email-p (email)
  "Return non-nil if EMAIL looks like a valid email address."
  (and (stringp email)
       (string-match-p "^[^@]+@[^@]+\\.[^@]+$" email)))

;;; Interactive Utilities

;;;###autoload
(defun library-insert-date ()
  "Insert current date at point."
  (interactive)
  (insert (format-time-string "%Y-%m-%d")))

;;;###autoload
(defun library-count-words-region (beg end)
  "Count words in region from BEG to END."
  (interactive "r")
  (let ((count (how-many "\\w+" beg end)))
    (message "Region has %d words" count)))

;;;###autoload
(defun library-duplicate-line ()
  "Duplicate the current line."
  (interactive)
  (let ((line (buffer-substring (line-beginning-position)
                                (line-end-position))))
    (end-of-line)
    (newline)
    (insert line)))

;;; Hook Utilities

(defun library-add-hooks (hooks function)
  "Add FUNCTION to all HOOKS."
  (dolist (hook hooks)
    (add-hook hook function)))

(defun library-remove-hooks (hooks function)
  "Remove FUNCTION from all HOOKS."
  (dolist (hook hooks)
    (remove-hook hook function)))

;;; Usage Examples

(defun library-demonstrate-usage ()
  "Demonstrate usage of library functions."
  (interactive)
  (with-output-to-temp-buffer "*Library Demo*"
    (princ "=== String Functions ===\n\n")
    (princ (format "Capitalize: %s\n"
                   (library-capitalize-words "hello world")))
    (princ (format "Reverse: %s\n"
                   (library-reverse-string "hello")))

    (princ "\n=== List Functions ===\n\n")
    (princ (format "Flatten: %s\n"
                   (library-flatten '(1 (2 3) ((4))))))
    (princ (format "Partition evens: %s\n"
                   (library-partition #'evenp '(1 2 3 4 5))))

    (princ "\n=== Functional Programming ===\n\n")
    (let ((add-two (library-compose #'1+ #'1+)))
      (princ (format "Composed function (add 2): %s\n"
                     (funcall add-two 5))))

    (princ "\n=== Stack Operations ===\n\n")
    (let ((stack (library-make-stack)))
      (library-stack-push stack 1)
      (library-stack-push stack 2)
      (library-stack-push stack 3)
      (princ (format "Pop: %s\n" (library-stack-pop stack)))
      (princ (format "Peek: %s\n" (library-stack-peek stack))))))

(provide 'library-functions)
;;; library-functions.el ends here
