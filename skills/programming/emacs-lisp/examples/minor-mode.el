;;; highlight-numbers.el --- Minor mode example  -*- lexical-binding: t; -*-

;; Copyright (C) 2025 Your Name

;; Author: Your Name <your.email@example.com>
;; Version: 1.0.0
;; Package-Requires: ((emacs "27.1"))
;; Keywords: convenience, faces
;; URL: https://github.com/user/highlight-numbers

;;; Commentary:

;; A simple minor mode that highlights numbers in the current buffer.
;;
;; This serves as an example of how to structure a minor mode with:
;; - Customization options
;; - Faces
;; - Automatic highlighting
;; - Enable/disable functionality
;; - Hook integration
;;
;; Usage:
;;
;;   M-x highlight-numbers-mode
;;
;; Or enable globally:
;;
;;   (add-hook 'prog-mode-hook #'highlight-numbers-mode)

;;; Code:

;;; Customization

(defgroup highlight-numbers nil
  "Highlight numbers in buffers."
  :prefix "highlight-numbers-"
  :group 'convenience)

(defcustom highlight-numbers-patterns
  '("\\<[0-9]+\\>"                    ; integers
    "\\<[0-9]+\\.[0-9]+\\>"           ; floats
    "\\<0x[0-9a-fA-F]+\\>")           ; hex numbers
  "List of regular expressions matching numbers to highlight."
  :type '(repeat regexp)
  :group 'highlight-numbers)

(defcustom highlight-numbers-ignore-modes
  '(calc-mode)
  "List of major modes where highlighting should be disabled."
  :type '(repeat symbol)
  :group 'highlight-numbers)

;;; Faces

(defface highlight-numbers-number
  '((t :inherit font-lock-constant-face))
  "Face for highlighted numbers."
  :group 'highlight-numbers)

;;; Variables

(defvar-local highlight-numbers--overlays nil
  "List of overlays created by highlight-numbers-mode in this buffer.")

(defvar highlight-numbers--update-timer nil
  "Timer for delayed highlighting updates.")

;;; Private Functions

(defun highlight-numbers--clear-overlays ()
  "Remove all number highlighting overlays in current buffer."
  (mapc #'delete-overlay highlight-numbers--overlays)
  (setq highlight-numbers--overlays nil))

(defun highlight-numbers--create-overlay (beg end)
  "Create highlighting overlay from BEG to END."
  (let ((ov (make-overlay beg end)))
    (overlay-put ov 'face 'highlight-numbers-number)
    (overlay-put ov 'highlight-numbers t)
    (push ov highlight-numbers--overlays)
    ov))

(defun highlight-numbers--highlight-region (start end)
  "Highlight numbers in region from START to END."
  (save-excursion
    (goto-char start)
    (dolist (pattern highlight-numbers-patterns)
      (save-excursion
        (goto-char start)
        (while (re-search-forward pattern end t)
          (let ((match-start (match-beginning 0))
                (match-end (match-end 0)))
            ;; Don't highlight if already in a string or comment
            (unless (nth 8 (syntax-ppss match-start))
              (highlight-numbers--create-overlay match-start match-end))))))))

(defun highlight-numbers--update ()
  "Update highlighting in visible portion of buffer."
  (when highlight-numbers-mode
    (let ((start (window-start))
          (end (window-end nil t)))
      ;; Clear existing overlays in visible region
      (dolist (ov highlight-numbers--overlays)
        (when (and (>= (overlay-start ov) start)
                   (<= (overlay-end ov) end))
          (delete-overlay ov)
          (setq highlight-numbers--overlays
                (delq ov highlight-numbers--overlays))))
      ;; Re-highlight visible region
      (highlight-numbers--highlight-region start end))))

(defun highlight-numbers--schedule-update ()
  "Schedule a delayed highlighting update."
  (when highlight-numbers--update-timer
    (cancel-timer highlight-numbers--update-timer))
  (setq highlight-numbers--update-timer
        (run-with-idle-timer 0.5 nil #'highlight-numbers--update)))

(defun highlight-numbers--after-change (beg end _len)
  "Hook function for `after-change-functions'.
Updates highlighting in region from BEG to END.
_LEN is ignored."
  (when highlight-numbers-mode
    (save-excursion
      ;; Extend region to cover complete lines
      (goto-char beg)
      (setq beg (line-beginning-position))
      (goto-char end)
      (setq end (line-end-position)))
    ;; Clear overlays in modified region
    (dolist (ov highlight-numbers--overlays)
      (when (and (>= (overlay-start ov) beg)
                 (<= (overlay-end ov) end))
        (delete-overlay ov)
        (setq highlight-numbers--overlays
              (delq ov highlight-numbers--overlays))))
    ;; Schedule update
    (highlight-numbers--schedule-update)))

(defun highlight-numbers--ignore-mode-p ()
  "Return non-nil if current major mode should be ignored."
  (apply #'derived-mode-p highlight-numbers-ignore-modes))

;;; Public Functions

;;;###autoload
(defun highlight-numbers-refresh ()
  "Refresh number highlighting in current buffer."
  (interactive)
  (when highlight-numbers-mode
    (highlight-numbers--clear-overlays)
    (highlight-numbers--highlight-region (point-min) (point-max))
    (message "Number highlighting refreshed")))

;;; Mode Definition

(defvar highlight-numbers-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map (kbd "C-c h r") #'highlight-numbers-refresh)
    map)
  "Keymap for `highlight-numbers-mode'.")

;;;###autoload
(define-minor-mode highlight-numbers-mode
  "Toggle highlighting of numbers in current buffer.

When enabled, all numbers matching patterns in
`highlight-numbers-patterns' are highlighted with
`highlight-numbers-number' face.

\\{highlight-numbers-mode-map}"
  :lighter " #"
  :keymap highlight-numbers-mode-map
  :group 'highlight-numbers
  (if highlight-numbers-mode
      (highlight-numbers--enable)
    (highlight-numbers--disable)))

(defun highlight-numbers--enable ()
  "Enable number highlighting in current buffer."
  (cond
   ((highlight-numbers--ignore-mode-p)
    (setq highlight-numbers-mode nil)
    (message "Number highlighting not available in %s" major-mode))
   (t
    ;; Add hooks
    (add-hook 'after-change-functions
              #'highlight-numbers--after-change nil t)
    ;; Initial highlighting
    (highlight-numbers--highlight-region (point-min) (point-max))
    (message "Number highlighting enabled"))))

(defun highlight-numbers--disable ()
  "Disable number highlighting in current buffer."
  ;; Remove hooks
  (remove-hook 'after-change-functions
               #'highlight-numbers--after-change t)
  ;; Clear overlays
  (highlight-numbers--clear-overlays)
  ;; Cancel pending timer
  (when highlight-numbers--update-timer
    (cancel-timer highlight-numbers--update-timer)
    (setq highlight-numbers--update-timer nil))
  (message "Number highlighting disabled"))

;;; Global Minor Mode

;;;###autoload
(define-globalized-minor-mode global-highlight-numbers-mode
  highlight-numbers-mode
  highlight-numbers--turn-on
  :group 'highlight-numbers)

(defun highlight-numbers--turn-on ()
  "Turn on `highlight-numbers-mode' if appropriate."
  (unless (or (minibufferp)
              (highlight-numbers--ignore-mode-p))
    (highlight-numbers-mode 1)))

(provide 'highlight-numbers)
;;; highlight-numbers.el ends here
