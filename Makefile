.PHONY: help install install-dev test coverage lint format clean validate-skills

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'

install:  ## Install package
	pip install -e .

install-dev:  ## Install package with dev dependencies
	pip install -e ".[dev]"

test:  ## Run tests
	pytest tests/ -v

coverage:  ## Run tests with coverage report
	pytest tests/ --cov=cli --cov-report=html --cov-report=term

lint:  ## Run linting checks
	@echo "Running Ruff..."
	ruff check .
	@echo "Checking with Black..."
	black --check .

format:  ## Auto-format code
	@echo "Formatting with Black..."
	black .
	@echo "Running Ruff fixes..."
	ruff check --fix .

clean:  ## Clean up build artifacts and cache
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

validate-skills:  ## Validate all skills and commands
	@echo "Validating skills..."
	@for skill_dir in skills/*/*/; do \
		if [ -f "$$skill_dir/SKILL.md" ]; then \
			echo "Checking $$skill_dir"; \
			python -m cli.main validate "$$skill_dir" || exit 1; \
		fi \
	done
	@echo "Validating commands..."
	@for cmd_file in commands/**/*.md; do \
		if [ -f "$$cmd_file" ]; then \
			echo "Checking $$cmd_file"; \
			python -m cli.main validate "$$cmd_file" || exit 1; \
		fi \
	done
	@echo "✓ All skills and commands are valid"

build:  ## Build distribution packages
	python -m build

publish:  ## Publish to PyPI (requires credentials)
	python -m twine upload dist/*

publish-test:  ## Publish to TestPyPI
	python -m twine upload --repository testpypi dist/*

list:  ## List available skills and commands
	python -m cli.main list

info:  ## Show info about a skill (usage: make info SKILL=python-ase)
	python -m cli.main info $(SKILL)

install-skill:  ## Install a skill (usage: make install-skill SKILL=python-ase)
	python -m cli.main install $(SKILL)

create-skill:  ## Create a new skill interactively
	python -m cli.main create --type skill

create-command:  ## Create a new command interactively
	python -m cli.main create --type command
