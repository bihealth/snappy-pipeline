.PHONY:
default: help

.PHONY: help
help:
	@echo help      -- display this help
	@echo fmt       -- apply code formatter
	@echo snakefmt  -- apply snakefmt code formatter
	@echo srcfmt    -- apply black and snakefmt formatters
	@echo lint      -- run linters
	@echo test      -- run tests through pytest
	@echo check     -- run code checker

.PHONY: fmt
fmt:
	ruff format .

.PHONY: snakefmt
snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)'

.PHONY: srcfmt
srcfmt: fmt snakefmt

.PHONY: lint
lint: check lint-snakefmt

.PHONY: check
check:
	ruff check .

.PHONY: lint-fmt
lint-fmt:
	ruff format --check .

.PHONY: lint-snakefmt
lint-snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)' --check


test:
	pytest

coverage:
	coverage report
	coverage html
