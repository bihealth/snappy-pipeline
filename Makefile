.PHONY:
default: help

.PHONY: help
help:
	@echo help      -- display this help
	@echo black     -- apply ruff code formatter
	@echo snakefmt  -- apply snakefmt code formatter
	@echo srcfmt    -- apply black and snakefmt formatters
	@echo lint      -- run linters
	@echo test      -- run tests through pytest

.PHONY: black
black:
	ruff .

.PHONY: snakefmt
snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)'

.PHONY: srcfmt
srcfmt: black snakefmt

.PHONY: lint
lint: lint-ruff lint-snakefmt

.PHONY: lint-ruff
lint-ruff:
	ruff check --fix .

.PHONY: lint-snakefmt
lint-snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)' --check

test:
	py.test

coverage:
	coverage report
	coverage html
