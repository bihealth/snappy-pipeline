.PHONY:
default: help

.PHONY: help
help:
	@echo help      -- display this help
	@echo rufffmt   -- apply ruff code formatter
	@echo black     -- apply black code formatter
	@echo snakefmt  -- apply snakefmt code formatter
	@echo srcfmt    -- apply black and snakefmt formatters
	@echo lint      -- run linters
	@echo test      -- run tests through pytest

.PHONY: rufffmt
rufffmt:
	ruff format .

.PHONY: black
black:
	black -l 100 .


.PHONY: snakefmt
snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)'

.PHONY: srcfmt
srcfmt: black rufffmt snakefmt

.PHONY: lint
lint: lint-ruff lint-snakefmt

.PHONY: lint-ruff
lint-ruff:
	ruff check .

.PHONY: lint-snakefmt
lint-snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)' --check

test:
	py.test

coverage:
	coverage report
	coverage html
