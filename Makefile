.PHONY:
default: help

.PHONY: help
help:
	@echo help      -- display this help
	@echo black     -- apply black code formatter
	@echo snakefmt  -- apply snakefmt code formatter
	@echo srcfmt    -- apply black and snakefmt formatters
	@echo lint      -- run linters
	@echo test      -- run tests through pytest
	@echo isort     -- run isort

.PHONY: black
black:
	black -l 100 .

.PHONY: snakefmt
snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)'

.PHONY: srcfmt
srcfmt: black snakefmt

.PHONY: lint
lint: flake8

.PHONY: isort
isort:
	isort --force-sort-within-sections --profile=black .

.PHONY: flake8
flake8:
	flake8

test:
	py.test

coverage:
	coverage report
	coverage html
