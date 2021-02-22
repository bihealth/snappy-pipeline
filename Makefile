.PHONY:
default: help

.PHONY: help
help:
	@echo help      -- display this help
	@echo black     -- apply black code formatter
	@echo snakefmt  -- apply snakefmt code formatter
	@echo lint      -- run linters
	@echo test      -- run tests through pytest

.PHONY: black
black:
	black -l 100 .

.PHONY: snakefmt
snakefmt:
	snakefmt -l 100 . --include '(\.smk$$|\.rules$$|^Snakefile)'

.PHONY: lint
lint: prospector

.PHONY: prospector
prospector:

test:
	py.test

coverage:
	coverage report
	coverage html
