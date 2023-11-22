.PHONY: all clean pytest coverage flake8 black mypy isort

CMD:=poetry run
PYMODULE:=src
TESTS:=tests

# Run all the checks which do not change files
all: isort black flake8 pytest

# Run the unit tests using `pytest`
pytest:
	$(CMD) pytest $(PYMODULE) $(TESTS)

# Lint the code using `flake8`
flake8:
	$(CMD) flake8 $(PYMODULE) $(TESTS)

# Automatically format the code using `black`
black:
	$(CMD) black $(PYMODULE) $(TESTS)

# Order the imports using `isort`
isort:
	$(CMD) isort $(PYMODULE) $(TESTS)

# Serve docs
serve:
	$(CMD) mkdocs serve

deploy:
	$(CMD) mkdocs gh-deploy --force
