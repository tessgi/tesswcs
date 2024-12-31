.PHONY: all clean pytest coverage flake8 black mypy isort

# Run all the checks which do not change files
all: isort black flake8 pytest

# Run the unit tests using `pytest`
pytest:
	poetry run pytest src tests

# Lint the code using `flake8`
flake8:
	poetry run flake8 src tests

# Automatically format the code using `black`
black:
	poetry run black src tests

# Order the imports using `isort`
isort:
	poetry run isort src tests

# Serve docs
serve:
	poetry run mkdocs serve

deploy:
	poetry run mkdocs gh-deploy --force
