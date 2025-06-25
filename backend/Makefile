.PHONY: setup lint format clean

setup:
	python -m venv venv
	. venv/bin/activate && pip install -r requirements.txt

lint:
	ruff check .

format:
	ruff check --fix-only
	ruff check --select I --fix
	ruff format

clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} \;
	find . -name ".ruff_cache" -type d -prune -exec rm -rf {} \;
