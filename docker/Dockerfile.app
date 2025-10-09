FROM python:3.11-slim
RUN apt-get update && apt-get install -y build-essential git curl && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir poetry
WORKDIR /workspace
COPY pyproject.toml /workspace/
RUN poetry config virtualenvs.create false && poetry install --no-interaction --no-ansi --all-extras
COPY . /workspace