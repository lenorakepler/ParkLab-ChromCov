# Pinned, reproducible environment for the chromcov coverage step.
# Two stages: resolve deps from the committed uv.lock, then a slim runtime.

FROM python:3.13-slim AS build

# uv installs from the lockfile, so the dependency set is byte-for-byte pinned.
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

WORKDIR /app
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy

# Install deps first (cached layer), then the project itself.
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-install-project --no-dev
COPY chromcov ./chromcov
COPY main.py README.md ./
RUN uv sync --frozen --no-dev

FROM python:3.13-slim AS runtime

# htslib CLI tools are handy for debugging CRAMs even though pysam bundles htslib.
RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=build /app /app
ENV PATH="/app/.venv/bin:$PATH"

# `chromcov` is the console entrypoint the CWL tool calls (define it in
# pyproject.toml [project.scripts]).
ENTRYPOINT ["chromcov"]
