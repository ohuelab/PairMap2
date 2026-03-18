FROM python:3.12-slim

WORKDIR /app

# System dependencies for RDKit build
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY pyproject.toml ./
COPY local/ ./local/
COPY pairmap2/ ./pairmap2/
COPY web/ ./web/

# Install dependencies
RUN pip install --no-cache-dir \
    "git+https://github.com/OpenFreeEnergy/Lomap@v3.0.1" \
    "git+https://github.com/ohuelab/PairMap.git" \
    "networkx>=3.0" "numpy>=1.24" "rdkit>=2023.3" "tqdm>=4.65" "pandas>=2.0" \
    "fastapi>=0.100" "uvicorn[standard]>=0.20" "python-multipart>=0.0.5" \
    "./local/gufe_stub" \
    -e .

# Writable jobs/cache directories
RUN mkdir -p /app/jobs /app/cache

EXPOSE 8000

CMD ["uvicorn", "web.backend.main:app", \
     "--host", "0.0.0.0", "--port", "8000"]
