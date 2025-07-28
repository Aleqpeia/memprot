FROM continuumio/miniconda3:latest

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libfftw3-dev \
    libopenmpi-dev \
    openmpi-bin \
    git \
    vim \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
RUN useradd -m -u 1000 mduser && \
    chown -R mduser:mduser /opt/conda

USER mduser
WORKDIR /home/mduser

# Install conda dependencies first (faster layer caching)
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    conda clean -afy

# Activate environment
SHELL ["conda", "run", "-n", "memprot", "/bin/bash", "-c"]

# Copy project files
COPY --chown=mduser:mduser . /home/mduser/memprot-analysis

# Install the package
WORKDIR /home/mduser/memprot-analysis
RUN conda run -n memprot pip install -e .

# Set the entry point
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "memprot"]
CMD ["memprot-analyze", "--help"]

# Expose port for Jupyter if needed
EXPOSE 8888