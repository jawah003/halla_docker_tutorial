# HAllA CLI runner image for the halla-docker-tutorial.


FROM continuumio/miniconda3:latest
# NOTE: continuumio/miniconda3 is deprecated upstream (its successor is
# `anaconda/miniconda`) and prints a deprecation banner on every run. It is kept
# here because this is the base that the recipe below is verified against. If you
# switch bases, re-verify the conda solve in the next step.

WORKDIR /app

# Build toolchain plus the shared libraries Kaleido needs to render figures in a
# container with no display.
RUN apt-get update && apt-get install -y \
    build-essential \
    libnss3 \
    libfontconfig1 \
    libdbus-1-3 \
    libx11-xcb1 \
    libasound2 \
    libxtst6 \
    libxss1 \
    libegl1 \
    libgl1 \
    libglib2.0-0 \
    libxrender1 \
    libxi6 \
    libxrandr2 \
    && rm -rf /var/lib/apt/lists/*


RUN conda create -n halla -c conda-forge -y \
    python=3.9 \
    r-base=4.3 \
    r-eva \
    r-xicor \
    rpy2 \
    python-kaleido
ENV PATH=/opt/conda/envs/halla/bin:$PATH

# HAllA is pinned so this tutorial keeps producing the same numbers.
#
# Pillow is pinned below 10 on purpose: Pillow 10 removed the deprecated `mode`
# argument that the matplotlib version in HAllA's dependency set still passes,
# which crashes PDF hallagram export at the very end of an otherwise successful
# run.
RUN pip install --no-cache-dir \
    halla==0.8.40 \
    "Pillow<10.0.0"

CMD ["halla", "--help"]