FROM jupyter/base-notebook

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    make \
    gcc \
    g++ \
    bash \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy and install python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt

USER ${NB_UID}
