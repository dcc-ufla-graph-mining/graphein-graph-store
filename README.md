# PDB Graph Store

PDB Graph Store is a high-productivity Python library for protein graph compression and storage optimization. The project focuses on reducing memory usage in protein graph modeling workflows by reorganizing graph structures into compact and shared representations.

The library integrates with common scientific computing and graph analysis ecosystems, allowing compressed graph representations to be used with minimal friction in research pipelines, experimental workflows, and large-scale graph processing tasks.

The compression workflow used by the tool is illustrated below:

<!-- workflow image here -->

---

## Overview

Protein graph datasets often contain large amounts of structural redundancy, leading to excessive memory consumption during graph loading, manipulation, and analysis.

PDB Graph Store addresses this problem by:

* Compressing protein graph representations
* Sharing repeated structures across graphs
* Reducing memory footprint during execution
* Providing reusable graph storage abstractions
* Supporting large-scale protein graph experimentation

The library is implemented in Python and designed to integrate naturally with existing graph-processing workflows.

---

# Repository Structure

```text
data/
├── Protein graph datasets

errors/
├── Execution error logs

jupyter-lab/
├── Docker Compose configuration for Jupyter Lab sessions

memory_footprint_results/
├── Memory usage plots generated with memory_profiler

min_max_results/
├── Descriptive statistics about graph structures (min/max/avg/total nodes and edges per dataset)

results/
├── Results generated from experiments 1–6

src/
├── Source code of the library
├── Experiment execution scripts

times/
├── Execution time measurements and profiling data

build
├── Script to execute experiments

lab
├── Script to launch Jupyter Lab
```

---

# Running Experiments

The project uses Docker Compose for reproducible execution.

## Running a Single Experiment

To execute a specific experiment:

1. Uncomment the desired experiment inside:

```text
src/compress/main.py
```

2. Set the dataset name in:

```text
./docker-compose.yaml
```

using the `DATASET` environment variable.

3. Execute:

```bash
sh build
```

---

## Running Experiments for All Datasets

To execute the experiment pipeline for all datasets inside `./data/`:

1. In `docker-compose.yml`:

* Uncomment:

```yaml
command: bash run_metadata_with_different_datasets.sh
```

* Comment:

```yaml
command: python src/compress/main.py
```

2. Execute:

```bash
sh build
```

---

# Jupyter Lab Environment

An interactive Jupyter environment is available for exploratory usage and notebook-based experimentation.

Start the environment with:

```bash
sh lab
```

After the container starts, access Jupyter Lab in your browser at:

```text
http://127.0.0.1:8888/lab/tree/src
```

The environment is already configured to interact with the PDB Graph Store library and datasets.

A notebook named lab.ipynb is included with the project and contains example workflows demonstrating how to use the compression and graph manipulation features provided by the library.

---

# Features

* Protein graph compression
* Shared graph representations
* Reduced memory consumption
* Experimental benchmarking support
* Memory profiling integration
* Docker-based reproducible execution
* Jupyter-based interactive experimentation

---

# Requirements

* Docker
* Docker Compose
* Bash shell

---

# Intended Usage

PDB Graph Store is intended for:

* Bioinformatics experimentation
* Graph compression research
* Large-scale graph processing workflows
* Memory-efficient graph storage studies

---
