FROM python:3.8-slim
WORKDIR /app

COPY ./requirements.txt /app/requirements.txt

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    bash

RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -U memory_profiler
RUN pip install --no-cache-dir -r requirements.txt

COPY ./run_metadata_with_different_datasets.sh /app/run_metadata_with_different_datasets.sh


