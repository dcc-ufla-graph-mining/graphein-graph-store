FROM python:3.8-slim
WORKDIR /app

COPY ./requirements.txt /app/requirements.txt

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make 

RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -U memory_profiler
RUN pip install --no-cache-dir -r requirements.txt


