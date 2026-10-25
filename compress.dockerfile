FROM python:3.10-slim
WORKDIR /app

COPY ./requirements.txt /app/requirements.txt

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    bash \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip "setuptools<82" wheel

RUN echo "setuptools<82" > /etc/pip-build-constraint.txt
ENV PIP_BUILD_CONSTRAINT=/etc/pip-build-constraint.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY ./run_metadata_with_different_datasets.sh /app/run_metadata_with_different_datasets.sh


