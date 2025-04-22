FROM python:3.10-slim
WORKDIR /app
COPY ./requirements.txt /app
COPY ./src /app/src
COPY ./data /app/data

RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make 

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

ENV DATA_DIR=/app/data

CMD ["python", "src/compress/metadata.py"]