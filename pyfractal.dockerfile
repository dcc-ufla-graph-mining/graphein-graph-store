FROM rayproject/ray:nightly-py310-cpu

USER root

ARG PYTHON_VERSION=3.10.12
ENV RAY_HOME=/home/ray

RUN apt-get update && apt-get install -y \
    wget \
    tar \
    gcc \
    g++ \
    make 


RUN $RAY_HOME/anaconda3/bin/conda init && \
    echo 'export PATH=$RAY_HOME/anaconda3/bin:$PATH' >> /home/ray/.bashrc && \
    $RAY_HOME/anaconda3/bin/conda install -y libgcc-ng python=$PYTHON_VERSION && \
    $RAY_HOME/anaconda3/bin/conda install -y -c conda-forge libffi=3.4.2 && \
    $RAY_HOME/anaconda3/bin/conda clean -y --all 

COPY ./requirements-all.txt /tmp/requirements-all.txt

RUN $RAY_HOME/anaconda3/bin/pip install --upgrade pip

RUN $RAY_HOME/anaconda3/bin/pip install -r /tmp/requirements-all.txt

ENV FRACTAL_HOME=/app/fractal
ENV SPARK_HOME=$FRACTAL_HOME/python/sparkbuild/spark


WORKDIR /app

EXPOSE 8265

RUN mkdir -p $SPARK_HOME && mkdir -p $FRACTAL_HOME

RUN rm -rf /var/lib/apt/lists/* \
    && wget https://download.java.net/java/ga/jdk11/openjdk-11_linux-x64_bin.tar.gz \
    && tar -xvf openjdk-11_linux-x64_bin.tar.gz \
    && mkdir /usr/lib/jvm/ \
    && mv jdk-11 /usr/lib/jvm/java-11-openjdk \
    && update-alternatives --install /usr/bin/java java /usr/lib/jvm/java-11-openjdk/bin/java 1 \
    && rm -rf openjdk-11_linux-x64_bin.tar.gz

ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk

COPY ./fractal $FRACTAL_HOME
COPY ./spark $SPARK_HOME
COPY ./requirements.txt $FRACTAL_HOME/requirements.txt

RUN cd $FRACTAL_HOME && chmod +x ./gradlew && ./gradlew jar 

RUN $RAY_HOME/anaconda3/bin/pip install --no-cache-dir fractal/python/sparkbuild/spark/python 2>&1 > installation.log
RUN $RAY_HOME/anaconda3/bin/pip install --no-cache-dir fractal/python 2>&1 >> installation.log
RUN $RAY_HOME/anaconda3/bin/pip install --no-cache-dir -r $FRACTAL_HOME/requirements.txt