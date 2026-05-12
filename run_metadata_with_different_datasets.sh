#!/bin/bash

if $(ls /app | grep -q "logs"); then
    echo "Directory /app/logs exists."
else
    mkdir /app/logs
fi

if $(ls /app/logs | grep -q "metadata.log"); then
    rm /app/logs/metadata.log
fi

touch /app/logs/metadata.log

# Check if the error directory exists
if $(ls /app | grep -q "errors"); then
    echo "Directory /app/errors exists."
else
    mkdir /app/errors
fi

# Check if the results directory exists
if $(ls /app | grep -q "results"); then
    echo "Directory /app/results exists."
else
    mkdir /app/results
fi

python -m memory_profiler src/compress/trash_measure_memory.py > "times/trash.txt"

EXCOUNT=0

# for dataset in data/*.txt
# do
    # for i in {1..5};
    # do
        # dataset_name=$(basename $dataset)
        # export DATASET=$dataset_name
        dataset_name="ligand_PLP.txt"
        export DATASET=$dataset_name
        export EXCOUNT
        echo $dataset_name
        python src/compress/main.py >> /app/logs/metadata.log 2>&1
        # mprof run --python python <script>
        mprof run --python python src/compress/measure_memory_v1.py > "times/${dataset_name}_v1.txt"
        mprof plot --output memory_usage_v1.pdf
        mprof run --python python src/compress/measure_memory_v2.py > "times/${dataset_name}_v2.txt"
        mprof plot --output memory_usage_v2.pdf

        EXCOUNT=$i
    # done
# done
