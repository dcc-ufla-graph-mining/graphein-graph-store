if $(ls /app | grep -q "logs"); then
    echo "Directory /app/logs exists."
else
    mkdir /app/logs
fi

if $(ls /app/logs | grep -q "metadata.log"); then
    echo "File /app/logs/metadata.log exists."
else
    touch /app/logs/metadata.log
fi

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

for dataset in data/*.txt
do
    dataset_name=$(basename $dataset)
    export DATASET=$dataset_name
    echo $dataset_name
    python src/compress/main.py >> /app/logs/metadata.log 2>&1
done