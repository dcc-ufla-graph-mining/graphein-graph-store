

--- para executar algum experimento, descomente-o em src/compress/main.py e entao set a variavel DATASET para o nome do dataset desejado em docker-compose.yml, ou descomente a linha `command: bash run_metadata_with_different_datasets.sh` e comente a linha `command: python src/compress/main.py` para executar com todos os datasets em ./data/

execute com o comando `sh build`
---

--- 

`sh lab`
execute lab/init para iniciar um container com uma sessão jupyter para interagir com a ferramenta via notebook python.
---