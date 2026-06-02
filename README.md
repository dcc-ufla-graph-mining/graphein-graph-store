--- 
In this work, we introduce PDB Graph Store, a high-productivity, integrated pro-
gramming library designed to optimize memory usage in protein graph modeling work-
flows. The proposed approach reorganizes protein graph data into compact, shared repre-
sentations that significantly reduce space requirements. Implemented as a Python pack-
age, PDB Graph Store integrates easily with existing scientific computing and graph
analysis tools, facilitating its adoption in research and experimental settings
---


--- para executar algum experimento, descomente-o em src/compress/main.py e entao set a variavel DATASET para o nome do dataset desejado em docker-compose.yml, ou descomente a linha `command: bash run_metadata_with_different_datasets.sh` e comente a linha `command: python src/compress/main.py` para executar com todos os datasets em ./data/

execute com o comando `sh build`
---

--- 

`sh lab`
execute lab/init para iniciar um container com uma sessão jupyter para interagir com a ferramenta via notebook python.
---


---
estrutura:

data/  -- diretorio com os datasets
errors/ -- diretorio para logar os erros de execucao
jupyter-lab/ -- diretorio com docker compose para iniciar o jupyter lab
memory_footprint_results/ -- diretorio com os graficos mostrando o uso de memoria durante a execucao da aplicacao usando memory profiler
min_max_results/ -- diretorio contendo dados descritivos sobre a estrutura dos grafos (tamanho min, max, avg e total de nodes e arestas por dataset)
results/ -- resultados dos experimentos 1-6 executados
src/ -- codigo fonte da ferramenta + codigo para geraçao dos experimentos
times/ -- diretorio com os dados descritivos mostrando o uso de memoria durante a execucao da aplicacao usando memory profiler
build -- script para rodar a aplicaçao em modo experimento (src/main.py)
lab -- script para rodar jupyter lab para teste da aplicacao

---