# sampling.analysis
## Função para análise do processo de amostragem em levantamentos fitossociológicos em função do número de indivíduos e da área basal.
#### Autor:  Pedro Higuchi                                   
 26/05/2019

* Como citar

* Higuchi, P. sampling.analysis: Função em linguagem de programação estatística R para análise do processo amostragem de levantamento fitossocioógicos em função do número de indivíduos e área basal. 2019. Disponvel em https://github.com/higuchip/sampling.analysis

# Observações:											                      
- a) O argumento x (planilha de dados) terá que conter as colunas parc (identificação das parcelas), spp (id. espécies),dap e estratos (obrigatório apenas no caso de amostragem estratificada);
- b) arquivo exemplo de entrada, disponível em https://raw.githubusercontent.com/higuchip/forest.din/master/dados_exemplo_amostragem.csv
- c) O argumento sys, representa o sistema de amostragem (AS = Amostragem Simples; EST = Estratificada e SIS = Sistemática);
- d) O argumento plot_size representa o tamanho de cada parcela em m2;
- e) O argumento forest_area representa o tamanho da área florestal em ha;
- f) O argumento strata_area deve ser usado apenas quando sys = EST e representa um vetor numérico com as áreas de cada estrato em ha;
- g) o argumento alfa representa o nível de significância pela tabela t de Student;
- h) o argumento LE representa o Limite de erro admissível.


Caso tenha dúvidas, sugestões ou queira contribuir, entre em contato: higuchip@gmail.com


# Referência:
-FELFILI, JM et al. Procedimentos e métodos de amostragem de vegetação. Fitossociologia no Brasil: métodos e estudos de casos. Viçosa, MG: Universidade Federal de Viçosa, p. 86-121, 2011.
