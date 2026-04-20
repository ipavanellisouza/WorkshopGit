# Análises de Montagem do Genoma de Sucupira

Este script (`analises.sh`) é um pipeline automatizado para a montagem do genoma da espécie Sucupira utilizando dados de sequenciamento Nanopore. Ele inclui etapas de cópia de dados, fusão de arquivos FASTQ, controle de qualidade (QC), filtragem de adaptadores, contagem de k-mers, montagem do genoma com hifiasm, avaliação da montagem e análise de completude com BUSCO.

## Data
17/03/2026

## Descrição das Etapas

1. **Cópia e Fusão de Dados Nanopore**:
   - Copia arquivos FASTQ de corridas específicas do Nanopore.
   - Funde os arquivos em um único arquivo comprimido por corrida.

2. **Controle de Qualidade Inicial**:
   - Executa FastQC e NanoPlot nos dados brutos.
   - Gera relatórios MultiQC.

3. **Filtragem de Dados**:
   - Remove adaptadores usando `porechop_abi`.

4. **Controle de Qualidade Pós-Filtragem**:
   - Executa NanoPlot nos dados filtrados.
   - Conta k-mers com `meryl` para análise de frequência.

5. **Montagem do Genoma**:
   - Usa `hifiasm` para montar o genoma a partir dos dados filtrados.
   - Converte os grafos de montagem em arquivos FASTA.

6. **Avaliação da Montagem**:
   - Executa QUAST para métricas de montagem.
   - Executa BUSCO para avaliação de completude baseada em genes conservados.

7. **Montagem com Telômeros**:
   - Realiza uma montagem adicional considerando sequências teloméricas (TTAGGG).

## Pré-requisitos

Antes de executar o script, certifique-se de que as seguintes ferramentas estão instaladas e disponíveis no PATH:

- `pigz` (para compressão/descompressão)
- `fastqc` (para controle de qualidade)
- `NanoPlot` (para visualização de dados Nanopore)
- `multiqc` (para relatórios agregados)
- `porechop_abi` (para remoção de adaptadores)
- `meryl` (para contagem de k-mers)
- `hifiasm` (para montagem de genoma)
- `quast.py` (para avaliação de montagem)
- `busco` (para análise de completude)

Além disso, os caminhos absolutos no script (como `/media/lgbio-nas1/...`) devem ser ajustados para o seu ambiente de sistema.

## Como Usar

1. Ajuste os caminhos no script para corresponder à sua estrutura de diretórios.
2. Torne o script executável: `chmod +x analises.sh`
3. Execute o script: `./analises.sh`

**Nota**: O script assume que você tem acesso aos diretórios de origem dos dados. Se necessário, monte os volumes ou copie os arquivos manualmente.

## Saídas Esperadas

- Arquivos FASTQ fundidos: `Ppu_*.fastq.gz`
- Relatórios de QC: `1.QC_Dadosbrutos/`, `4.QC_dados_filtrados/`
- Dados filtrados: `3.Filtragem_dadosbrutos/`
- Montagem: `5.Montagem/ONT_Sucupira.*`
- Montagem com telômeros: `6.Montagem-telomero/`
- Avaliações: `output_quast_Sucupira/`, `busco_sucupira/`

## Observações

- O script foi projetado para um ambiente específico e pode requerer modificações para outros sistemas.
- Verifique os logs de saída para erros durante a execução.
- A montagem resultou em um L90 de 8, indicando um bom número de cromossomos para a espécie haploide.

Para mais detalhes, consulte os comentários no script ou os logs gerados.