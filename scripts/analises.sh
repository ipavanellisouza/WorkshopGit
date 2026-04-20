##Montagem do genoma de Sucupira usando dados Nanopore

##17/03/2026

#Corrida Nanopore Ppu
#Nome da corrida: Cbr-Ppu-26-01-2026
#n50 ~8 kb
cp /media/lgbio-nas1/sequenciamentos/PromethION/Cbr-Ppu-26-01-2026/Ppu/20260126_1712_P2S-01745-B_PBA04053_0bbb3224/fastq_pass/PBA04053_pass_0bbb3224_8c2adb75_* .
zcat PBA04053_pass_0bbb3224_8c2adb75_* | pigz > Ppu_26_01_26_merged.fastq.gz

#Corrida Nanopore Ppu
#Nome da corrida: Ppu-Edy-29-01-2026
#N50 - ~6 kb
cp /media/lgbio-nas1/sequenciamentos/PromethION/Ppu-Edy-29-01-2026/Ppu/20260129_1641_P2S-01745-A_PBA04053_6c0ebec5/fastq_pass/PBA04053_pass_6c0ebec5_7e006c3f_* .
zcat PBA04053_pass_6c0ebec5_7e006c3f_* | pigz > Ppu_29_01_26_merged.fastq.gz

#Corrida Nanopore Ppu
#Nome da corrida: Ppu-09-03-2026
#N50 - 23 kb
cp /media/lgbio-nas1/sequenciamentos/PromethION/Ppu09-03-2026/Ppu/20260309_1617_P2S-01745-B_PBA04053_ade413f1/fastq_pass/PBA04053_pass_ade413f1_d444857b_* .
zcat PBA04053_pass_ade413f1_d444857b_* | pigz > Ppu_09_03_26_merged.fastq.gz

#Corrida Nanopore Ppu
#Nome da corrida: Ppu-Sol-10-03-2026
#N50 - ~23 kb
cp  /media/lgbio-nas1/sequenciamentos/PromethION/Ppu-Sol-10-03-2026/Ppu/20260310_1643_P2S-01745-A_PBC92083_5cfd19d1/fastq_pass/PBC92083_pass_5cfd19d1_167fa0b6_* .
zcat PBC92083_pass_5cfd19d1_167fa0b6_* | pigz > Ppu_10_03_26_merged.fastq.gz
 
#Corrida Illumina Miseq
#Nome da corrida: Pterodon (MA1)
#/media/lgbio-nas1/sequenciamentos/Backup_Drive/Pterodon-genomico

##FastQC
fastqc -t 32 Ppu_* -o 1.QC_Dadosbrutos/
NanoPlot --fastq ./../Ppu_*  -o ./ 
multiqc 1.QC_dadosbrutos -n multiqc_fastqc_report -o 1.QC_dadosbrutos/

#Filtragem
porechop_abi -abi -i Ppu_26_01_26_merged.fastq.gz -o 3.Filtragem_dadosbrutos/Ppu_26_01_26_trim.fq -t 20
##Saved result to /media/lgbio-nas1/isabela.pavanelli/ONT_Sucupira/3.Filtragem_dadosbrutos/Ppu_26_01_26_trim.fq

porechop_abi -abi -i Ppu_29_01_26_merged.fastq.gz -o 3.Filtragem_dadosbrutos/Ppu_29_01_26_trim.fq -t 20
##Saved result to /media/lgbio-nas1/isabela.pavanelli/ONT_Sucupira/3.Filtragem_dadosbrutos/Ppu_29_01_26_trim.fq

porechop_abi -abi -i Ppu_09_03_26_merged.fastq.gz -o 3.Filtragem_dadosbrutos/Ppu_09_03_26_trim.fq -t 20
##Saved result to /media/lgbio-nas1/isabela.pavanelli/ONT_Sucupira/3.Filtragem_dadosbrutos/Ppu_09_03_26_trim.fq

porechop_abi -abi -i Ppu_10_03_26_merged.fastq.gz -o 3.Filtragem_dadosbrutos/Ppu_10_03_26_trim.fq -t 20
#Saved result to /media/lgbio-nas1/isabela.pavanelli/ONT_Sucupira/3.Filtragem_dadosbrutos/Ppu_10_03_26_trim.fq

#QC_dados filtrados
4.QC_dados_filtrados$ NanoPlot --fastq ../3.Filtragem_dadosbrutos/Ppu_* -o . --threads 20

meryl k=31 memory=64 threads=12 count ../3.Filtragem_dadosbrutos/Ppu_* output Ppu_trim.meryl

meryl histogram Ppu_trim.meryl | sed 's/\t/ /g' > ppu_trim.meryl.hist

##Montagem
mkdir 5.Montagem
cd 5.Montagem
hifiasm -o ONT_Sucupira -t 24 --primary --write-ec --ont ../3.Filtragem_dadosbrutos/Ppu_* 2> Sucupira_hifiasm.log

awk '/^S/{print ">"$2;print $3}' ONT_Sucupira.p_ctg.gfa > ONT_Sucupira.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ONT_Sucupira.a_ctg.gfa > ONT_Sucupira.a_ctg.fasta

quast.py 5.Montagem/ONT_Sucupira.p_ctg.fasta -o output_quast_Sucupira

##Ficou lindissimo L90 de 8 - numero de cromossomos da especie haploide

## Fazer montagem de telomeros
mkdir 6.Montagem-telomero
cd 6.Montagem-telomero
/ONT_Sucupira/6.Montagem-telomero$ hifiasm -o sucupira_telomero -t 20 --telo-m TTAGGG ../3.Filtragem_dadosbrutos/Ppu_* 2> sucupira_telo.log

#transformando em fasta
awk '/^S/{print ">"$2;print $3}' sucupira_telomero.bp.hap1.p_ctg.gfa > sucupira_telomero.bp.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' sucupira_telomero.bp.hap2.p_ctg.gfa > sucupira_telomero.bp.hap2.p_ctg.fasta


pigz -p 32 sucupira_telomero.bp.hap1.p_ctg.fasta

##BUSCO

busco -i ONT_Sucupira.p_ctg.fasta -l fabales_odb12 -o busco_sucupira -m genome -c 2
