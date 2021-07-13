# Genome annotation for *Penstemon barbatus* 

## Contents
1. Introduction
2. Software links & citations
3. Annotation

## 1. Introduction
### Plant & Genome 
  *Penstemon barbatus* is a plant species found within Plantaginaceae. This species, like most within *Penstemon* is diploid with 8 pairs of chromosomes (Freeman, 1983). 
  
### References
1. Freeman, C.C. Chromosome numbers in Great Plains species of Penstemon (Scrophulariaceae). *Brittonia* 35, 232–238 (1983). https://doi.org/10.2307/2806022

## 2. Software links
### Used in multiple places
1. BUSCO/3.0.2
   - https://busco.ezlab.org/
   - Manni M, Berkeley MR, Seppey M, Simao FA, Zdobnov EM. 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. arXiv:2106.11799 [q-bio] [Internet]. 

### Transcriptome Assembly
1. Fastp/0.14.1
   - https://github.com/OpenGene/fastp
   - Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
3. Trinity/2.8.3
   - https://github.com/trinityrnaseq/trinityrnaseq/wiki
   - Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.
5. Bowtie2/2.3.4.3
   - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
   - Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

### Genome Annotation
1. RepeatModeler/2.0.1
   - http://www.repeatmasker.org/RepeatModeler/
   - Smit, AFA, Hubley, R. RepeatModeler Open-1.0. 2008-2015 <http://www.repeatmasker.org>. 
3. MAKER/2.31.10
   - http://www.yandell-lab.org/software/maker.html
   - Holt, C., Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491 (2011). https://doi.org/10.1186/1471-2105-12-491
5. MPICH/3.3.2
   - https://www.mpich.org/

## 3. Genome Annotation
### Transcriptome assembly
  We used three tissue types to assemble the transcriptome: floral (separated into lateral stamens, ventral stamens, and rest of flower), root, and shoot. Each sample was sequenced using NovaSeq 150 bp paired-end sequencing at KU Medical Center Genomics Core (https://www.kumc.edu/genomics.html) after going through RNA extraction with the Qiagen RNeasy mini prep kit and NEBNext Ultra II directional RNA library prep kit for Illumina + Poly(A) selection.
  
1. Trimming
   - Trimmed first 12 bp off each sequence due to quality and excluded any sequences with low quality using fastp for each sample. 
```
fastp -f 12 -i ./Forward/1_S36_L002_R1_001.fastq.gz -I ./Reverse/1_S36_L002_R2_001.fastq.gz -o 1F.fastq.gz -O 1R.fastq.gz -w 1
``` 
   - Result of trimming (1 sample):
```
Filtering result:
reads passed filter: 62920484
reads failed due to low quality: 106614
reads failed due to too many N: 754
reads failed due to too short: 13744
reads with adapter trimmed: 2891176
bases trimmed due to adapters: 31983123
``` 

2. Assembly with Trinity
   - Trinity run on all samples 
```
module load trinity/2.8.3
module load bowtie2/2.3.4.3

Trinity --seqType fq --SS_lib_type FR --left barbForwardfastp2.fastq.gz  --right barbReversefastp2.fastq.gz --max_memory 500G --CPU 24
``` 
   - Results:
```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  117271
Total trinity transcripts:      199546
Percent GC: 38.31

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 4574
        Contig N20: 3555
        Contig N30: 2932
        Contig N40: 2473
        Contig N50: 2084

        Median contig length: 571
        Average contig: 1113.20
        Total assembled bases: 222134494


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 4092
        Contig N20: 3075
        Contig N30: 2447
        Contig N40: 1940
        Contig N50: 1445

        Median contig length: 355
        Average contig: 736.71
        Total assembled bases: 86394222
```

3. BUSCO score
   - Code
```
module load busco/3.0.2

run_BUSCO.py -i Trinity.fasta -o barb_transcriptome_trinity_eudicot10 -l eudicotyledons_odb10 -m tran -c 5  
```
   - Results: C:96.8% [S:41.0%, D:55.8%], F:1.8%, M:1.4%, n:2121

### Genome annotation

  

