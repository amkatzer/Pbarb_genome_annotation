# Genome annotation for *Penstemon barbatus* 

## Contents
1. Introduction
2. Software links & citations
3. Annotation

## 1. Introduction
### Plant & Genome 
  *Penstemon barbatus* is a plant species found within Plantaginaceae. This species, like most within *Penstemon* is diploid with 8 pairs of chromosomes (Freeman, 1983). The expected length of the genome is 750 Mbps (Broderick *et al.* 2011).
  
### References
1. Freeman, C.C. Chromosome numbers in Great Plains species of Penstemon (Scrophulariaceae). *Brittonia* 35, 232–238 (1983). https://doi.org/10.2307/2806022
2. Broderick, Shaun R., et al. "A survey of Penstemon’s genome size." Genome 54.2 (2011): 160-173. https://doi.org/10.1139/G10-106 

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
### Helpful blogs & websites
1. https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2 **This is what a lot of my code follows when getting into running MAKER**
2. https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html#gsc.tab=0
3. MAKER tutorial wiki http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018
4. Botany 2020 workshop: https://github.com/bcbc-group/Botany2020NMGWorkshop

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
1. Run RepeatModeler to obtain repeat library
```
module load repeatmodeler/2.0.1

BuildDatabase -name pbarb_rep_mod -engine ncbi /panfs/pfs.local/scratch/hileman/a681k477/genome_barb/2020_barb_genome.fasta

RepeatModeler -engine ncbi -database pbarb_rep_mod -pa 20
```

3. MAKER
**RUN 1**
Initiate configuration files:
```
maker -CTL
```

Edit maker_opts.ctl: 
```
#-----Genome (these are always required)
genome=2020_barb_genome_A_ref_C_D_99_sealer_rnd1_scaffold.unmasked.fa #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=Trinity.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=uniprot_sprot_plants.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=consensi.fa.classified #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/panfs/pfs.local/software/7/install/maker/2.31.10/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

Run 1 MAKER: Use -fix_nucleotides flag if genome contains any non-ATCG nucledotide letters
```
module load maker/2.31.10

mpirun maker -base 20210621_maker6 -fix_nucleotides
```

Merge output files from all scaffolds:
```
module load maker/2.31.10

gff3_merge -s -d 20210621_maker6_master_datastore_index.log > rd_20210621_maker6.all.maker.gff
gff3_merge -d 20210621_maker6_master_datastore_index.log -g -n -o rd_20210621_maker6.genemodels.nofasta.noevidence$
gff3_merge -d 20210621_maker6_master_datastore_index.log -g -o rd_20210621_maker6.genemodels.fasta.noevidence.gff
gff3_merge -d 20210621_maker6_master_datastore_index.log -n -o rd_20210621_maker6.genemodels.nofasta.evidence.gff
gff3_merge -d 20210621_maker6_master_datastore_index.log -o rd_20210621_maker6.genemodels.fasta.evidence.gff
fasta_merge -d 20210621_maker6_master_datastore_index.log
```

Use output files from run 1 to train SNAP
   - filter to AED 0.25> and amino acid length 50+
   - Following this blog: https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
```
module load maker/2.31.10

maker2zff -x 0.25 -l 50 ../20210621_maker6.maker.output/20210621_maker6_master_datastore_index.log 

fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1

fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

forge export.ann export.dna > forge.log 2>&1

hmm-assembler.pl barb . > ./barb.hmm
```
Gene-stats output:
```
MODEL11014 skipped due to errors
MODEL16543 skipped due to errors
MODEL17223 skipped due to errors
MODEL33208 skipped due to errors
MODEL33207 skipped due to errors
MODEL33209 skipped due to errors
MODEL32959 skipped due to errors
MODEL32969 skipped due to errors
27 sequences
0.346517 avg GC fraction (min=0.323717 max=0.373614)
25372 genes (plus=12511 minus=12861)
1631 (0.064283) single-exon
23741 (0.935717) multi-exon
211.657654 mean exon (min=1 max=6556)
357.998138 mean intron (min=20 max=9991)
```

MAKER run 2:
Use output from run 1 in run 2. No longer need to incorporate transcriptome fasta or RepeatModeler library. 
*Ab initio* gene predictors trained: SNAP trained with run 1; Augustus trained with default tomato 

Generate configuration files
```
maker -CTL
```

Separate out est2genome, protein2genome, & repeats from aggregated gff file
```
awk '{ if ($2 == "est2genome") print $0 }' 20210621_maker6.all.maker.gff > 20210621_maker6.all.maker.est2genome.gff

awk '{ if ($2 == "protein2genome") print $0 }' 20210621_maker6.all.maker.gff > 20210621_maker6.all.maker.protein2genome.gff

awk '{ if ($2 ~ "repeat") print $0 }' 20210621_maker6.all.maker.gff > 20210621_maker6.all.maker.repeats.gff
```

Edit maker_opts.ctl
```
#-----Genome (these are always required)
genome=2020_barb_genome_A_ref_C_D_99_sealer_rnd1_scaffold.unmasked.fa #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=20210621_maker6.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=20210621_maker6.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=20210621_maker6.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=maker6_aed0.25_length50_genome.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=tomato #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=3 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

Run 2 MAKER
```
module load maker/2.31.10

mpirun maker -base 20210624_maker7 -fix_nucleotides
```

Merge output 
```
module load maker/2.31.10

gff3_merge -s -d 20210624_maker7_master_datastore_index.log > 20210624_maker7.all.maker.gff
gff3_merge -d 20210624_maker7_master_datastore_index.log -g -n -o 20210624_maker7.genemodels.nofasta.noevidence.gff
gff3_merge -d 20210624_maker7_master_datastore_index.log -g -o rd_20210624_maker7.genemodels.fasta.noevidence.gff
gff3_merge -d 20210624_maker7_master_datastore_index.log -n -o rd_20210624_maker7.genemodels.nofasta.evidence.gff
gff3_merge -d 20210624_maker7_master_datastore_index.log -o 20210624_maker7.genemodels.fasta.evidence.gff
fasta_merge -d 20210624_maker7_master_datastore_index.log
```
- note that there will be more files from this code than in the first MAKER run

Run SNAP on run 2 to look at gene stats
```
module load maker/2.31.10

maker2zff -x 0.25 -l 50 ../20210624_maker7.maker.output/20210624_maker7_master_datastore_index.log 

fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
```

SNAP gene-stats output:
```
MODEL11803 skipped due to errors
MODEL7466 skipped due to errors
27 sequences
0.346517 avg GC fraction (min=0.323717 max=0.373614)
13842 genes (plus=6881 minus=6961)
132 (0.009536) single-exon
13710 (0.990464) multi-exon
195.739777 mean exon (min=1 max=6538)
376.232086 mean intron (min=4 max=14776)
```

**To be continued!**

4. Weird MAKER troubleshooting solutions
   - Running MAKER with MPI will sometimes cause issues with the master_datastore_index.log file. To fix delete the log file and run:
   ```
   maker -dsindex -fix_nucleotides -base [same name as run output] 
   ```
   - If you are having issues with a scaffold failing due to this error "Can't kill a non-numeric process ID ... line 1050" this is a memory problem. You will need to re-run the scaffold with more memory. If you have run out of attempts, change the number of attempts in the maker_opts.ctl file and re-start the run. Only that scaffold will run. http://yandell-lab.org/pipermail/maker-devel_yandell-lab.org/2017-September/002944.html 
   - If at any point you hit a time limit but MAKER is not finished, you can stop and re-start. There are checkpoints that MAKER will use so it will not start from where it stopped last instead of the beginning.   
