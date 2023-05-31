# Contrasting new and available reference genomes to highlight uncertainties in assemblies and areas for future improvement: an example with Monodontid species

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        NA   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288

- [Abstract](#abstract)
- [Summary](#summary)
- [Overview of workflow](#overview-of-workflow)
- [Prepping files for input](#prepping-files-for-input)
- [Assembly of read datasets using Flye](#assembly-of-read-datasets-using-flye)
- [Polishing the assemblies using short reads](#polishing-the-assemblies-using-short-reads)
- [Scaffold using Juicer and 3D-DNA workflows](#scaffold-using-juicer-and-3d-dna-workflows)
- [Assemble RNA-seq data for EST predictions](#assemble-rna-seq-data-for-est-predictions)
- [Assess genome quality](#assess-genome-quality)
- [Predicting repeat elements](#predicting-repeat-elements)
- [Predicting gene models](#predicting-gene-models)
- [Citations](#citations)

# Abstract
Reference genomes provide a foundational framework for evolutionary investigations, ecological analysis, and conservation science, yet uncertainties in the assembly of reference genomes are difficult to assess, and by extension rarely quantified. We generated three beluga (Delphinapterus leucas) and one narwhal (Monodon monoceros) reference genomes and contrasted these with published chromosomal scale assemblies for each species to quantify discrepancies associated with genome assemblies. The new reference genomes achieved chromosomal scale assembly using a combination of PacBio long reads, Illumina short reads, and Hi-C scaffolding data. For beluga, we identified discrepancies in the order and orientation of contigs in 2.2-3.7% of the total genome depending on the pairwise comparison of references. In addition, unsupported higher order scaffolding was identified in published reference genomes. In contrast, we estimated 8.2% of the compared narwhal genomes featured discrepancies, with inversions being notably abundant (5.3%). Discrepancies were linked to repetitive elements in both species, indicating additional layers of data providing information on ultra-long genomic distances are needed to resolve persistent errors in reference genome construction. We present a conceptual summary for improving the accuracy of reference genomes with relevance to end-user needs and how they relate to levels of assembly quality and uncertainty.

# Summary 
Commands were carried out on Compute Canada's Cedar cluster using high performance computing job submissions/interactive sessions; command lines are provided, but note other SBATCH and module information was provided for each job (third party program version numbers provided here). Sample S_20_00708 (narwhal) is provided as an example here; the same set of commands were used for other samples (i.e. S_20_00693, S_20_00702, S_20_00703). Depending on the complexity of assembly, end-users will need to provide tailored resource requests (i.e. cores/threads+wall time). A major barrier to implimenting this workflow is the proper installation and compatability of dependencies. This was facilitated, in large part, through Compute Canada's Cedar cluster; users may experience a different level of frustration depending on what level of support they have, or how inherently stubborn they are. Note also, there is ample opportunity to pipe commands. As an initial run, my preference was to write intermediate files to disk in order to troubleshoot.

# Overview of workflow
![Mermaid_flow_assembly_diagram2_24v23](https://github.com/tbringloe/Monodontid_assemblies_2023/assets/79611349/59d47e14-fe95-4289-bd16-d63285bed722)

# Prepping files for input
### Download HI-C data. 
Files were downloaded July 4th, 2022.
```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/000/SRR8568870/SRR8568870_2.fastq.gz -o SRR8568870_Dovetail_Hi-C_of_the_narwhal_juvenile_male_liver_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR856/000/SRR8568870/SRR8568870_2.fastq.gz -o SRR8568870_Dovetail_Hi-C_of_the_narwhal_juvenile_male_liver_2.fastq.gz
```
### Download EST evidence (RNA-seq data). 
Files were downloaded July 5th, 2022.
```
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/003/SRR8578703/SRR8578703_1.fastq.gz -o SRR8578703_RNA-Seq_of_Monodon_monoceros_adult_male_spleen_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/003/SRR8578703/SRR8578703_2.fastq.gz -o SRR8578703_RNA-Seq_of_Monodon_monoceros_adult_male_spleen_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/002/SRR8578702/SRR8578702_1.fastq.gz -o SRR8578702_RNA-Seq_of_Monodon_monoceros_adult_male_lymph_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/002/SRR8578702/SRR8578702_2.fastq.gz -o SRR8578702_RNA-Seq_of_Monodon_monoceros_adult_male_lymph_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/005/SRR8578705/SRR8578705_1.fastq.gz -o SRR8578705_RNA-Seq_of_Monodon_monoceros_adult_male_liver_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/005/SRR8578705/SRR8578705_2.fastq.gz -o SRR8578705_RNA-Seq_of_Monodon_monoceros_adult_male_liver_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/004/SRR8578704/SRR8578704_1.fastq.gz -o SRR8578704_RNA-Seq_of_Monodon_monoceros_adult_male_kidney_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/004/SRR8578704/SRR8578704_2.fastq.gz -o SRR8578704_RNA-Seq_of_Monodon_monoceros_adult_male_kidney_2.fastq.gz
```
### Check quality of reads and trim prior to mapping/assembly
Here, a sample.list file with prefixes was created. This is fed into the xargs command, which allows for parallel processing of files depending on the number of threads available, specified with -P; the command to execute is listed in parentheses. For trimmomatic, illumina adapters were removed using specified ill_adap_uni.fa file with adapter sequences (provided here in the repository), trailing (3') sequences were trimmed once quality dipped below 25 in a sliding 10 bp window (SLIDINGWINDOW), the first ten basepairs were trimmed (HEADCROP), sequences with average quality less than 25 were discarded (AVGQUAL), and sequences less then 75 bp were discarded (MINLEN). The last 5 basepairs were also cropped through another round of trimmomatic using CROP (note, crop removes 3' end to achieve specified length, so value is desired read length not amount to crop!). These parameters were tailored somewhat depending on the read dataset (some had poorer quality reads, some had different lengths, ect). Results for fastqc were summarized using multiqc.
```
module load fastqc/0.11.9
module load trimmomatic/0.39
module load mugqic/MultiQC/1.12
cat sample.list | xargs -I {} -n 1 -P 4 sh -c "fastqc {}_1.fastq.gz {}_2.fastq.gz -o ."
multiqc .
##check results and adjust trimmomatic parameters as needed
cat sample.list | xargs -I {} -n 1 -P 4 sh -c "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {}_1.fastq {}_2.fastq trimmed/{}_1_trimmed.fq.gz trimmed/{}_1_tup.fq.gz trimmed/{}_2_trimmed.fq.gz trimmed/{}_2_tup.fq.gz ILLUMINACLIP:ill_adap_uni.fa:2:30:10 SLIDINGWINDOW:10:25 HEADCROP:10 AVGQUAL:25 MINLEN:75"
cat sample.list | xargs -I {} -n 1 -P 4 sh -c "fastqc trimmed/{}_trimmed_1.fq.gz trimmed/{}_trimmed_2.fq.gz -o trimmed"
multiqc trimmed
```

### Convert raw long read BAM files to fastq
BAM files provided by genome quebec were initially sorted then converted to fastq for analysis using samtools; -@ specifies number of threads to use, -m specifies memory allocation for threads, -n specifies to sort by name, -o specifies output file. In some cases, files were moved to cedar cluster from DFO servers using scp (secure copy protocol)

```
module load samtools/1.13
samtools sort -@ 32 -m 4G -n -o 0422_sorted.bam *422.subreads.bam 
samtools bam2fq -@ 32 0422_sorted.bam > 0422.fastq
##If moving files, compress then use scp and enter password promt
gzip 0422.fastq
scp 0422.fastq.gz tbringlo@cedar.computecanada.ca:~/MOBELS-WGS.2
```

### Evaluate read length distributions and subsample
Read lengths were subsampled for the longest reads totalling 200 billion base pairs, representing ~50x coverage, ranging from 9560 kbp and 51% of the dataset, to 11210 kbp and 80% of the dataset (i.e. samples S_20_00708, S_20_00693, S_20_00702, S_20_00703). For seqkit, seq function specifies to transform sequences, -m specifies minimum read length to filter and output

```
module load bbmap/38.86
module load seqkit/0.15.0
readlength.sh in=0422.fastq.gz out=0422_lengths # use bbmap to generate read length distribution information, output is evaluated and informs -m argument for seqkit
seqkit seq -m 10960 0422.fastq.gz > 0422_10k.fastq # removes reads less than specified read length
readlength.sh in=0422_10k.fastq out=0422_10k_lengths # use bbmap to verify read file was properly subsampled and is ready for assembly; alternatively check flye log file once assembly initiates
```
# Assembly of read datasets using Flye
Flye is a long read assembler that creates an initial repeat graph with inexact pathways through contigs (called disjointigs), which is corrects through read mapping to disentangle pathways through the assembly. In order to construct the graph, contigs are extended by overlapping reads; because the minimum read length used for assembly was ~10kbp, this is the minimum overlap length set by Flye for extending contig sequences. The result is (ideally) a highly contiguous assembly at the contig level. The assembler is appropriate for small to mammalian sized genomes, and takes raw sequence data as input. See the publication in nature for more details: https://www.nature.com/articles/s41587-019-0072-8; and manual: https://github.com/fenderglass/Flye. --pacbio-raw specifies raw read dataset, -t specifies threads, -i specifies number of polishing iterations to perform (none here since we will use polish downstream), --scaffold specifies to allow scaffolding (very few scaffolds produced in output and could be removed from command line), --out-dir specifies output directory. Note, some assemblies were restarted due to wall timeout errors; these are reinitiated at last checkpoint using --resume

```
#Using local install of Flye v.2.9
Flye/bin/flye --pacbio-raw 0422_10k.fastq -t 32 -i 0 --scaffold --out-dir Flye_output_v1_0422_raw
```

# Polishing the assemblies using short reads
The assembly of raw pacbio reads still contains errors, mostly indels and incorrect substitutions (see: https://www.nature.com/articles/s41587-018-0004-z). These can be corrected with accurate short reads. The program pilon considers pileups of short reads at all genome positions to correct the sequence data. Only reads with a single mapping position are considered by the pileup process, and changes are only made when the evidence overwhelmingly supports it. Because corrected errors can change ambiguously mapped reads to a single location, multiple iterations of polishing are recommended to "converge" on a corrected assembly (however, the return on efforts diminishes rapidly after ~2 rounds). See pilon publication (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963) and manual (https://github.com/broadinstitute/pilon/wiki) for more details. Files were first prepped for Pilon by mapping short reads to the Flye assembly. Following this, files were converted to bam format, sorted, then indexed. The same set of commands were applied to the polished genome (i.e. 2 rounds of polishing performed). For bowtie2 commands, --score-min sets parameters for mapping thresholds (see online manual), -p specifies number of threads, --no-unal specifies to not output reads that did not map, -x specifies the index to use (i.e. reference genome), -1 and -2 specifies forward and reverse reads, respectively, -s specifies output; --For Pilon commands, --diploid specifies reference is from a diploid individual, --iupac specifies to use ambiguity codes (not used here), --output specifies output header, --outdir specifies output directory. Note, the introduction of ambiguity codes through the --iupac flag does introduce complications downstream, as many programs will switch non ACTG characters to N; consider removing this flag if heterogeneity is not of interest.
```
module load bowtie2/2.4.4
module load samtools/1.13
#Build bowtie2 index
bowtie2-build assemblies/Flye_output_v2_0422_raw/0422_Flye_V2_assembly.fasta short_Reads/IDX/Narwhal
#map reads to reference genome using standard thresholds (~10% divergence in high quality reads)
bowtie2 --score-min L,0,-0.6 -p 8 --no-unal -x short_Reads/IDX/Narwhal -1 short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035_R1.fastq.gz -2 short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035_R2.fastq.gz -S short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035.sam

#prep bowtie2 output for pilon using samtools
samtools view -S -@ 8 -b short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035.sam > short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035.bam
samtools sort -@ 8 -o short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035_s.bam short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035.bam
samtools index -@ 8 short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035_s.bam

#Pilon commands
module load pilon/1.24
java -Xmx150g -jar $EBROOTPILON/pilon.jar --genome assemblies/Flye_output_v2_0422_raw/0422_Flye_V2_assembly.fasta --frags short_Reads/NS.1801.IDT_i7_162---IDT_i5_162.ADN_20_01035_s.bam --diploid --output Flye_output_v2_0422_raw_pilon_polished --outdir assemblies/0422_Flye_V2_assembly/polished_pilon
```

# Scaffold using Juicer and 3D-DNA workflows
Scaffolding is arguably the most important but least straightforward step towards building a reference genome. I tested multiple scaffolding strategies, including reference based scaffolding using RagTag, Hi-C contacts using pin_hic, and k-mer based approaches using RAILS. In the end, there was no comparison to the rigor and transparent documentation of the Juicer+3D-dna workflows. At its core, juicer+3d-DNA are simply layers of bash scripts corresponding to several steps. Those steps take hi-c reads, align them to a reference genome, remove duplicate and chimeric reads, perform iterative scaffolding based on contract densities, detect and correct misjoin errors (regardless of where they were introduced in the workflow), polish scaffolding based on known errors (such as spurious high densities at telomeres), collaps heterozygous regions, and produce a final set of chromosomes; any non-scaffolded sequences are put aside and appended to the fasta file as debris. The user can view .hic files at each step by loading into the java program juicebox. Users can optionally manually edit their assembly, which was done here to correct obvious misjoins and non-concatenated scaffolds missed by the workflow. The workflow was developed under the mandate to produce 1000$ genomes, that is, using cheap short reads to produce draft assemblies and scaffold with the juicer+3D-dna workflows; an "assembly" cookbook has even been offered for prospective users (https://aidenlab.org/assembly/manual_180322.pdf). See also, the preprint here: https://www.biorxiv.org/content/10.1101/254797v1; and the official paper with supplemental material here: https://www.science.org/doi/10.1126/science.aal3327, (Dudchenko et al. 2017; 2018). Github pages here: https://github.com/aidenlab/juicer, https://github.com/aidenlab/3d-dna, and https://github.com/aidenlab/Juicebox.

```
#Users will need samtools, java, and bwa in their filepath; if leveraging gpu capabilities, other dependencies may also be required.
module load StdEnv/2020
module load gcc/9.3.0
module load bwa/0.7.17
module load java/17.0.2
module load nvhpc/22.7
module load cuda/11.7
module load samtools/1.13
module load python/3.10.2
#Define this variable to avoid changing throughout, use prefix before .fasta of assembly
asm=S_20_00708_PilonPolished2_contigs
#The workflow leverages cleave site information based on the restriction enzyme used, so create a sites file, along with a file specifying lengths of contigs/scaffolds/chromosomes being fed into the workflow; generate_site_positions.py <restriction enzyme used> <prefix for new files> <draft reference genome>
python generate_site_positions.py MboI "$asm" references/"$asm".fasta
#generate chrom.sizes file
awk 'BEGIN{OFS="\t"}{print $1, $NF}' "$asm"_MboI.txt > "$asm".chrom.sizes
#index the reference using bwa
bwa index reference/"$asm".fasta
#Run Juicer workflow to generate input files for 3D-dna; see github page for info on how to setup the directory and filepaths; once setup the sed command can replace a modified ASM_DIR in the specified file path, which again helps swapping out all relevant code for specified sample.
sed -i "s/ASM_DIR/"$asm"/g" juicer_CPU_modified.sh
#I used the CPU script, this threw some errors downstream, but safely made its way through the mapping, chimeric flagging and deduplication steps. The output should be a merged_nodups.txt file, which is fed into 3d-DNA along with the reference
#-g specifies genome prefix for files, -s is the restriction enzyme used to generate hi-c library, -z is path to reference+index, -t specifies threads, -p specifies path to size file, -y specifies path to the restriction sites file, generated using python, --assembly should tell the workflow to produce the merge_dedups.txt file and exit early; -S flag specifies the stage to start/end at, probably redundant with previous flag.
bash juicer_CPU_modified.sh -g "$asm" -s MboI -z references/"$asm".fasta -t 32 -p references/"$asm".chrom.sizes -y restriction_sites/"$asm"_MboI.txt --assembly -S early

#Next run 3D-dna workflow using juicer output
bash 3d-dna/run-asm-pipeline.sh "$asm".fasta "$asm"_merged_nodups.txt
#Since the above juicer workflow will exit early, if needed one can use the visualize script with 3D-dna to produce .hic and .assembly files to load into juicebox.
#First generate the .assembly file
awk -f 3d-dna/utils/generate-assembly-file-from-fasta.awk "$asm".fasta > "$asm"_draft.assembly
#Generate .hic file; the specifies commands are defaults, see run-asm-pipeline.sh file for descriptors
bash 3d-dna/visualize/run-assembly-visualizer.sh -p true -q 1 -i -c $asm".assembly "S_20_00703_PilonPolished2_contigs"$asm"_merged_nodups.txt
##Load raw chromosome assembly into juicebox and manually edit obvious scaffolding errors, save as a review.assembly file, then use post-review script to produce final assembly, and use -g flag to insert 100bp Ns in gaps
bash 3d-dna/run-asm-pipeline-post-review.sh -g 100 --sort-output -r "$asm".rawchrom.review.assembly "$asm".fasta "$asm"_merged_nodups.txt

#generate .agp file for ncbi submission using juicebox_assembly_converter.py script provided by https://github.com/phasegenomics/juicebox_scripts
```
*It is highly advisable at this point that you name your sequences according to NCBI standards. Namely, any relevant information should be incorporated into the sequence header (e.g. HiC_scaffold_19 [organism=Delphinapterus leucas] [isolate=S_20_00703] [chromosome=9]). Failure to adhere to NCBI standards will create errors in your submission that require correction before approval.*

# Assemble RNA-seq data for EST predictions
### Assembly transcripts using RNA-spades
Publically available RNAseq data for beluga and narwhal were assembled de novo using SPAdes-rna. --rna specifies to call on rna assembler, -t specifies number of threads, -o specifies output directory (see rna-SPAdes manual here: https://cab.spbu.ru/files/release3.11.1/rnaspades_manual.html). Trinity is another industry standard and likely does a better job at resolving isoforms; because the RNAseq data is simply being used as EST evidence, isoforms were not of specific interest. The transcripts.fa output file was used for further analysis. 
```
#Run SPAdes-rna, looping over sample IDs stored in text file "list"
module load StdEnv/2020
module load spades/3.15.4
cat list | while read line
do
mkdir rna_spades_"$line"
spades.py --rna -s "$line".fq.gz -t 8 -o rna_spades_"$line"
done
```
### Clean transcripts for contaminant vectors
Sequence vectors were also cleaned from the assembled transcripts using seqclean and the UINVEC database. 
```
##Download UniVec database and screen transcripts to remove vector contamination (e.g. viral inserts, but UniVec also contains adapters, linkers and primers)
##Seqclean is distributed through the pasa pipeline, but here I simply downloaded the specific script here:https://sourceforge.net/projects/seqclean/files/
module load blast+/2.12.0
mkdir UNIVEC/
cd UNIVEC/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv

makeblastdb  -in UniVec -dbtype nucl
makeblastdb  -in UniVec_Core -dbtype nucl

##Because I had multiple assembled datasets I looped the command over a list containing all fasta file names
cat narwhal.list | while read line
do
seqclean-x86_64/seqclean $line -v UNIVEC/UniVec_Core -c 4
done
```
### Deduplicate the transcripts using mmseqs2
In order to reduce computation downstream when aligning transcripts to the reference genome for predicting genes, transcripts from all the assemblies were pooled and clustered using mmseqs2 (https://mmseqs.com/latest/userguide.pdf), such that a single representative sequence was derived for annotations. Clustering mode will aggregate sequences that overlap with query by 80% (-c flag) and have 99% similarity (--min-seq-id flag). --cov-mode 1 sets coverage mode to specifically apply to the target sequence (must overlap with query be specified amount, reciprocal does not apply), --cluster-mode prioritizes larger sequences as query (i.e. hopefully full length transcripts in this case), --threads specifies number of threads to use. If computational resources are limiting, a similar set of commands can be run by swapping out mmseqs cluster with mmseqs linclust, which is a less sensitive clustering algorithm, but orders of magnitude quicker. In my experience, both commands produced similar sized datasets, linclust being a bit larger.
```
##Use mmseqs2 to cluster and derive representative transcripts from the pooled assemblies
module load mmseqs2/13-45111
#concatenate all the transcript output from rna_SPAdes runs
cat Narwhal&ast;_transcripts.fasta > Narwhal_combined_transcripts_13ix22.fasta
mmseqs createdb Narwhal_combined_transcripts_13ix22.fa Narwhal_transcripts_db
#Store temporary files in temporary directory
mkdir tmp_dir
mmseqs cluster Narwhal_transcripts_db Narwhal_transcripts_db_clu tmp_dir --threads 4 --cov-mode 1 -c 0.80 --cluster-mode 2 --min-seq-id 0.99
mmseqs result2repseq Narwhal_transcripts_db Narwhal_transcripts_db_clu Narwhal_transcripts_db_clu_rep
mmseqs result2flat Narwhal_transcripts_db Narwhal_transcripts_db Narwhal_transcripts_db_clu_rep Narwhal_transcripts_rep.fasta --use-fasta-header
```

# Assess genome quality
### Run BUSCO to assess geneome completeness and for downstream ab initio gene predictions
BUSCO screens assembled genomes for single copy genes from specified databases, and uses these hits as a proxy for the "completeness" of the genome. Several relevant databases are available for Narhwal and Beluga. Cetartiodactyla includes even-toed ungulates and cetaceans, and contains ~13,000 BUSCOs at its latest version (dated 2019-11-20); a mammalian database also contains ~10k BUSCOs. I found however, busco scores were low (~50%) when screening against these databases, regardless of whether I screened newly assembled or already published genomes for narwhal and beluga. So I specified the vertebrata database (~3500 BUSCOS) and scores were 90+%, suggesting the current mammalian and cetartiodactyla databases are not representative of cetaceans. For busco commands, -c specifies number of threads, -f force overwrites a previous run at output path, -r restarts a partially completed run, -i specifies input genome, -o specifies output prefix/folder, -l specifies lineage to use, --out_path specifies output path
```
#seach for BUSCOs using the cetartiodactyla database
module load busco/5.2.2
busco -c 8 -m genome -f -r -i 0422_Flye_V2_assembly.fasta -o 0422_flye_v2_raw_busco -l vertebrata_odb10 --out_path ~/MOBELS-WGS.2/BUSCO/0422_flye_v2_raw_busco
```

### Compare new assemblies to current published ones
Part of the motivation to assemble new reference genomes was to overcome limitations in the currently published genomes. Beluga was assembled from short read data and iteratively scaffolded in the absence of information to bridge repeat regions, while the narwhal genome assembly was guided by the (likely erroneous) beluga reference genome. Quast is a program that allows user to compare an assembly to reference data, and determine the number of inconsistencies between the two assemblies. Misassemblies are categorized according to inversions (segments of inverted DNA), relocations (segments of DNA broken and spaced apart on the same reference scaffold), and transpositions (segments of DNA that appear on different reference scaffolds). The total length of misassembled genome is also reported, but it simply adds the length of all contigs/scaffolds with misassemblies identified, the true amount of misassembled genome would be a lot less. See Quast manual here: http://quast.sourceforge.net/docs/manual.html. For Quast flags: -t specifies number of threads, -r specifies the reference genome to compare to (here, the newly assembled genome was more contiguous, so I specified as the new narwhal genome), --large specifies the genomes are large (i.e. mammalian) genomes, --fragmented specifies the reference genome is fragmented (quast will account for this when identifying misassemblies), -o specifies prefix for output folder, the last bit points to the genome to compare to the reference. 
```{bash, eval=FALSE}
module load quast/
#compare new assembly to reference genome of Westbury et al. 2019
quast -t 8 -r S_20_00703_MOBELS-0202_Beluga_final_2ix22.fasta --large --fragmented -o quast_Beluga_0202_final_v_Jones GCA_002288925.3_ASM228892v3_genomic.fna
```

# Predicting repeat elements
### Ab initio predict repeat elements
Before predicting genes, it is imperative to predict and mask repeat areas (otherwise this will lead to spurious alignments and gene predictions). Repeat modeler is a workflow that leverages several repeat prediction algorithms, including repeat scout, RECON, tandem repeat finder, and optional long tandem repeat finder LTR retriever. Several other dependencies are also required (and a slight pain, but not impossible, to setup); see online information here for Repeat Modeler:https://www.repeatmasker.org/RepeatModeler/. At its core, repeat modeler is an ab initio predictor. It scans subsets of the genome in multiple rounds to predict repeat elements/families, under the assumption that repeat elements are randomly dispersed throughout the genome; for more complex and larger genomes, tuning parameters to subset a larger portion of the assembly should, in theory, improve predictions. Note, however, that because a given round informs the subsequent round, the algorithm bloats quickly with temporary files and used disk space, and round 6 can take days to complete.A concensus file of the predicted repeat elements is output at each round, so if the algorithm is not allowed to complete, you can use what was predicted up until that point. The LTR predictions are tacked on as an additional analysis at the end of the workflow, however. -pa specifies number of threads to use, -name specifies database prefix, -database points to the created database, -LTRStruct specifies to enable LTR predictions at end of round 6.

```
module load perl/5.30.2
module load bedtools/2.30.0

#Define some variables for flags below
export GENOME_NAME=S_20_00708_MOBELS-0422_Narwhal_final_2ix22.fasta
export GENOME_SIZE=2336652991
export SLURM_NTASKS=4 # repeat modeler does not use multiple threads efficiently due to some limitations with dependencies used, 4 threads will speed up computation without wasting a lot of resources. Make sure to set a large time wall limit; for large mammalian genomes repeat modeler run time could easily exceed 5 days.

#Running Repeat Modeler, first building the database, then running for ab initio predictions
RepeatModeler-2.0.3/BuildDatabase -name ${GENOME_NAME}_db ${GENOME_NAME} > BuildDatabase_narwhal.out
RepeatModeler-2.0.3/RepeatModeler -pa ${SLURM_NTASKS} -database ${GENOME_NAME}_db -LTRStruct 1> RepeatModeler.sdtout 2> RepeatModeler.sdterr
```

### Run Repeat Masker
Repeat masker comes with the Dfam database out of the box, but more comprehensive libraries, including those from RepBase (if you have an arm and a leg to pay for access) can be specified during configuration. Here, repeats are identified from known sequences (dfam in this case), including newly predicted ones specific to the genome/species of interest (i.e. the predictions from repeat modeler, i.e. the consensi.fa output file). Repeat masker will output a gff file which can be used for softmasking with bedtool functions, as in, repeat regions will be switched to lower case in the genome fasta file, which some programs will recognise to not consider (alignments may still be allowed to flow into or overlap repeat regions). Perl scripts provided by Repeat Masker will also convert the repeat annotations to gff3 and create models of the repeat landscape the user can peruse if needed. -lib specifies the repeat library to consider (custom one in this case, as in library built at configuration+repeat modeler results), -e ncbi specifies seach engine, -gff specifies to output gff, -no_is skips bacterial insertion element check, -e specifies search engine (this command worked, but online the flag is -engine), -x returns repetitive regions masked with Xs rather than Ns (but remember, here, masking is applied with bedtools and MAKER downstream), -a shows alignments in a .align output file, -pa specifies number of threads to use. See repeat masker documentation for more details (https://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html). For Bedtools, -fi is the input genome, -soft initiates soft masking, -fo is the output genome, -bed specifies the bed file for masking, generated from the repeat masker analysis (https://bioweb.pasteur.fr/docs/modules/bedtools/2.17.0/content/maskfastafromBed.html).
```
export GENOME_NAME=S_20_00708_MOBELS-0422_Narwhal_final_2ix22.fasta
export GENOME_SIZE=2336652991
export SLURM_NTASKS=4 # repeat modeler does not use multiple threads efficiently due to some limitations with dependencies used, 4 threads will speed up computation without wasting a 
# Combine repeat modeler results (final output and LTR module). LTR results will be embedded in RM output folder. The user could consider combining with the built RepeatMask database (default dfam db), but I only used species specific results from RepeatModeler to reduce runtimes. As a consequence, there is likely a few percentages of repeats not detected in this analysis (confirmed checking previous narwhal assembly repeat results).
cat ${GENOME_NAME}_db-families.fa LtrRetriever-redundant-results.fa > ${GENOME_NAME}_CombinedRepeatLib.lib

# Run repeat masker
RepeatMasker/RepeatMasker -lib RepeatMasker/Libraries/RepeatMasker_custom.lib -e ncbi -gff -x -no_is -a -pa ${SLURM_NTASKS} ${GENOME_NAME} > RepeatMasker.sdtout 2> RepeatMasker.sdterr

# SoftMask Genome # This is probably the genome to consider publishing, with all repeat regions soft masked. The bed file is also useful for other analyses where users may wish to exclude repetive regions from their analysis ##Note, output does not appear to be softmasked, revisit this.
maskFastaFromBed -soft -fi ${GENOME_NAME} -fo ${GENOME_NAME}.softmasked -bed ${GENOME_NAME}.out.gff

# Get repeat features in gff3 format (used for EVM later on)
perl RepeatMasker/util/rmOutToGFF3.pl ${GENOME_NAME}.out > ${GENOME_NAME}.out.gff3

# Model repeat landscape (not needed for annotation)
perl RepeatMasker/util/calcDivergenceFromAlign.pl -s ${GENOME_NAME}.divsum ${GENOME_NAME}.align
perl RepeatMasker/util/createRepeatLandscape.pl -div ${GENOME_NAME}.divsum -g ${GENOME_SIZE} > ${GENOME_NAME}.html
```
### Prep repeat element files for MAKER2 predictions
Here, I follow the workflow of darencard (https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2). MAKER2, which will be used for gene annotations, will hard mask complex repeats; this can be detected internally, or it will automatically hard mask any repeat regions specified with a gff3. Because we already ran RepeatMasker, it makes sense not to double on computational efforts. So the below commands isolate complex repeats in gff3, and simple repeats in fasta format, which will be soft masked my MAKER2. Note, as observed here (https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html#gsc.tab=0), MAKER2 won't take RepeatMasker generated gff3 without some reformating. Below code is courtesy of darencard and his clever workarounds.
```
# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" ${GENOME_NAME}.out.gff3 > ${GENOME_NAME}.out.complex.gff3
# reformat to work with MAKER
module load perl/5.30.2
cat ${GENOME_NAME}.out.complex.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > ${GENOME_NAME}.out.complex.reformat.gff3

##isolate simple repeats in fasta format. In order for this to work the fasta file must have sequences on a single line. I imported into geneious and exported to correct this.
grep -A 1 -F -e "DNA" -e "Unknown" -e "rRNA" -e "tRNA" -e "Satellite" -e "Simple" ${GENOME_NAME}_db-families_single_line.fasta > ${GENOME_NAME}_db-families_tmp.fasta
#removes -- lines created using the -A grep flag
grep -v -e "--" ${GENOME_NAME}_db-families_tmp.fasta > ${GENOME_NAME}_db-families_simple.fasta
```

# Predicting gene models
### MAKER round 1: Generate gene predictions using MAKER2 workflow
Finally, we are ready to predict genes. The workflow detailed here is a mix of insights from these helpful tutorials :https://reslp.github.io/blog/My-MAKER-Pipeline/; https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html#gsc.tab=0; https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2. Similar to repeat modeler, the MAKER workflow integrates evidence across numerous datasets and third party programs, and makes predictions weighting all the evidence provided. The program also integrates repeat masker; the strategy used here took some detective work. Previously, I have annotated repeat regions outside the MAKER workflow (using repeat modeler+repeat masker). MAKER will classify repeats as simple or complex, and soft or hard mask these regions respectively. This is important because complex repeats (e.g. retrotransposons) can wreak havoc on gene predictions, so they should be removed from alignment altogther. On the other hand, simple repeats are in genes of interest, and should be allowed to align when making gene predictions. I cannot supply these regions to MAKER in gff, however, because all the repeats regions will be hard masked (confirmed here:https://groups.google.com/g/maker-devel/c/7UbOIvwaaRM). So ideally I would supply complex repeats in gff, and allow MAKER to identify simple repeats internally. Following this tutorial (https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2), I used repeat masker output from above to identify complex repeats and convert to gff with formatting appropriate for MAKER, and supplied the custom repeat library for Repeat Masker to consult within the MAKER workflow (thus allowing MAKER to classify the simple repeats). Another option is to specify simple as a RepBase library to consult (but I did not have a subsription). Or just specify the custom library and allow MAKER to classify the complex repeats. Looking towards our true goal of annotating genes, the first step is to get an initial set of gene predictions we can use to make ab initio predictions using several third part programs. Here, I am using the single copy BUSCOs identified previously as protein evidence, and a concatenated file of transcripts previously assembled using SPAdes-rna, and purged of vectors using seqclean. A first round of the MAKER workflow is initiated, and the output is used for ab initio predictions in other third party programs. MAKER documentation can be found here:http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Repeat_Masking
```
#module load maker/3.01.03
#First step is to modify control files used to specify inputs for the MAKER run. Generate a new set of control files with the following command
maker -CTL

#Next, change settings in the files generated. In particular, several settings should be set in the maker_opts.ctl file; the settings I specified are shown below.

#-----Genome (these are always required)
#genome=S_20_00703_MOBELS-0202_Beluga_final_2ix22.fasta #genome sequence (fasta file or fasta embeded in GFF3 file) # note, not the softmasked version
#organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Repeat Masking (leave values blank to skip repeat masking)
#rmlib=RepeatMasker/Libraries/RepeatMasker_custom.lib #provide an organism specific repeat library in fasta format for RepeatMasker
#rm_gff=mask.out.complex.reformat.gff3 #pre-identified repeat elements from an external GFF3 file
#softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----EST Evidence (for best results provide a file for at least one)
#est=Beluga_combined_transcripts_13ix22.fa #set of ESTs or assembled mRNA-seq in fasta format

#-----Protein Homology Evidence (for best results provide a file for at least one)
#protein=Beluga_final_2ix22_SCBUSCO.faa #protein sequence file in fasta format (i.e. from mutiple organisms)

#-----Gene Prediction
#est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
#protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no

#-----MAKER Behavior Options
#alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no

#Now run MAKER2 # Fix nucleotides switches any ambiguity codes to Ns, since MAKER won't interpret these. We shouldn't have any...The first round will align all protein and EST evidence using blast searches and tidy up alignments and gene predictions using exonerate. Repeat masker will soft mask simple repeats in these searches, while complex repeats are removed altogether.
maker -fix_nucleotides

#Convert output files into a format accepted by SNAP
gff3_merge -d my_assembly_datastore_index.log
```

### Generate ab-initio gene predictions by training SNAP
Using the output from an initial round of MAKER, we can provide hints to make ab initio gene models, starting with SNAP. SNAP documentation is thin, but it has been benchmarked among several other ab initio predictors, and it appears to perform well; see github information here:https://github.com/KorfLab/SNAP. 
```{bash, eval=FALSE}
#module to load
module load StdEnv/2020
module load gcc/9.3.0
module load snap/2017-05-17
module load perl/5.30.2
module load maker/3.01.03

#job commands
#create files for snap, added filters for gene models with length at least 50 amino acids and AED scores of at least 0.25 or lower (0 means perfect concordance with available evidence, 1 indicates no concordance); these filters toss out poor models, which should improve accuracy of SNAP predictions
maker2zff -x 0.25 -l 50 -d 00703_round1_master_datastore_index.log
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/fathom genome.ann genome.dna -gene-stats > snap_gene_stats_output.txt
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/fathom genome.ann genome.dna -validate > snap_validate_output.txt
cat snap_validate_output.txt | grep "error" > snap_validate_errors.txt
##create some line to grep and remove errors based on output
#grep -vwE "some text from error output" genome.ann > genome.ann2
grep -vwE "MODEL26093|MODEL31041|MODEL32008|MODEL16958|MODEL44350|MODEL37154|MODEL27287|MODEL27288" genome.ann > genome.ann2
#Rerun fathom and confirm there are no errors
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/fathom genome.ann2 genome.dna -validate > snap_validate2_output.txt

#Now create input files for training snap. These commands will take gene identified by SNAP as unique and complete, that is, no isoforms and the sequence (nucleotide/protein) contains start and stop codons. The commands also export these gene regions and the bordering 1000 bp on either side, which will be used downstream for training augustus.
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/fathom genome.ann2 genome.dna -categorize 1000 > categorize.log 2>&1
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/snap/2017-05-17/bin/forge export.ann export.dna > forge.log 2>&1

#Now run snap on input files
hmm-assembler.pl 00703_snapr1 . > 00703_snap1.hmm
```

### Generate ab-initio gene predictions by training Augustus
We will also use Augustus to make gene predictions. Augustus is a generalized hidden Markov model that provides probability distributions for sequences, including introns, exons and intergenic regions. It is built into several third party programs, such as BUSCO, and is a leading standard for ab initio gene prediction models. The model output from Augustus can be fed into the MAKER workflow. First we need to generate the model. Here is a tutorial for training augustus: https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html
```
module load augustus/3.4.0
module load bioperl/1.7.8
#First generate a genebank format file to input into Augustus, running a handy perl script in the directory with the export.ann and export.dna files produced by snap fathon; script here:https://github.com/hyphaltip/genome-scripts/blob/master/gene_prediction/zff2augustus_gbk.pl
perl zff2augustus_gbk.pl > augustus.gbk

#split the august.gbk file into training and test datasets. This script (and others specified below) comes with Augustus and is in the scripts folder
perl 3.4.0/scripts/randomSplit.pl augustus.gbk 100

#create a new species for training Augustus
perl 3.4.0/scripts/new_species.pl --species=narwhal_00708_r2 --AUGUSTUS_CONFIG_PATH=3.4.0/config/

#Train Augustus with training set file, then evaluate results
etraining --species=narwhal_00708_r2 augustus.gbk.train --AUGUSTUS_CONFIG_PATH=3.4.0/config/
augustus --species=narwhal_00708_r2 augustus.gbk.test --AUGUSTUS_CONFIG_PATH=3.4.0/config/ | tee first_training.out

#Optimize prediction parameters
#Default kfold value is 8, setting CPU limit to 8, so if specifying more threads be sure to increase number of kfolds (--kfold <value>); optimization takes about 1-2 days
perl 3.4.0/scripts/optimize_augustus.pl --cpus=8 --species=beluga_00703 augustus.gbk.train --AUGUSTUS_CONFIG_PATH=3.4.0/config/

#Retrain Augustus and compare results to the first round. This can take a while; on my second round of retraining Augustus this ran for over a week on 8 cpus...
etraining --species=narwhal_00708_r2 augustus.gbk.train --AUGUSTUS_CONFIG_PATH=3.4.0/config/
augustus --species=narwhal_00708_r2 augustus.gbk.test --AUGUSTUS_CONFIG_PATH=3.4.0/config/ | tee second_training.out

#If confident with model, move to next round of MAKER predictions. Augustus webpage suggests models with gene prediction sensitivity <20% are probably not good enough to use, though there is plenty of debate to be had about that, especially within the context of the maker workflow (which provides hints to Augustus about gene location based on other data, and weighs different annotation types to make final predictions, meaning Augustus is one source of information among many). 
```

### MAKER round 2: generate a second round of gene predictions
Now we can run MAKER again, this time specifying the output from round 1 and the trained models for SNAP and Augustus. Because we already used the est and protein evidence in round 1, those results can be specified in the round1 gff, rather than the fasta files (thus avoiding costly realignments), so swap out this evidence.
```
#get annotations into a single file to input into second round

#Change the following lines in the maker_opts.ctl file
#-----Re-annotation Using MAKER Derived GFF3
maker_gff=maker.gff.round1 #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #removing previous input here
#-----Protein Homology Evidence (for best results provide a file for at least one)
protein= #also removing evidence from round 1
#-----Gene Prediction
snaphmm=~/MOBELS-WGS.2/phase4/Final_assemblies/maker_workflow/beluga/00703_round2/00703_snap_r1_24i23.hmm
augustus_species=beluga_00703 #Augustus gene prediction species model
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no

#Once you setup the new environment with new evidence and fresh control files, run MAKER2 for a second round
maker -fix_nucleotides
```

## Make a second round of models with SNAP and AUGUSTUS
MAKER can be run iteratively, in theory improving gene predictions with each round of model training and evidence weighting by MAKER. Different tutorials suggest a different number of iterations (2-4 rounds of MAKER). Here, I used iterative MAKER results to train two rounds of SNAP and Augustus models, making for three rounds of MAKER predictions. Simply repeat the above steps to get a second round of gene models and a final round of MAKER predictions. Returns on gene prediction accuracy should diminish after 3 rounds of MAKER, and there is the danger of overtraining SNAP and Augustus.

### Evaluate gene model results
We can determine if the predictions have converged using three strategies: 1) by looking at gene numbers and lengths after each round of MAKER; 2) look at AED distributions, where 0 to 1 quatifies the confidence in a gene model, and 3) visualize the gene models in a genome browser.
```
#First we must extract annotations and transcripts for downstream functional annotations, launched within output directory
#Fasta_merge will create protein and nucleotide sequences of gene models
gff3_merge -d <sample>_datastore_index.log
fasta_merge -d <sample>_datastore_index.log
```
I take some leads from this [tutorial](https://github.com/Joseph7e/MAKER-genome-annotations-tutorial). Next we can count the number of gene models and lengths after each round of maker. These should converge somewhat after 3 rounds. In our case there was a difference of 4,497 predicted genes between rounds 1 and 2, and 407 between rounds 2 and 3.
```
cat <roundN.full.gff> | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
```

Next the distribution of AED scores are visualised with a with perl script. Lower scores indicate more empirical evidence is consistent with a predicted gene. At the expense of repeating what I have seen online, 95% of scores should be below 0.50 for a well annotated genome. For our work here, 93.5% of gene models had an AED score of 0.5 or less. The script to calculate the accumulation values can be found [here](https://github.com/mscampbell/Genome_annotation/blob/master/AED_cdf_generator.pl).
```
perl AED_cdf_generator.pl -b 0.025 <roundN.full.gff>
```
[JBrowse](https://jbrowse.org/jb2/) can also be used to visualize annotation tracks.

### Adding functional information to annotations
Begin by extracting maker gene model annotations. The full annotation file will still containt information on mapped transcripts and repeat regions, which are not needed anymore.
```
cat S_20_00703_all_round3_25v23.gff3 | awk '{ if ($2 == "maker") print $0}' > 703.gene_models.gff3
```
Next we rename the bulky maker IDs using maker built-in scripts.
```
#Use maker built in script to create a map for the new names
maker_map_ids --prefix MOBELS_S_20_00703- --justify 5  703.gene.models.gff3 > 703.gene.models.name.map
#Make copies of the files to rename, just to be safe
cp 703.gene_models.gff3 703.gene_models.rename.gff3
cp 703.all.maker.transcripts.fasta 703.all.maker.transcripts.rename.fasta
cp 703.all.maker.proteins.fasta 703.all.maker.proteins.rename.fasta
# replace names in GFF files
map_gff_ids 703.gene_models.name.map 703.gene.models.rename.gff3
# replace names in FASTA headers
map_fasta_ids 703.gene_models.name.map 703.all.maker.transcripts.rename.fasta
map_fasta_ids 703.gene_models.name.map 703.all.maker.proteins.rename.fasta
```
Next we need to generate blast and interproscan output, the functional information to be associated with the gene annotations. Here the non-redundant and curated [Swiss-Prot/UniProt](https://www.uniprot.org/help/downloads) library is used as a blast database.
```
#Get blast results
blastp -query 703.all.maker.proteins.rename.fasta -db uniprot_sprot.fasta -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp
#Get protein domains using interproscan
interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i 703.all.maker.proteins.rename.fasta -o output.iprscan
```
Now we can add the functional information to the annotation file
```
maker_functional_gff uniprot_sprot.fasta output.blastp 703.gene.models.rename.gff3 > 703.gene.models.rename.putative_function.gff
maker_functional_fasta uniprot_sprot.fasta output.blastp 703.all.maker.proteins.rename.fasta > 703.all.maker.proteins.renamed.putative_function.fasta
maker_functional_fasta uniprot_sprot.fasta output.blastp 703.all.maker.transcripts.rename.fasta > 703.all.maker.transcripts.rename.putative_function.fasta
```
# DONE :)

# Citations

```
Alonge, M., Soyk, S., Ramankrishnan, S., Wang, X., Goodwin, S., Sedlazeck, F. J., Lippman, Z. B., Schatz, M. 2019. RaGOO: fast and accurate reference-guided scaffolding of draft genomes. Genome Biology, 20, 224
Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
Benson, G. (1999). Tandem repeats finder: a program to analyse DNA sequences. Nucleic Acids Research, 27, 573-580.
Boa, Z., Eddy, S. R. (2002). Automated de novo identification of repeat sequence families in sequenced genomes. Genome Research, 12, 1269-1276.
Bolger, A. M., Lohse, M., & Usadel, B. (2015). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 30, 2114-20.
Bushmanova, E., Antipov, D., Lapidus, A., Prjibelski, A. D. (2019). ranSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience, 8, giz100.
Bushnell, Brian. (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. https://sourceforge.net/projects/bbmap/.
Dudchenko, O., Batra, S.S, Omer, A.D., Nyquist, S.K, Hoeger, M., Durand, N.C., Shamim, M.S., Machol, I., Lander, E.S., Aiden, A.P., Aiden, E.L 2017. De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds. Science, 356, 92-95.
Dudchenko, O., Shamim, M.S., Batra, S.S., Durand, N.C., Musial, N.T., Mostofa, R., Pham, M., St. Hilaire, B.G., Yao, W., Stamenova, E., Hoeger, M., Nyquist, S.K., Korchina, V., Pletch, K., Flanagan, J.P., Tomaszewicz, A., McAloose, D., Estrada, C.P., Novak, B.J., Omer, A.D., Aiden, E.L. 2018. The Juicebox Assembly Tools module facilitates de novo assembly of mammalian genomes with chromosome-length scaffolds for under $1000. Biorxiv, https://doi.org/10.1101/254797.
Camacho, C., Coulouris, G., Avagyuan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T. (2009). BLAST+:architecture and applications. BMC Bioinformatics, 10, 421.
Ellinghaus, D., Kurtz, S., Willhoeft, U. (2008). LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics, 9, 18.
Ewels, P., Magnusson, M., Lundin, S., Käller, M., (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32, 3047-3048.
Flynn, J. M., Hubley, R., Goubert, C., Rosen, J., Clark, A. G., Feschotte, C., Smit, A. F. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences, 117, 9451-9457.
Fu, L., Niu, B., Zhu, Z., Wu. S., Li, W. (2012). CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, 28, 3150-3152.
Guan, D., McCarthy, S. A., Ning, Z., Wang, G., Wang, Y., Durbin, R. (2021). Efficient iterative Hi-C scaffolder based on N-best neighbors. BMC Bioinformatics, 22, 569.
Holt, C., Yandell, M. 2011. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics, 12, 491.
Katoh, K, Standley, D. M. (2013). MAFFT Multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution, 30, 772-780.
Kolmogorov, M., Yuan, J., Lin, Y, Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37, 540-546.
Korf, I. (2004). Gene finding in novel genomes. BMC Bioinformatics, 5, 59.
Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9.
Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9.
Li, H., Durbin, R. (2009). Fast and accurate short read alignment with Burrow-Wheeler transform. Bioinformatics, 25, 1754-1760.
Nattestad, M., Schatz, M.C. (2016) Assemblytics: a web analytics tool for the detection of variants from an assembly. Bioinformatics, 32, 3021-3023.
Ou, S, Jiang, N. (2018). LTR_retriever: A highly accurate and sensitive program for identification of long terminal repeat retrotransposons. Plant Physiology, 176, 1410-1422.
Price, A. L., Jones, N. C., Pevzner, P. A. (2005). De novo identification of repeat families in large genomes. Bioinformatics, 21, i351-i358.
Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUCSO: Assessing genome assembly and annotation completeness: In: Kollmar M. (ed.) Gene Prediction. Methods in Molecular Biology. Humana New York. 
Shen, W., Le, S., Li, Y., Hu, F. 2016. SeqKit: A cross-platform and ultrafast toolkit for fasta/q file manipulation. PLoS ONE, 11, e0163962.
Smit, A. F. A., Hubley, R., Green, P. (2013-2015). RepeatMasker Open-4.0. http://www.repeatmasker.org.
Stanke, M., Morgenstern, B. (2005). AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Research, 33, W465-467.
Strinegger, M, Söding, J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive datasets. Nature Biotechnology, 35, 1026-1028.
Walker, B. J., Abeel, T., Shea, T., Priest, M., Abouelliel, A., Sakthikumar, S., Cuomo, C. A., Zeng, Q., Wortman, J., Young, S. K., Earl, A. (2014). Pilon: An integrated tool for comprehensive microbial variant detection and genome assembly improvement. PLoS ONE, 9, e112963.
Wheeler, T.J. (2009). Large-scale neighbor-joining with NINJA. In S.L. Salzberg and T. Warnow (Eds.), Proceedings of the 9th Workshop on Algorithms in Bioinformatics. WABI 2009, pp. 375-389. Springer, Berlin.
```
