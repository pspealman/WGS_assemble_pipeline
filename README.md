## WGS_assemble_pipeline
Pipeline for de novo WGS assembly beginning with a hybrid (long and short-read) method. 

The intent of this script is not to have a complete and annotated genome by the end, it is intended to produce high quality WGS contigs (as described by the NCBI [here](https://www.ncbi.nlm.nih.gov/genbank/wgs/)) which are the first step in assembling a genome. By the end of this pipeline you should have a FASTA file of sufficient quality for a submission. 

This pipeline relies on several tools:

1. [Minimap2](https://github.com/lh3/minimap2) - long read aligner (v2.22)
2. [MaSuRCA](https://github.com/alekseyzimin/masurca) - de novo hybrid assembler (v4.0.9)
3. [Ragtag](https://github.com/malonge/RagTag) - Genome assembly toolkit (v2.1.0)
4. [Liftoff](https://github.com/agshumate/Liftoff) - Genome annotation mapping (lift over with no chain hassle)
5. [BUSCO](https://busco.ezlab.org/) - Genome assembly QC (v5.3.0)
6. [mugio.py](https://github.com/pspealman/mugio) - Long read toolkit (v1.7)

The flow of information is: 
```
[Long + short reads] -> (MaSuRCA) -> [Raw Scaffolds] -> (Ragtag) -> [Corrected Scaffolds], [QC]
  [Corrected Scaffolds] -> (Ragtag) -> [QC]
  [Corrected Scaffolds] -> (BUSCO) -> [QC]
  [Corrected Scaffolds] -> (mugio.py) -> [Relabelled Scaffolds]
```
### Setting default parameters:
Here we set declare the locations for both MaSuRCA and mugio.py, it is assumed the other programs can be invoked in your environment by name if correctly intalled.
```
masurca=<path_to_masurca>/MaSuRCA-4.0.9/bin/masurca
mugio=<path_to_mugio>/mugio.py

ref_fa=<path_to_reference_fasta>/GCF_000146045.2_R64_genomic_GAP1.fna
gff=<path_to_reference_gff>/GCF_000146045.2_R64_genomic_GAP1_edit.gff
```
### Example script
One example script would be:
```
# set variables for a single genome
strain=DGY1734 #for the name in ont and outpath
	read_1=<path_to_illumina_fastq>/n01.fastq.gz
	read_2=<path_to_illumina_fastq>/n02.fastq.gz
	ont=<path_to_ont_fastq>${strain}gDNA.fastq
	outpath=${strain}_gDNA/
  # masurca always outputs the primary scaffold file in this location
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
  # ragtag always outputs the corrected scaffolds in this location
  scf_fix=${outpath}/ragtag/ragtag.correct.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
  #Masurca builds assemblies
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
  #
	#Some assemblies will have errors that cross the wrong centromere, telomeres, or transposons
  # ragtag correct will identify these and break them back down to smaller contigs
	ragtag.py correct -t 4 ${ref_fa} ${scf_pri} \
		-R ${ont} \
		-T ont --remove-small -u -w \
		-o ${outpath}/ragtag
	#
  # For submission mitochondria and chromosomes need to be identified
  # Liftoff will generate a new GFF using a minimap2 enabled lookup from the reference genome
	liftoff ${scf_fix} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/liftoff.txt
	#
  # mugio then parses the liftoff gff and corrects the ragtag corrected scaffolds
  # it will relabel mitochondrial sequences and chromosomes with 90% chromosome match
  # if a chromosome only has a single contig it will be correctly designated 
  # contig names will be abbreviated for submission (<= 50 characters)
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/liftoff.txt \
		-f ${scf_fix} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_submission
	#
  # 
  # busco will produce a score of deeply conserved CDS present within the assembly
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
  # Ragtag's asmstats will calculate the number of contigs (n), the total bases, and N50
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt
```
### Results

QC should be evaluated for each assembly as the complexity of the CNV structures will dictate the utility of the results. However, a good hueristic would be a BUSCO score > 90%, median depth > 15x, and an N50 > genome_size/300 (adapted from [Jung et al. 2019](https://doi.org/10.1016/j.tplants.2019.05.003) and [Jung et al. 2020](https://doi.org/10.1371/journal.pcbi.1008325)) 

Final scaffolds will be available in the `${outpath}/${strain}_submission` location

