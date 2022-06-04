#conda install -c bioconda ragtag

module load minimap2/2.22
module load busco/5.3.0

masurca=/scratch/ps163/masurca/MaSuRCA-4.0.9/bin/masurca
chromosome_scaffolder=/scratch/ps163/masurca/MaSuRCA-4.0.9/bin/chromosome_scaffolder.sh
mugio=/scratch/ps163/Project_Grace/mugio.py

ref_fa=/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna
gff=/scratch/cgsb/gresham/pieter/genome_annotations/Inhouse_GAP1/GCF_000146045.2_R64_genomic_GAP1_edit.gff

strain=DGY1734
nid=79
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/${strain}_gDNA/pass/${strain}_gDNA.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	#
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
	#
	ragtag.py correct -t 4 ${ref_fa} ${scf_pri} \
		-R ${ont} \
		-T ont --remove-small -u -w \
		-o ${outpath}/ragtag
		
	#
	liftoff ${outpath}/ragtag/ragtag.correct.fasta \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/liftoff.txt \
		-f ${outpath}/ragtag/ragtag.correct.fasta \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_submission
	#
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt
#
strain=DGY1736
nid=82
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/${strain}_gDNA/pass/${strain}_gDNA.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
		
	ragtag.py correct -t 4 ${ref_fa} ${scf_pri} \
		-R ${ont} \
		-T ont --remove-small -f 10000 -u -w \
		-o ${outpath}/ragtag
		
	
	liftoff ${outpath}/ragtag/ragtag.correct.fasta \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/liftoff.txt
	
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/liftoff.txt \
		-f ${outpath}/ragtag/ragtag.correct.fasta \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_submission
	#
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt
#
strain=DGY1740
nid=52
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/${strain}_gDNA/pass/${strain}_gDNA.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
	#
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt
#
strain=DGY1747
nid=49
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/${strain}_gDNA/pass/${strain}_gDNA.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
	#
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt

strain=DGY1751
nid=78
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/${strain}_gDNA/pass/${strain}_gDNA.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
	
	busco -i ${outpath}/ragtag/ragtag.correct.fasta -l saccharomycetes_odb10 -o ${outpath}/busco_${strain} -m genome > ${outpath}/busco_${strain}_results.txt
	ragtag.py asmstats ${outpath}/ragtag/ragtag.correct.fasta > ${outpath}/asmstats_${strain}_results.txt

strain=DGY1657
nid=40
	read_1=/scratch/cgsb/gencore/out/Gresham/2016-10-26_HG7CVAFXX/merged/HG7CVAFXX_n01_mini02_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2016-10-26_HG7CVAFXX/merged/HG7CVAFXX_n02_mini02_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/erisapfel_nanopore/fastq/${strain}.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
	
strain=DGY1728
nid=50
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/erisapfel_nanopore/fastq/${strain}.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission

strain=DGY1744
nid=83
	read_1=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n01_mini02_partii_${nid}.fastq.gz
	read_2=/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged//HKFYTBGX2_n02_mini02_partii_${nid}.fastq.gz
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/erisapfel_nanopore/fastq/${in_name}.fastq
	ont=/scratch/cgsb/gresham/LABSHARE/Data/nanopore/erisapfel_nanopore/fastq/${strain}.fastq
	outpath=/scratch/ps163/Project_Grace/masurca/${strain}_gDNA/
	scf_alt=${outpath}/CA.mr.49.17.15.0.02/alternative.genome.scf.fasta
	scf_pri=${outpath}/CA.mr.49.17.15.0.02/primary.genome.scf.fasta
	#
	mkdir -p ${outpath}
	cd ${outpath}
	#
	${masurca} \
		-t 32 \
		-i ${read_1},${read_2} \
		-r ${ont} \
		-o ${outpath}
		
	liftoff ${scf_alt} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_alt_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_alt_liftoff.txt \
		-f ${scf_alt} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_alt_submission
	
	liftoff ${scf_pri} \
		${ref_fa} -g ${gff} \
		-copies -sc 0.75 \
		-o ${outpath}/scf_pri_liftoff.txt
	#
	python ${mugio} --assign_from_lift_off \
		-gff ${gff} \
		-lift_gff ${outpath}/scf_pri_liftoff.txt \
		-f ${scf_pri} \
		--minimum 90 -organism "Saccharomyces cerevisiae" -isolate ${strain} -mito NC_001224.1 \
		-o ${outpath}/${strain}_scf_pri_submission
