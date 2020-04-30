#!/bin/bash
#create -n common bwa samtools strelka fastqc
source activate common
export PATH=~/dev/bbtools/bbmap/:$PATH
# fastqc
if [[ $(ls -1 qc/*fastqc.html 2>/dev/null | wc -l) -eq 0 ]] ; then
	fastqc -t $(nproc) -o qc fastq/*.fastq.gz
fi
mkdir -p gridss bams qc variants
# escherichia_coli_o127_h27_str_c43_90.c43_90_irspv01.dna.toplevel.fa 
# Escherichia_coli.HUSEC2011CHR1.dna.toplevel.fa
for ref in Escherichia_coli_bl21_de3_.ASM956v1.dna.toplevel.fa C43DE3_Genome.fasta BL21DE3_Genome.fasta ; do
	prefix=$ref
	prefix=${prefix/_Genome.fasta/}
	prefix=${prefix/escherichia_coli./}
	prefix=${prefix/Escherichia_coli_/}
	prefix=${prefix/.dna.toplevel.fa/}
	echo $prefix
	for FQ1 in fastq/*_R1_001.fastq.gz ; do
		sample=$(echo $(basename $FQ1) | cut -b 3-4)
		FQ2=${FQ1/R1/R2}
		# BWA
		bam=bams/${prefix}_${sample}.bam
		if [[ ! -f $bam ]] ; then
			bwa mem -t $(nproc) $ref $FQ1 $FQ2 | samtools sort -@ $(nproc) - -o $bam
		fi
		if [[ ! -f $bam.bai ]] ; then
			samtools index $bam &
		fi
		# BBMap
		FQ2=${FQ1/R1/R2}
		FQS=_bbduk_adapter_trim.fastq.gz
		if [[ ! -f $FQ1$FQS ]] ; then
			bbduk.sh -Xmx1g in1=$FQ1 in2=$FQ2 out1=$FQ1$FQS out2=$FQ2$FQS ref=~/dev/bbtools/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
		fi
		bam=bams/bb_${prefix}_${sample}.bam
		if [[ ! -f $bam ]] ; then
			bbmap.sh -Xmx4g ref=$ref in=$FQ1$FQS in2=$FQ2$FQS sam=1.3 out=$bam.unsorted.bam && \
			samtools sort -@ $(nproc) $bam.unsorted.bam -o $bam && \
			samtools index $bam && \
			rm $bam.unsorted.bam
		fi
		# QC
		if [[ ! -f qc/${prefix}_${sample}.alignment_summary_metrics ]] ; then
			picard CollectMultipleMetrics I=$bam R=$ref O=qc/${prefix}_${sample} \
				PROGRAM=CollectAlignmentSummaryMetrics \
				PROGRAM=CollectBaseDistributionByCycle \
				PROGRAM=CollectInsertSizeMetrics \
				PROGRAM=MeanQualityByCycle \
				PROGRAM=QualityScoreDistribution \
				PROGRAM=CollectGcBiasMetrics \
				PROGRAM=CollectQualityYieldMetrics \
				&
		fi
		if [[ ! -f qc/${prefix}_${sample}.wgs_metics ]] ; then
			picard CollectWgsMetrics I=$bam R=$ref O=qc/${prefix}_${sample}.wgs_metics & 
		fi
	done
	# SNV calling
	if [[ ! -f variants/$prefix.bb.vcf ]] ; then
		callvariants.sh -Xmx8g ref=$ref multisample=t in=$(ls -1 bams/bb_$prefix*.bam | sort | tr '\n' ',') sample=S1,S2,S3,S4,S5,S6 vcf=variants/$prefix.bb.vcf # zhist=$prefix.bb.zhist shist=$prefix.bb.shist
	fi
	# sequenza
	mkdir -p sequenza
	echo source activate sequenza
	
	if [[ ! -f $ref.gc50Base.wig.gz ]] ; then
		echo sequenza-utils gc_wiggle --fasta $ref -o $ref.gc.wig.gz
	fi
	echo sequenza-utils bam2seqz -n bams/${prefix}_S1.bam -t bams/${prefix}_S2.bam --fasta $ref -gc $ref.gc.wig.gz -o sequenza/S2.seqz.gz "&"
	echo sequenza-utils bam2seqz -n bams/${prefix}_S1.bam -t bams/${prefix}_S3.bam --fasta $ref -gc $ref.gc.wig.gz -o sequenza/S3.seqz.gz "&"
	echo sequenza-utils bam2seqz -n bams/${prefix}_S4.bam -t bams/${prefix}_S5.bam --fasta $ref -gc $ref.gc.wig.gz -o sequenza/S5.seqz.gz "&"
	echo sequenza-utils bam2seqz -n bams/${prefix}_S4.bam -t bams/${prefix}_S6.bam --fasta $ref -gc $ref.gc.wig.gz -o sequenza/S6.seqz.gz "&"
	echo wait
	# sequenza
	mkdir -p cnvkit/S23 cnvkit/S56
	echo source activate cnvkit
	echo cnvkit.py batch -m wgs -n bams/${prefix}_S1.bam -f $ref --output-reference cnvkit/S23/S23_ref.cnn -d cnvkit/S23 bams/${prefix}_S2.bam bams/${prefix}_S3.bam "&" 
	echo cnvkit.py batch -m wgs -n bams/${prefix}_S4.bam -f $ref --output-reference cnvkit/S56/S56_ref.cnn -d cnvkit/S56 bams/${prefix}_S5.bam bams/${prefix}_S6.bam  "&"
	echo wait
	source activate common
	if [[ ! -f gridss/${prefix}_gridss.vcf ]] ; then
		~/dev/gridss/scripts/gridss.sh \
			-j ~/dev/gridss/target/gridss-2.9.0-gridss-jar-with-dependencies.jar \
			-w gridss \
			-o gridss/${prefix}_gridss.vcf \
			-a gridss/${prefix}_assembly.bam \
			-r $ref \
			--jvmheap 8g \
			--sanityCheck \
			bams/${prefix}_*.bam
			#--keepTempFiles
		if [[ ! -f gridss/${prefix}_gridss.vcf ]] ; then
			echo "GRIDSS failure"
			exit 1
		fi
	fi
	if [[ ! -f gridss/bb_${prefix}_gridss.vcf ]] ; then
		~/dev/gridss/scripts/gridss.sh \
			-j ~/dev/gridss/target/gridss-2.9.0-gridss-jar-with-dependencies.jar \
			-w gridss \
			-o gridss/bb_${prefix}_gridss.vcf \
			-a gridss/bb_${prefix}_assembly.bam \
			-r $ref \
			--jvmheap 8g \
			--keepTempFiles \
			--sanityCheck \
			bams/bb_${prefix}_*.bam
		if [[ ! -f gridss/bb_${prefix}_gridss.vcf ]] ; then
			echo "GRIDSS failure"
			exit 1
		fi
	fi
done
wait

