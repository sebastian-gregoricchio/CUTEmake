import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import math
import numpy
import pandas as pd
from itertools import chain

# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
        raise ValueError("constraint_to(): Expected a list, got str instead")

    return "({})".format("|".join(values))

# working diirectory
home_dir = os.path.join(config["workflow_configuration"]["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir


## Directories
MACSDIR = "04_MACS3_peaks/"
GOPEAKSDIR = "05_GoPeaks_peaks/"
SUMMARYDIR = "06_Overall_quality_and_info/"

## Global parameters
MAPQ = str(config["bam_features"]["MAPQ_threshold"])
genome_fasta = str(config["genomic_annotations"]["genome_fasta"])
BLACKLIST = str(config["genomic_annotations"]["blacklist"])



# loading the sample table
sample_metadata = pd.read_csv(str(config["workflow_configuration"]["sample_config_table"]),  sep='\t+', engine='python')   # target_id | input_id | broad
#sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'input_id', 'broad'], axis=1, inplace=False)
sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'control_id', 'broad'], axis=1, copy=False)
TARGETNAMES = list(numpy.unique(list(sample_metadata.target_id)))
INPUTNAMES = list(numpy.unique(list(sample_metadata.control_id)))
SAMPLENAMES = list(numpy.unique(TARGETNAMES + INPUTNAMES))


# generation of global wildcard_constraints
wildcard_constraints:
    SAMPLES = constraint_to(SAMPLENAMES),
    TARGET = constraint_to(TARGETNAMES),
    INPUT = constraint_to(INPUTNAMES)


# Correlations and heatmaps outputs
correlation_outputs = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf")) # PCA 1-2
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf")) # PCA 2-3
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf")) # heatamap_spearman
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf")) # hetamap_pearson
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf")) # scatterplot_spearman
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")) # scatterplot_pearson


# Chromosome remove chr_remove_pattern
if (len(config["bam_features"]["remove_other_chromosomes_pattern"]) > 0):
    chr_remove_pattern = '^chrM|^M|'+config["bam_features"]["remove_other_chromosomes_pattern"]
else:
    chr_remove_pattern = '^chrM|^M'

# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================

# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        BAM_shifted_sorted = expand(os.path.join("01_BAM_filtered/shifted", ''.join(["{sample}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])), sample = SAMPLENAMES),
        multiQC_BAM_html = os.path.join(SUMMARYDIR, ''.join(["multiQC_bams/multiQC_report_BAMs.html"])),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        raw_bigWig = expand(os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{sample}_mapq", MAPQ, "_raw.coverage_bs", str(config["normalization"]["bigWig_binSize"]), ".bw"])), sample=SAMPLENAMES),
        normalized_bigWig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", MAPQ, "_RPGC.normalized_bs", str(config["normalization"]["bigWig_binSize"]), ".bw"])), sample=SAMPLENAMES),
        macs_peaks = expand("04_MACS3_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES),
        gopeaks = expand("05_GoPeaks_peaks/{target}.filtered.BAMPE_peaks.bed", target = TARGETNAMES),
        correlation_outputs = correlation_outputs,
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),
        aggregated_qc_macs = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_MACS.peaks.tsv"),
        aggregated_qc_gopeaks = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_GoPeaks.peaks.tsv")
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================


# -------------------------------------------------------------------------------------------

if not (os.path.exists(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))):
    # ----------------------------------------------------------------------------------------
    # Reads alignement
    rule generate_genome_index:
        input:
            genome = ancient(genome_fasta),
        output:
            genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
        threads: 1
        benchmark:
            "benchmarks/generate_genome_index/generate_genome_index---benchmark.txt"
        shell:
            """
            $CONDA_PREFIX/bin/samtools faidx {input.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------


if (eval(str(config["bam_features"]["skip_bam_filtering"])) == False):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, ".bam"]))),
            bam_mapq_only_sorted_toFix = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_toFix.bam"]))),
            bam_mapq_only_sorted_index_toFix = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_toFix.bai"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bai"]))),
            idxstats_file = os.path.join(SUMMARYDIR, "reads_per_chromosome/{SAMPLES}_idxstats_read_per_chromosome.txt"),
            bam_mapq_only_sorted_woMT = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
            bam_mapq_only_sorted_woMT_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"]))
        params:
            sample = "{SAMPLES}",
            MAPQ_threshold = MAPQ,
            chr_remove_pattern = chr_remove_pattern
        threads:
            workflow.cores
        log:
            fixmate_log = "01_BAM_filtered/FixMateInformation_logs/{SAMPLES}_FixMateInformation.log"
        benchmark:
            "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

            $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted_toFix}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_toFix} {output.bam_mapq_only_sorted_index_toFix}

            $CONDA_PREFIX/bin/gatk FixMateInformation \
            --INPUT {output.bam_mapq_only_sorted_toFix} \
            --OUTPUT {output.bam_mapq_only_sorted} \
            --ASSUME_SORTED false \
            --ADD_MATE_CIGAR true \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT &> {log.fixmate_log}

            $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} > {output.idxstats_file}

            printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '{params.chr_remove_pattern}' | xargs ${{CONDA_PREFIX}}/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} > {output.bam_mapq_only_sorted_woMT}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mapq_only_sorted_woMT} > {output.flagstat_filtered}
            """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
            bai_mdup = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"]))
        params:
            sample = "{SAMPLES}"
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        benchmark:
            "benchmarks/bam_link__skip_filtering/bam_link__skip_filtering---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s $BAM_REAL {output.bam_mdup}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """


# ----------------------------------------------------------------------------------------
# Generate ATAC-like shifted bams
rule bam_shifting:
    input:
        BAM = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"]))
    output:
        BAM_shifted_toSort = temp(os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_shifted.ToSort.bam"]))),
        BAM_shifted_sorted = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_shifted_sorted_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        BAM_shifted_sorted_flagstat = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted_flagstat.txt"]))
    params:
        sample = "{SAMPLES}",
        minFragmentLength = str(config["read_shifting"]["minFragmentLength"]),
        maxFragmentLength = str(config["read_shifting"]["maxFragmentLength"])
    threads:
        workflow.cores
    shell:
        """
        printf '\033[1;36m{params.sample}: Shifting reads in BAM...\\n\033[0m'
        alignmentSieve \
        -p {threads} \
        --ATACshift \
        --bam {input.BAM} \
        --outFile {output.BAM_shifted_toSort} \
        --minFragmentLength {params.minFragmentLength} \
        --maxFragmentLength {params.maxFragmentLength}

        printf '\033[1;36m{params.sample}: Sorting shifted BAM...\\n\033[0m'
        samtools sort -@ {threads} {output.BAM_shifted_toSort} -o {output.BAM_shifted_sorted}
        samtools index -@ {threads} -b {output.BAM_shifted_sorted} {output.BAM_shifted_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from shifted BAM...\\n\033[0m'
        samtools flagstat {output.BAM_shifted_sorted} -@ {threads} > {output.BAM_shifted_sorted_flagstat}
        """
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule fastQC_BAMs:
    input:
        BAM = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"]))
    output:
        html = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_fastqc.html"])),
        zip = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join(config["workflow_configuration"]["output_directory"], "02_BAM_fastQC/"),
        sample = "{SAMPLES}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fastQC_BAMs/fastQC_BAMs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
rule multiQC_BAMs:
    input:
        BAM_fastqc_zip = expand(os.path.join("02_BAM_fastQC/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_fastqc.zip"])), sample=SAMPLENAMES),
        macs_peaks = expand("04_MACS3_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES),
        gopeaks = expand("05_GoPeaks_peaks/{target}.filtered.BAMPE_peaks.bed", target = TARGETNAMES)
    output:
        multiqcReportBAM = os.path.join(SUMMARYDIR, ''.join(["multiQC_bams/multiQC_report_BAMs.html"]))
    params:
        fastQC_BAM_reports_dir = "02_BAM_fastQC/",
        BAM_flagstat_dir = "01_BAM_filtered/flagstat/",
        macs_dir = MACSDIR,
        gopeaks_dir = GOPEAKSDIR,
        multiQC_BAM_outdir = os.path.join(config["workflow_configuration"]["output_directory"], SUMMARYDIR, ''.join(["multiQC_bams/"]))
    benchmark:
        "benchmarks/multiQC_BAMs/multiQC_BAMs---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'

        $CONDA_PREFIX/bin/multiqc -f \
        --outdir {params.multiQC_BAM_outdir} \
        -n multiQC_report_BAMs.html \
        --dirs \
        {params.fastQC_BAM_reports_dir} \
        {params.BAM_flagstat_dir} \
        {params.macs_dir} \
        {params.gopeaks_dir}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution:
    input:
        BAM = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"]))
    output:
        plot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{SAMPLES}_fragment_size_distribution.pdf"),
        fragmentSize_RawFragmentLengths = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_RawFragmentLengths.txt"),
        fragmentSize_metrics = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_metrics.txt")
    params:
        build_summary_directory = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log"),
        build_summary_directory_table = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize"),
        sample = "{SAMPLES}",
        plotFormat = "pdf",
        binSize = str(config["quality_controls"]["fragmentSize_window_length"]),
        blacklist = BLACKLIST,
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    log:
        out = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.out"),
        err = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.err")
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fragment_size_distribution/fragment_size_distribution---{SAMPLES}_benchmark.txt"
    priority: -1
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        mkdir -p {params.build_summary_directory}
        mkdir -p {params.build_summary_directory_table}

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        -p {threads} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.maxFragmentLength} \
        --blackListFileName {params.blacklist} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        --table {output.fragmentSize_metrics} \
        --histogram {output.plot} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution_report:
    input:
        plots = expand(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES)
    output:
        replot_script = temp(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/replot_script.R")),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        report_ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf")
    params:
        distribution_plots_pattern = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf"),
        dir = os.path.join(home_dir,""),
        summary_dir = SUMMARYDIR,
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    threads: 1
    log:
      ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/ggplot_replotting.log"),
      pdfcombine = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/pdfcombine.log")
    benchmark:
        "benchmarks/fragment_size_distribution_report/fragment_size_distribution_report---benchmark.txt"
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'
        $CONDA_PREFIX/bin/pdfcombine {params.distribution_plots_pattern} -o {output.report_pdf} -sf &> {log.pdfcombine}


        printf '\033[1;36mReplotting fragmentSizeDistribution reports in R (ggplot version)...\\n\033[0m'
        echo "tb = do.call(rbind, lapply(list.files('{params.dir}{params.summary_dir}fragmentSizeDistribution_plots/table_and_fragmentSize', pattern = 'RawFragmentLengths', full.names = T), function(x)(read.delim(x, h=T, skip=1))))" > {output.replot_script}
        echo "n.samples = length(unique(tb[,3]))" >> {output.replot_script}
        echo "plot = ggplot2::ggplot(data = tb, ggplot2::aes(x = Size, y = Occurrences, color = Sample)) + ggplot2::geom_smooth(method = 'loess', formula = y ~ x, span = 0.05, show.legend = F, se = F, color = 'navyblue', linewidth = 0.5) + ggplot2::xlim(c(1,{params.maxFragmentLength})) + ggplot2::theme_classic() + ggplot2::facet_wrap(~Sample, scale='free', ncol = floor(sqrt(n.samples))) + ggplot2::theme(axis.ticks = ggplot2::element_line(color ='black'), axis.text = ggplot2::element_text(color = 'black'), strip.background = ggplot2::element_blank())" >> {output.replot_script}
        echo "pdf(file = '{params.dir}{output.report_ggplot}', width = floor(sqrt(n.samples)) * 2.7, height = ceiling(n.samples / floor(sqrt(n.samples))) * 1.5)" >> {output.replot_script}
        echo "print(plot)" >> {output.replot_script}
        echo "invisible(dev.off())" >> {output.replot_script}

        $CONDA_PREFIX/bin/Rscript {output.replot_script} &> {log.ggplot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Extract chromosome sizes fro converstion to bigWig
rule extract_chromosome_sizes:
    input:
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt"
    benchmark:
        "benchmarks/extract_chromosome_sizes/extract_chromosome_sizes---benchmark.txt"
    shell:
        """
        printf '\033[1;36mExtracting chromosome sizes from genome index (.fai)...\\n\033[0m'
        cut -f1,2 {input.genome_fai} > {output.chrSizes}
        """
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# PeakCalling on shifted bams (MACS3 callpeak)
rule peakCalling_MACS3:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        input_bam_all = expand(os.path.join("01_BAM_filtered/shifted", ''.join(["{input}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])), input = INPUTNAMES)
    output:
        peaksPE = "04_MACS3_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
    params:
        sample = "{TARGET}",
        macs_version = "macs3",
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        input_suffix = "_mapq"+MAPQ+"_woMT_shifted_sorted.bam",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        macs_qValue_cutoff = config["peak_calling"]["macs_qValue_cutoff"],
        blacklist = BLACKLIST
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = "04_MACS3_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.out",
        err = "04_MACS3_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.err"
    benchmark:
        "benchmarks/peakCalling_MACS3/peakCalling_MACS3---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

        INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)
        CALL_BROAD=$(grep {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

        if [ $CALL_BROAD == "false" ]; then
            BROAD=""
            PEAKEXT=narrowPeak
            SUMMITS='--call-summits'
        else
            BROAD="--broad"
            PEAKEXT=broadPeak
            SUMMITS=''
        fi


        if [ $INPUT_ID == {params.sample} ]; then
            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -f BAMPE \
            -g {params.genomeSize} \
            -q {params.macs_qValue_cutoff} \
            --keep-dup all \
            --outdir 04_MACS3_peaks \
            --name {params.sample}.filtered.BAMPE ${{BROAD}} ${{SUMMITS}} > {log.err} 2> {log.out}
        else
            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/shifted/${{INPUT_ID}}{params.input_suffix} \
            -f BAMPE \
            -g {params.genomeSize} \
            -q {params.macs_qValue_cutoff} \
            --keep-dup all \
            --outdir 04_MACS3_peaks \
            --name {params.sample}.filtered.BAMPE ${{BROAD}} ${{SUMMITS}} > {log.err} 2> {log.out}
        fi

        # add chr to peak files
        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a 04_MACS3_peaks/{params.sample}.filtered.BAMPE_peaks.${{PEAKEXT}} -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > 04_MACS3_peaks/{params.sample}.filtered.BAMPE_peaks_chr.${{PEAKEXT}}
        """


##########################################
###            PEAK CALLING            ###
##########################################
# ----------------------------------------------------------------------------------------
# PeakCalling on shifted bams (GoPeaks callpeak)
rule peakCalling_GoPeaks:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        input_bam_all = expand(os.path.join("01_BAM_filtered/shifted", ''.join(["{input}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])), input = INPUTNAMES),
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt"
    output:
        peaksPE = "05_GoPeaks_peaks/{TARGET}.filtered.BAMPE_peaks.bed"
    params:
        sample = "{TARGET}",
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        input_suffix = "_mapq"+MAPQ+"_woMT_shifted_sorted.bam",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        gopeaks_pValue_cutoff = config["peak_calling"]["GoPeaks_pValue_cutoff"],
        blacklist = BLACKLIST
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = "05_GoPeaks_peaks/logs/{TARGET}_GoPeaks_log.out",
        err = "05_GoPeaks_peaks/logs/{TARGET}_GoPeaks_log.err"
    benchmark:
        "benchmarks/peakCalling_GoPeaks/peakCalling_GoPeaks---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: calling peaks (GoPeaks)...\\n\033[0m'

        INPUT_ID=$(grep {params.sample} {params.sample_config_table} | cut -f 2)
        CALL_BROAD=$(grep {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

        if [ $CALL_BROAD == "false" ]; then
            BROAD=""
            MDIST=1000
        else
            BROAD="--broad"
            MDIST=3000
        fi


        if [ $INPUT_ID == {params.sample} ]; then
            $CONDA_PREFIX/bin/gopeaks \
            --bam {input.target_bam} \
            --chromsize {input.chrSizes} \
            --mdist $MDIST \
            --pval {params.gopeaks_pValue_cutoff} \
            --prefix 05_GoPeaks_peaks/{params.sample}.filtered.BAMPE ${{BROAD}} > {log.err} 2> {log.out}
        else
            $CONDA_PREFIX/bin/gopeaks \
            --bam {input.target_bam} \
            --control 01_BAM_filtered/shifted/${{INPUT_ID}}{params.input_suffix} \
            --chromsize {input.chrSizes} \
            --mdist $MDIST \
            --pval {params.gopeaks_pValue_cutoff} \
            --prefix 05_GoPeaks_peaks/{params.sample}.filtered.BAMPE ${{BROAD}} > {log.err} 2> {log.out}
        fi

        # add chr to peak files
        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a 05_GoPeaks_peaks/{params.sample}.filtered.BAMPE_peaks.bed -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > 05_GoPeaks_peaks/{params.sample}.filtered.BAMPE_chr_peaks.bed
        """
# ----------------------------------------------------------------------------------------


##########################################
###      Normalization (bigWigs)       ###
##########################################

rule raw_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"]))
    output:
        raw_bigWig = os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{SAMPLES}_mapq", MAPQ, "_raw.coverage_bs", str(config["normalization"]["bigWig_binSize"]), ".bw"]))
    params:
        sample = "{SAMPLES}",
        blacklist = BLACKLIST,
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        ignore_for_normalization = config["genomic_annotations"]["ignore_for_normalization"],
        bw_binSize = config["normalization"]["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLES}_raw_bamCoverage_log.out",
        err = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLES}_raw_bamCoverage_log.err"
    benchmark:
        "benchmarks/raw_bigWig/raw_bigWig---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating raw coverage bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/raw_coverage/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.raw_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing None \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --extendReads \
        -p {threads} > {log.out} 2> {log.err}
        """


rule RPGC_normalized_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"]))
    output:
        RPGC_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{SAMPLES}_mapq", MAPQ, "_RPGC.normalized_bs", str(config["normalization"]["bigWig_binSize"]), ".bw"]))
    params:
        sample = "{SAMPLES}",
        blacklist = BLACKLIST,
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        ignore_for_normalization = config["genomic_annotations"]["ignore_for_normalization"],
        bw_binSize = config["normalization"]["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLES}_RPGC.normalized_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLES}_RPGC.normalized_log.err"
    benchmark:
        "benchmarks/RPGC_normalized_bigWig/RPGC_normalized_bigWig---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC-normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.RPGC_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --extendReads \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------









##########################################
###          Quality controls          ###
##########################################

# ----------------------------------------------------------------------------------------
# Generation of matrix scores
rule multiBigwigSummary:
    input:
        norm_bw = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", MAPQ, "_RPGC.normalized_bs", str(config["normalization"]["bigWig_binSize"]), ".bw"])), sample=SAMPLENAMES)
    output:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz"),
        rawCounts = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.txt")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "/Sample_comparisons/"])),
        labels = ' '.join(SAMPLENAMES),
        window = config["quality_controls"]["multiBigwigSummary_binning_window_size"],
        blacklist = BLACKLIST
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/multiBigwigSummary/multiBigwigSummary---benchmark.txt"
    shell:
        """
        printf '\033[1;36mComparing the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        $CONDA_PREFIX/bin/multiBigwigSummary bins \
        -p {threads} \
        -b {input.norm_bw} \
        --labels {params.labels} \
        --binSize {params.window} \
        --blackListFileName {params.blacklist} \
        -o {output.matrix} \
        --outRawCounts {output.rawCounts}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule PCA_and_samples_correlation:
    input:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz")
    output:
        PCA_OneTwo = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf"),
        PCA_TwoThree = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf"),
        PCA_tab = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA_on_BigWigs_wholeGenome_values.txt"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"),
        matrix_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_matrix_on_BigWigs_wholeGenome_spearmanMethod.txt"),
        matrix_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_matrix_on_BigWigs_wholeGenome_pearsonMethod.txt"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")
    params:
        make_directory = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/"),
        heatmap_color = config["quality_controls"]["correlation_heatmap_color"]
    benchmark:
        "benchmarks/PCA_and_samples_correlation/PCA_and_samples_correlation---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        printf '\033[1;36m    - plotting PCAs (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_OneTwo} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf' --outFileNameData {output.PCA_tab}
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_TwoThree} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf' --PCs 2 3


        printf '\033[1;36m    - plotting Spearman correlation heatmap (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        --outFileCorMatrix {output.matrix_spearman} \
        -o {output.hetamap_spearman}

        printf '\033[1;36m    - plotting Pearson correlation heatmap (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        --outFileCorMatrix {output.matrix_pearson} \
        -o {output.hetamap_pearson}



        printf '\033[1;36m    - plotting Spearman correlation scatterplot (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_spearman}

        printf '\033[1;36m    - plotting Pearson correlation scatterplot (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_pearson}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of Lorenz curves
rule Lorenz_curve:
    input:
        BAM_sorted = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_sorted_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"]))
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{SAMPLES}_Lorenz_curve_deeptools.plotFingreprint.pdf"),
        lorenz_metrics = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_metrics/{SAMPLES}_Lorenz_quality.metrics_deeptools.plotFingreprint.txt"),
        lorenz_counts = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_counts/{SAMPLES}_Lorenz_raw.counts_deeptools.plotFingreprint.txt")
    params:
        #all_bams = ' '.join(expand(os.path.join(''.join(["01_BAM_filtered/{SAMPLES}_mapq", MAPQ ,"_sorted_woMT.bam"])), sample=SAMPLENAMES)),
        #labels = ' '.join(SAMPLENAMES),
        labels = str("{SAMPLES}"),
        blacklist = BLACKLIST,
        binSize = config["quality_controls"]["plotFingerprint"]["binSize"],
        sampledRegions = config["quality_controls"]["plotFingerprint"]["sampledRegions"],
        extra_params = config["quality_controls"]["plotFingerprint"]["extra_parameters"]
    threads:
        max(math.floor((workflow.cores-1)/2), 1)
    log:
        out = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/{SAMPLES}_deeptools_plotFingreprint.log")
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve---{SAMPLES}_benchmark.txt"
    priority: -10
    shell:
        """
        printf '\033[1;36m{params.labels}: plotting Lorenz curves-Fingerprint...\\n\033[0m'

        $CONDA_PREFIX/bin/plotFingerprint \
        --bamfiles {input.BAM_sorted} \
        --plotFile {output.lorenz_plot} \
        --labels {params.labels} \
        --blackListFileName {params.blacklist} \
        --binSize {params.binSize} \
        --numberOfSamples {params.sampledRegions} \
        --outQualityMetrics {output.lorenz_metrics} \
        --outRawCounts {output.lorenz_counts} \
        -p {threads} {params.extra_params} &> {log.out}
        """


rule Lorenz_curve_merge_plots:
    input:
        lorenz_plots = expand(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{sample}_Lorenz_curve_deeptools.plotFingreprint.pdf"), sample=SAMPLENAMES)
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples_combined.pdf"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf")
    params:
        lorenz_plots_pattern = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/*_Lorenz_curve_deeptools.plotFingreprint.pdf"),
        dir = os.path.join(home_dir,""),
        summary_dir = SUMMARYDIR,
    threads: 1
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve_merge_plots---benchmark.txt"
    priority: -5
    shell:
        """
        printf '\033[1;36mCombine Lorenz curves-Fingerprint for all samples...\\n\033[0m'
        $CONDA_PREFIX/bin/pdfcombine {params.lorenz_plots_pattern} -o {output.lorenz_plot} -sf

        printf '\033[1;36mMake combined Lorenz curves-Fingerprint plot...\\n\033[0m'
        $CONDA_PREFIX/bin/Rscript \
        -e "require(dplyr)" \
        -e "tables = list.files(path = '{params.dir}06_Overall_quality_and_info/LorenzCurve_plotFingreprint/lorenz_counts', pattern = '.plotFingreprint.txt', full.names = T)" \
        -e "combined_table = data.frame()" \
        -e "for (i in 1:length(tables)) (combined_table = rbind(combined_table, dplyr::mutate(read.delim(tables[i], skip = 2, h=F), sample = gsub('_Lorenz_raw[.]counts_deeptools[.]plotFingreprint[.]txt','',basename(tables[i]))) %>% dplyr::rename(counts = V1) %>% dplyr::arrange(counts) %>% dplyr::mutate(cumulative_sum = cumsum(counts), rank = (1:nrow(.))/nrow(.)) %>% dplyr::mutate(cumulative_sum = cumulative_sum/max(cumulative_sum))))" \
        -e "pdf('{params.dir}{output.lorenz_plot_ggplot}', width = 8, height = 6.5)" \
        -e "ggplot2::ggplot(data = combined_table, ggplot2::aes(x = rank, y = cumulative_sum, color = sample)) + ggplot2::geom_line() + ggplot2::ggtitle('Fingerprints (Lorenz curves) all samples') + ggplot2::xlim(c(0,1)) + ggplot2::xlab('Normalized rank') + ggplot2::ylab('Fraction with reference to the bin with highest coverage') + ggplot2::theme_classic() + ggplot2::theme(axis.text = ggplot2::element_text(color = 'black'), axis.ticks = ggplot2::element_line(color = 'black'))" \
        -e "invisible(dev.off())"
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
## Compute FRiP for MACS peaks
rule MACS_peak_QC:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        peaks = "04_MACS3_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
    output:
        qc = temp(os.path.join(SUMMARYDIR, "peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc_MACS.txt"))
    params:
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        peak_prefix = "04_MACS3_peaks/{TARGET}.filtered.BAMPE_peaks",
        blacklist = BLACKLIST,
        target = "{TARGET}",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    benchmark:
        "benchmarks/MACS_peak_QC/MACS_peak_QC---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: computing peak stats (MACS)...\\n\033[0m'

        # define peak file
        CALL_BROAD=$(grep {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

        if [ $CALL_BROAD == "false" ]; then
            CALLING_MODE="narrow"
            PEAK="{params.peak_prefix}.narrowPeak"
        else
            CALLING_MODE="broad"
            PEAK="{params.peak_prefix}.broadPeak"
        fi

        # get the number of peaks
        peak_count=$(wc -l < $PEAK)

        # get the number of mapped reads
        mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

        # calculate Fraction of Reads In Peaks
        frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

        # compute peak genome coverage
        peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
        genome_size={params.genomeSize}
        genomecov=$(bc -l <<< "$peak_len/$genome_size")

        # rounding fractions
        genomecov_round=$(printf "%.5f\n" "$genomecov")
        frip_round=$(printf "%.3f\n" "$frip")

        # write peak-based QC metrics to output file
        printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}
        """

rule aggregate_FRiP_MACS:
    input:
        qc = expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_MACS.txt"), target = TARGETNAMES)
    output:
        aggregated_qc = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_MACS.peaks.tsv")
    params:
        all_qc = ' '.join(expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_MACS.txt"), target = TARGETNAMES))
    threads: 1
    benchmark:
        "benchmarks/aggregate_FRiP_MACS/aggregate_FRiP_MACS---allTargets_benchmark.txt"
    shell:
        """
        # print header of FRiP report
        printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
        cat {params.all_qc} >> {output.aggregated_qc}
        """

# ------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
## Compute FRiP for MACS peaks
rule GoPeaks_peak_QC:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        peaks = "05_GoPeaks_peaks/{TARGET}.filtered.BAMPE_peaks.bed"
    output:
        qc = temp(os.path.join(SUMMARYDIR, "peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc_GoPeaks.txt"))
    params:
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        blacklist = BLACKLIST,
        target = "{TARGET}",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    benchmark:
        "benchmarks/GoPeaks_peak_QC/GoPeaks_peak_QC---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: computing peak stats (GoPeaks)...\\n\033[0m'

        # define peak file
        CALL_BROAD=$(grep {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

        if [ $CALL_BROAD == "false" ]; then
            CALLING_MODE="narrow"
        else
            CALLING_MODE="broad"
        fi

        PEAK={input.peaks}

        # get the number of peaks
        peak_count=$(wc -l < $PEAK)

        # get the number of mapped reads
        mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

        # calculate Fraction of Reads In Peaks
        frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

        # compute peak genome coverage
        peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
        genome_size={params.genomeSize}
        genomecov=$(bc -l <<< "$peak_len/$genome_size")

        # rounding fractions
        genomecov_round=$(printf "%.5f\n" "$genomecov")
        frip_round=$(printf "%.3f\n" "$frip")

        # write peak-based QC metrics to output file
        printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}
        """

rule aggregate_FRiP_GoPeaks:
    input:
        qc = expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_GoPeaks.txt"), target = TARGETNAMES)
    output:
        aggregated_qc = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_GoPeaks.peaks.tsv")
    params:
        all_qc = ' '.join(expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_GoPeaks.txt"), target = TARGETNAMES))
    threads: 1
    benchmark:
        "benchmarks/aggregate_FRiP_GoPeaks/aggregate_FRiP_GoPeaks---allTargets_benchmark.txt"
    shell:
        """
        # print header of FRiP report
        printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
        cat {params.all_qc} >> {output.aggregated_qc}
        """

# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
#                                 END PIPELINE                                 #
# ------------------------------------------------------------------------------
