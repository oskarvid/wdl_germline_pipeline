workflow GermlineVarCall {
# Input files
  Array[File] known_indels_sites_indices
  Array[File] known_indels_sites_VCFs
  Array[File] scatter_interval_list
  Array[File] fastqfiles
  File input_fastq1
  File input_fastq2

# Reference files 
  File ref_fasta_index
  File ref_fasta
  File fasta_bwt
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  
# Reference VCF files
  File v1000g_vcf
  File omni_vcf
  File dbsnp_vcf
  File hapmap_vcf
  File mills_vcf

# Reference VCF file indexes
  File v1000g_vcf_index
  File omni_vcf_index
  File dbsnp_vcf_index
  File hapmap_vcf_index
  File mills_vcf_index

# Scatter Gather shard group configuration file
  File groups

# Jar files  
  File picard
  File gatk3

# String names  
  String Base_Name
  String final_gvcf_name
  String recalibrated_bam_basename = Base_Name + ".aligned.duplicates_marked.recalibrated"
  String outputfolder = "/wdl_pipeline/"
  String unmapped_basename = "unmapped_bam"

call CreateSequenceGroupingTSV {
  input:
    Groups = groups,
}

scatter (fastq in fastqfiles) {
  call FastqToSam {
   input:
      PICARD = picard,
      Input_Fastq1 = fastq,
     Unmapped_Basename = unmapped_basename,
  }

  call BwaMem {
    input:
      fasta_bwt = fasta_bwt,
      ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_bwt = ref_bwt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      Input_Fastq1 = fastq,
      Base_Name = Base_Name + ".bwa",
  }

  call MergeBamAlignment {
    input:
      PICARD = picard,
      ref_fasta_index = ref_fasta_index,
      Unmapped_Bam = FastqToSam.outputbam,
      Aligned_Bam = BwaMem.outputfile,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      Output_Bam_Basename = unmapped_basename,
  }
}
call MarkDup {
  input:
    PICARD = picard,
    Base_Name = Base_Name + ".markdup.sortsam.bwa",
    Input_File = MergeBamAlignment.output_bam,
}
    
  # Perform Base Quality Score Recalibration (BQSR) and call variants on the sorted BAM in parallel
  scatter (subgroup in scatter_interval_list) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        GATK3 = gatk3,
        Input_Bam = MarkDup.MarkDupOutputBam,
        Input_Bam_Index = MarkDup.MarkDupOutputBai,
        Recalibration_Report_Filename = Base_Name + ".recal_data.grp",
        Sequence_Group_Interval = subgroup,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        v1000g_vcf = v1000g_vcf,
        mills_vcf = mills_vcf,
        v1000g_vcf_index = v1000g_vcf_index,
        mills_vcf_index = mills_vcf_index,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
    # Apply the recalibration model by interval
    call PrintReads {
      input:
        GATK3 = gatk3,
        Input_Bam = MarkDup.MarkDupOutputBam,
        Input_Bam_Index = MarkDup.MarkDupOutputBai,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        Recalibration_Report = BaseRecalibrator.Recalibration_Report,
        Output_Bam_Basename = recalibrated_bam_basename,
        Sequence_Group_Interval = subgroup,
    }
    
    # Generate GVCFs
    call HaplotypeCaller {
      input:
        GATK3 = gatk3,
        Input_Bam = PrintReads.recalibrated_bam,
        Input_Bam_Index = PrintReads.recalibrated_bam_index,
        Sequence_Group_Interval = subgroup,
        Gvcf_Basename = Base_Name + ".haplotypecaller.bqsr.baserecal.markdup.sortsam.bwa",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      Input_Bqsr_Reports = BaseRecalibrator.Recalibration_Report,
      Output_Report_Filename = Base_Name + ".recal_data.grp",
      GATK3 = gatk3,
  }

  call GatherBamFiles {
    input:
      PICARD = picard,
      Input_Bams = PrintReads.recalibrated_bam,
      Output_Bam_Basename = Base_Name + ".bqsr.baserecal.markdup.sortsam.bwa",
  }

  # Combine GVCFs into a single sample GVCF file
  call GatherVCFs {
    input:
      PICARD = picard,
      Input_Vcfs = HaplotypeCaller.output_gvcf,
      Input_Vcfs_Indexes = HaplotypeCaller.output_gvcf_index,
      Output_Vcf_Name = final_gvcf_name,
  }

  call GenotypeGVCFs {
    input:
      GATK3 = gatk3,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf = HaplotypeCaller.output_gvcf,
      Input_Vcf_Index = HaplotypeCaller.output_gvcf_index,
      Output_Name = final_gvcf_name + ".genotypegvcf",
  }

  call VariantRecalibratorSNP {
    input:
      GATK3 = gatk3,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf_Index = GenotypeGVCFs.output_vcf_index,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      v1000g_vcf = v1000g_vcf,
      v1000g_vcf_index = v1000g_vcf_index,
      omni_vcf = omni_vcf,
      omni_vcf_index = omni_vcf_index,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      hapmap_vcf = hapmap_vcf,
      hapmap_vcf_index = hapmap_vcf_index,
      Mode = "SNP",
      Output_Vcf_Name = final_gvcf_name + ".varRec.SNP",
  }

  call VariantRecalibratorINDEL {
    input:
      GATK3 = gatk3,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf_Index = GenotypeGVCFs.output_vcf_index,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      mills_vcf = mills_vcf,
      mills_vcf_index = mills_vcf_index,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      Mode = "INDEL",
      Output_Vcf_Name = final_gvcf_name + ".varRec.INDEL",
  }

  call ApplyRecalibrationSNP {
    input:
      GATK3 = gatk3,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      TranchesFile = VariantRecalibratorSNP.tranchesFile,
      RecalFile = VariantRecalibratorSNP.recalFile,
      Mode = "SNP",
      Output_Vcf_Name = final_gvcf_name + ".applyRec.SNP",
  }

  call ApplyRecalibrationINDEL {
    input:
      GATK3 = gatk3,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      TranchesFile = VariantRecalibratorINDEL.tranchesFile,
      RecalFile = VariantRecalibratorINDEL.recalFile,
      Mode = "INDEL",
      Output_Vcf_Name = final_gvcf_name + ".applyRec.INDEL",
  }

  # Outputs that will be retained when execution is complete
  
  output {
    MarkDup.*
    GatherBqsrReports.*
    GatherVCFs.*
    VariantRecalibratorSNP.*
    VariantRecalibratorINDEL.*
    ApplyRecalibrationSNP.*
    ApplyRecalibrationINDEL.*
  }
}
# Generate sets of intervals for scatter-gathering over chromosomes

task CreateSequenceGroupingTSV {
  File Groups

    command {
      cat ${Groups} > /dev/stdout
    }

  output {
    Array[Array[String]] sequence_grouping = read_tsv(stdout())
  }
}

task FastqToSam {
  File PICARD
  File Input_Fastq1
  String Unmapped_Basename

    command {
      time java -Xmx3G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      /Jar/picard.jar \
      FastqToSam \
      FASTQ=${Input_Fastq1} \
      O=${Unmapped_Basename}.bam \
      READ_GROUP_NAME=G \
      SAMPLE_NAME=test \
      LIBRARY_NAME=RH \
      SORT_ORDER=coordinate \
      PLATFORM=ILLUMINA
    }
  
  output {
    File outputbam = "${Unmapped_Basename}.bam"
  }
}

task BwaMem {
  File ref_fasta
  File fasta_bwt
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  File Input_Fastq1 
  String Base_Name
  
    command {
      bwa mem -t 2 -R "@RG\tID:G\tSM:test\tLB:RH\tPL:ILLUMINA\tPU:NotDefined" -M ${ref_fasta} -p ${Input_Fastq1} > ${Base_Name}.sam
    }

  output {
    File outputfile = "${Base_Name}.sam"
  }
}

task MergeBamAlignment {
  File PICARD
  File ref_fasta_index
  File Unmapped_Bam
  File Aligned_Bam
  File ref_fasta
  File ref_dict
  String Output_Bam_Basename

    command {
      java -Xmx3G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      /Jar/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${Aligned_Bam} \
      UNMAPPED_BAM=${Unmapped_Bam} \
      OUTPUT=${Output_Bam_Basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="coordinate" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=200000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="0.7.12-r1039" \
      PROGRAM_GROUP_COMMAND_LINE="bwa mem -t 18 -R -M Input1 Input2 > output.sam" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAP_CONTAMINANT_READS=true
    } 
  output {
    File output_bam = "${Output_Bam_Basename}.bam"
  }
}

task MarkDup {
  Array[File] Input_File
  File PICARD
  String Base_Name
  
    command {
      java -Xmx3G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      /Jar/picard.jar \
      MarkDuplicates \
      I=${sep=' I=' Input_File} \
      O=${Base_Name}.bam \
      VALIDATION_STRINGENCY=LENIENT \
      METRICS_FILE=${Base_Name}.metrics \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=200000 \
      CREATE_INDEX=true
    }
  output {
    File MarkDupOutputBam = "${Base_Name}.bam"
    File MarkDupOutputBai = "${Base_Name}.bai"
    File MetricsFile = "${Base_Name}.metrics"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File GATK3
  File Input_Bam
  File Input_Bam_Index
  String Recalibration_Report_Filename
  Array[File] Sequence_Group_Interval
  File dbsnp_vcf
  File dbsnp_vcf_index
  File v1000g_vcf
  File mills_vcf
  File v1000g_vcf_index
  File mills_vcf_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -XX:ParallelGCThreads=1 -Xmx3G \
      -jar ${GATK3} \
      -T BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${Input_Bam} \
      -o ${Recalibration_Report_Filename} \
      -knownSites ${dbsnp_vcf} \
      -knownSites ${v1000g_vcf} \
      -knownSites ${mills_vcf} \
      -L ${sep=" -L " Sequence_Group_Interval} \
      -XL hs37d5 \
      -cov ContextCovariate \
      -cov CycleCovariate
  }
  output {
    File Recalibration_Report = "${Recalibration_Report_Filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task PrintReads {
  File GATK3
  File Input_Bam
  File Input_Bam_Index
  File Recalibration_Report
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Array[File] Sequence_Group_Interval
  String Output_Bam_Basename

  command {
    java -Xmx3G -XX:ParallelGCThreads=4 \
      -jar ${GATK3} \
      -T PrintReads \
      -R ${ref_fasta} \
      -I ${Input_Bam} \
      -o ${Output_Bam_Basename}.bam \
      -BQSR ${Recalibration_Report} \
      -L ${sep=" -L " Sequence_Group_Interval} \
      -L unmapped \
      -XL hs37d5
  }
  output {
    File recalibrated_bam = "${Output_Bam_Basename}.bam"
    File recalibrated_bam_index = "${Output_Bam_Basename}.bai"

  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
#       INPUT=${Input_Unmapped_Reads_Bam} \
task GatherBamFiles {
  File PICARD
  Array[File] Input_Bams
  String Output_Bam_Basename

  command {
    java -Xmx3G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
      /Jar/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' Input_Bams} \
      OUTPUT=${Output_Bam_Basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
  }
  output {
    File output_bam = "${Output_Bam_Basename}.bam"
    File output_bam_index = "${Output_Bam_Basename}.bai"
    File output_bam_md5 = "${Output_Bam_Basename}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  File GATK3
  Array[File] Input_Bqsr_Reports
  String Output_Report_Filename

  command {
    java -Xmx3G -cp \
      ${GATK3} \
      org.broadinstitute.gatk.tools.GatherBqsrReports \
      I=${sep=' I=' Input_Bqsr_Reports} \
      O=${Output_Report_Filename}
  }
  output {
    File output_bqsr_report = "${Output_Report_Filename}"
  }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File GATK3
  File Input_Bam
  File Input_Bam_Index
  Array[File] Sequence_Group_Interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String Gvcf_Basename

  command {
    java -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3G \
      -jar ${GATK3} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${Gvcf_Basename}.g.vcf \
      -I ${Input_Bam} \
      -L ${sep=" -L " Sequence_Group_Interval} \
      -XL hs37d5 \
      -ERC GVCF
  }

  output {
    File output_gvcf = "${Gvcf_Basename}.g.vcf"
    File output_gvcf_index = "${Gvcf_Basename}.g.vcf.idx"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task GatherVCFs {
  File PICARD
  Array[File] Input_Vcfs
  Array[File] Input_Vcfs_Indexes
  String Output_Vcf_Name

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx3G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
    /Jar/picard.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' Input_Vcfs} \
    OUTPUT=${Output_Vcf_Name}.g.vcf
  }
  output {
    File output_vcfs = "${Output_Vcf_Name}.g.vcf"
    File output_vcfs_index = "${Output_Vcf_Name}.g.vcf.idx"
  }
}

task GenotypeGVCFs {
  File GATK3
  Array[File] Input_Vcf
  Array[File] Input_Vcf_Index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Name

  command {
    java -Xmx8G -XX:ParallelGCThreads=4 -jar \
    ${GATK3} \
    -T GenotypeGVCFs \
    -nt 4 \
    -R ${ref_fasta} \
    -o ${Output_Name}.g.vcf \
    --variant ${sep=' --variant ' Input_Vcf}
  }
  output {
    File output_vcf = "${Output_Name}.g.vcf"
    File output_vcf_index = "${Output_Name}.g.vcf.idx"
  }
}

task VariantRecalibratorSNP {
  File GATK3
  File Input_Vcf
  File Input_Vcf_Index
  File v1000g_vcf
  File omni_vcf
  File dbsnp_vcf
  File hapmap_vcf
  File v1000g_vcf_index
  File omni_vcf_index
  File dbsnp_vcf_index
  File hapmap_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Vcf_Name
  String Mode

  command {
    java -Xmx3G -XX:ParallelGCThreads=4 -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -nt 6 \
    -R ${ref_fasta} \
    -input ${Input_Vcf} \
    -mode ${Mode} \
    -resource:v1000G,known=false,training=true,truth=false,prior=10.0 ${v1000g_vcf} \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni_vcf} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf} \
    -an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
    -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 \
    -tranche 97.0 -tranche 90.0 \
    -recalFile ${Output_Vcf_Name}.recal \
    -tranchesFile ${Output_Vcf_Name}.tranches \
    -rscriptFile ${Output_Vcf_Name}.plots.R
  }
  output {
    File recalFile = "${Output_Vcf_Name}.recal"
    File tranchesFile = "${Output_Vcf_Name}.tranches"
    File rscriptFile = "${Output_Vcf_Name}.plots.R"
  }
}

task VariantRecalibratorINDEL {
  File GATK3
  File Input_Vcf
  File Input_Vcf_Index
  File dbsnp_vcf
  File mills_vcf
  File mills_vcf_index
  File dbsnp_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Vcf_Name
  String Mode

  command {
    java -Xmx3G -XX:ParallelGCThreads=4 -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -nt 6 \
    -R ${ref_fasta} \
    -input ${Input_Vcf} \
    -mode ${Mode} \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_vcf} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf} \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
    -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 \
    -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -recalFile ${Output_Vcf_Name}.recal \
    -tranchesFile ${Output_Vcf_Name}.tranches \
    -rscriptFile ${Output_Vcf_Name}.plots.R
  }
  output {
    File recalFile = "${Output_Vcf_Name}.recal"
    File tranchesFile = "${Output_Vcf_Name}.tranches"
    File rscriptFile = "${Output_Vcf_Name}.plots.R"
  }
}

task ApplyRecalibrationSNP {
  File GATK3
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File Input_Vcf
  File TranchesFile
  File RecalFile
  String Mode
  String Output_Vcf_Name

  command {
    java -jar -Xmx3G \
    ${GATK3} \
    -T ApplyRecalibration \
    -nt 6 \
    -input ${Input_Vcf} \
    -R ${ref_fasta} \
    -mode ${Mode} \
    --ts_filter_level 99.6 \
    -tranchesFile ${TranchesFile} \
    -recalFile ${RecalFile} \
    -o ${Output_Vcf_Name}.g.vcf
  }
  output {
    File output_vcf = "${Output_Vcf_Name}.g.vcf"
    File output_vcf_index = "${Output_Vcf_Name}.g.vcf.idx"
  }
}

task ApplyRecalibrationINDEL {
  File GATK3
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File Input_Vcf
  File TranchesFile
  File RecalFile
  String Mode
  String Output_Vcf_Name

  command {
    java -jar -Xmx3G \
    ${GATK3} \
    -T ApplyRecalibration \
    -nt 6 \
    -input ${Input_Vcf} \
    -R ${ref_fasta} \
    -mode ${Mode} \
    --ts_filter_level 95.0 \
    -tranchesFile ${TranchesFile} \
    -recalFile ${RecalFile} \
    -o ${Output_Vcf_Name}.g.vcf
  }
  output {
    File output_vcf = "${Output_Vcf_Name}.g.vcf"
    File output_vcf_index = "${Output_Vcf_Name}.g.vcf.idx"
  }
}
