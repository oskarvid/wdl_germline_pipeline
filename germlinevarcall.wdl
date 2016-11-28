workflow GermlineVarCall {
# Input files
  Array[File] scattered_calling_intervals
  Array[File] known_indels_sites_indices
  Array[File] known_indels_sites_VCFs
  File dbSNP_vcf_index
  File dbSNP_vcf
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
  File vrresource1
  File vrresource2
  File vrresource3
  File vrresource4
  File vrresource5
  File vrresource1_index
  File vrresource2_index
  File vrresource3_index
  File vrresource4_index
  File vrresource5_index

# Jar files  
  File picard
  File gatk3
  File gatk4

# String names  
  String Base_Name
  String final_gvcf_name
  String recalibrated_bam_basename = Base_Name + ".aligned.duplicates_marked.recalibrated"

call CreateSequenceGroupingTSV {
  input:
    ref_dict = ref_dict,
}

call BwaMem {
  input:
    Ref_Fasta = ref_fasta,
    Fasta_Bwt = fasta_bwt,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    ref_bwt = ref_bwt,
    ref_amb = ref_amb,
    ref_ann = ref_ann,
    ref_pac = ref_pac,
    ref_sa = ref_sa,
    Input_Fastq1 = input_fastq1,
    Input_Fastq2 = input_fastq2,
    Base_Name = Base_Name
}

call SortSam {
  input:
    PICARD = picard,
    Base_Name = Base_Name +".sortsam.bwa",
    Input_File = BwaMem.outputfile,
    Ref_Fasta = ref_fasta
}

call SetNm {
  input:
    Ref_Fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    PICARD = picard,
    Base_Name = Base_Name +".setnm.sortsam.bwa",
    Input_Bam = SortSam.SamOutputBam,
}

call MarkDup {
  input:
    PICARD = picard,
    Base_Name = Base_Name + ".markdup.setnm.sortsam.bwa",
    Input_File = SortSam.SamOutputBam
}
    
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        GATK4=gatk4,
        Input_Bam = MarkDup.MarkDupOutputBam,
        Input_Bam_Index = MarkDup.MarkDupOutputBai,
        Recalibration_Report_Filename = Base_Name + ".recal_data.csv",
        Sequence_Group_Interval = subgroup,
        DbSNP_Vcf = dbSNP_vcf,
        DbSNP_Vcf_Index = dbSNP_vcf_index,
        Known_Indels_Sites_VCFs = known_indels_sites_VCFs,
        Known_Indels_Sites_Indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }  
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        GATK4=gatk4,
        Input_Bam = MarkDup.MarkDupOutputBam,
        Input_Bam_Index = MarkDup.MarkDupOutputBai,
        Output_Bam_Basename = recalibrated_bam_basename,
        Recalibration_Report = GatherBqsrReports.output_bqsr_report,
        Sequence_Group_Interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
  }
  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports{
    input:
      Input_Bqsr_Reports = BaseRecalibrator.Recalibration_Report,
      Output_Report_Filename = Base_Name + ".recal_data.csv",
      GATK4 = gatk4,
  }
    
  # Do an additional round of recalibration on the unmapped reads (which would otherwise 
  # be left behind because they're not accounted for in the scatter intervals). This is 
  # done by running ApplyBQSR with "-L unmapped".
  Array[String] unmapped_group_interval = ["unmapped"]
  call ApplyBQSR as ApplyBQSRToUnmappedReads {
    input:
      GATK4=gatk4,
      Input_Bam = MarkDup.MarkDupOutputBam,
      Input_Bam_Index = MarkDup.MarkDupOutputBai,
      Output_Bam_Basename = recalibrated_bam_basename,
      Recalibration_Report = GatherBqsrReports.output_bqsr_report,
      Sequence_Group_Interval = unmapped_group_interval,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
  }
    
  # Merge the recalibrated BAM files resulting from by-interval recalibration
  # TODO: when we have capability of adding elements to arrays, can just have one array 
  # as an input and add the output of the above task to the scattered printreads bams
  call GatherBamFiles {
    input:
      PICARD=picard,
      Input_Bams = ApplyBQSR.recalibrated_bam,
      Input_Unmapped_Reads_Bam = ApplyBQSRToUnmappedReads.recalibrated_bam,
      Output_Bam_Basename = Base_Name + ".bqsr.baserecal.markdup.setnm.sortsam.bwa",
  }

  # Call variants in parallel over WGS calling intervals
  scatter (subInterval in scattered_calling_intervals) {
  
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        GATK3 = gatk3,
        Input_Bam = GatherBamFiles.output_bam,
        Input_Bam_Index = GatherBamFiles.output_bam_index,
        Interval_List = subInterval,
        Gvcf_Basename = Base_Name + ".haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
    }
  }
  
  # Combine by-interval GVCFs into a single sample GVCF file
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
      Input_Vcf = GatherVCFs.output_vcfs,
      Input_Vcf_Index = GatherVCFs.output_vcfs_index,
      Output_Name = final_gvcf_name + ".genotypegvcf.haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
  }

  call VariantRecalibratorSNP {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf_Index = GenotypeGVCFs.output_vcf_index,
      GATK3 = gatk3,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      VrResource1 = vrresource1,
      VrResource2 = vrresource2,
      VrResource3 = vrresource3,
      VrResource4 = vrresource4,
      VrResource1_Index = vrresource1_index,
      VrResource2_Index = vrresource2_index,
      VrResource3_Index = vrresource3_index,
      VrResource4_Index = vrresource4_index,
      Mode = "SNP",
      Output_Vcf_Name = final_gvcf_name + ".SNP.genotypegvcf.haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
  }

  call VariantRecalibratorINDEL {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      Input_Vcf_Index = GenotypeGVCFs.output_vcf_index,
      Input_Vcf = GenotypeGVCFs.output_vcf,
      GATK3 = gatk3,
      VrResource5 = vrresource5,
      VrResource5_Index = vrresource5_index,
      Mode = "INDEL",
      Output_Vcf_Name = final_gvcf_name + ".INDEL.genotypegvcf.haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
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
      Output_Vcf_Name = final_gvcf_name + ".applyrecal.SNP.genotypegvcf.haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
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
      Output_Vcf_Name = final_gvcf_name + ".applyrecal.INDEL.genotypegvcf.haplotypecaller.bqsr.baserecal.markdup.setnm.sortsam.bwa",
  }

  # Outputs that will be retained when execution is complete
  
  output {
    MarkDuplicates.duplicate_metrics
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
  File ref_dict

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.  It outputs to stdout
  # where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    print tsv_string
    CODE
  >>>
  output {
    Array[Array[String]] sequence_grouping = read_tsv(stdout())
  }
}

task BwaMem {
  File Ref_Fasta
  File Fasta_Bwt
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  File Input_Fastq1 
  File Input_Fastq2
  String Base_Name
  
    command {
      bwa mem -t 18 -R "@RG\tID:G\tSM:test\tLB:RH\tPL:ILLUMINA\tPU:NotDefined" -M ${Ref_Fasta} ${Input_Fastq1} ${Input_Fastq2} > ${Base_Name}.sam
    }
  output {
    File outputfile = "${Base_Name}.sam"
  }
}

task SortSam {
  File Input_File
  File PICARD
  File Ref_Fasta
  String Base_Name

    command {
      java -Xmx8G -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${PICARD} \
      SortSam \
      INPUT=${Input_File} \
      OUTPUT=${Base_Name}.bam \
      SORT_ORDER="coordinate" \
      MAX_RECORDS_IN_RAM=1000000 \
      CREATE_MD5_FILE=false \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=LENIENT
    }
  output {
    File SamOutputBam = "${Base_Name}.bam"
    File SamOutputBamIndex = "${Base_Name}.bai"
  }
}

task SetNm {
  File PICARD
  File Input_Bam
  File Ref_Fasta
  File ref_dict
  File ref_fasta_index
  String Base_Name
  
    command {
      java -Xmx8G -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${PICARD} \
      SetNmAndUqTags \
      INPUT=${Input_Bam} \
      OUTPUT=${Base_Name}.setnm.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      REFERENCE_SEQUENCE=${Ref_Fasta}
    }

  output {
    File SamOutputBam = "${Base_Name}.bam"
    File output_bam_index = "${Base_Name}.bai"
    File output_bam_md5 = "${Base_Name}.bam.md5"
  }
}

task MarkDup {
  File Input_File
  File PICARD
  String Base_Name
  
    command {
      java -Xmx8G -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${PICARD} \
      MarkDuplicates \
      I=${Input_File} \
      O=${Base_Name}.bam \
      VALIDATION_STRINGENCY=LENIENT \
      METRICS_FILE=samples.metrics \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
      CREATE_INDEX=true
    }
  output {
    File MarkDupOutputBam = "${Base_Name}.bam"
    File MarkDupOutputBai = "${Base_Name}.bai"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File GATK4
  File Input_Bam
  File Input_Bam_Index
  String Recalibration_Report_Filename
  Array[String] Sequence_Group_Interval
  File DbSNP_Vcf
  File DbSNP_Vcf_Index
  Array[File] Known_Indels_Sites_VCFs
  Array[File] Known_Indels_Sites_Indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx8G \
      -jar ${GATK4} \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${Input_Bam} \
      --useOriginalQualities \
      -O ${Recalibration_Report_Filename} \
      -knownSites ${DbSNP_Vcf} \
      -knownSites ${sep=" -knownSites " Known_Indels_Sites_VCFs} \
      -L ${sep=" -L " Sequence_Group_Interval}
  }
  output {
    File Recalibration_Report = "${Recalibration_Report_Filename}"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File GATK4
  File Input_Bam
  File Input_Bam_Index
  String Output_Bam_Basename
  File Recalibration_Report
  Array[String] Sequence_Group_Interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  command {
    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8G \
      -jar ${GATK4} \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${Input_Bam} \
      --useOriginalQualities \
      -O ${Output_Bam_Basename}.bam \
      -bqsr ${Recalibration_Report} \
      -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40 \
      --emit_original_quals \
      -L ${sep=" -L " Sequence_Group_Interval}
  }
  output {
    File recalibrated_bam = "${Output_Bam_Basename}.bam"
    File recalibrated_bam_checksum = "${Output_Bam_Basename}.bam.md5"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  File GATK4
  Array[File] Input_Bqsr_Reports
  String Output_Report_Filename

  command {
    java -Xmx8G -jar \
      ${GATK4} \
      GatherBQSRReports \
      -I ${sep=' -I ' Input_Bqsr_Reports} \
      -O ${Output_Report_Filename}
  }
  output {
    File output_bqsr_report = "${Output_Report_Filename}"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  File PICARD
  Array[File] Input_Bams
  File Input_Unmapped_Reads_Bam
  String Output_Bam_Basename

  command {
    java -Xmx8G -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${PICARD} \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' Input_Bams} \
      INPUT=${Input_Unmapped_Reads_Bam} \
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

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File GATK3
  File Input_Bam
  File Input_Bam_Index
  File Interval_List
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String Gvcf_Basename

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8G \
      -jar ${GATK3} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${Gvcf_Basename}.g.vcf \
      -I ${Input_Bam} \
      -L ${Interval_List} \
      -ERC GVCF
  }
#      -variant_index_parameter 128000 \
#      -variant_index_type LINEAR

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
    java -Xmx8G -Djava.io.tmpdir=`pwd`/tmp -jar \
    ${PICARD} \
    MergeVcfs \
    INPUT=${sep=' INPUT=' Input_Vcfs} \
    OUTPUT=${Output_Vcf_Name}.vcf
  }
  output {
    File output_vcfs = "${Output_Vcf_Name}.vcf"
    File output_vcfs_index = "${Output_Vcf_Name}.vcf.idx"
  }
}

task GenotypeGVCFs {
  File GATK3
  File Input_Vcf
  File Input_Vcf_Index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Name

  command {
    java -Xmx8G -jar \
    ${GATK3} \
    -T GenotypeGVCFs \
    -nt 18 \
    -R ${ref_fasta} \
    -o ${Output_Name}.vcf \
    --variant ${Input_Vcf}
  }
  output {
    File output_vcf = "${Output_Name}.vcf"
    File output_vcf_index = "${Output_Name}.vcf.idx"
  }
}

task VariantRecalibratorSNP {
  File GATK3
  File Input_Vcf
  File Input_Vcf_Index
  File VrResource1
  File VrResource2
  File VrResource3
  File VrResource4
  File VrResource1_Index
  File VrResource2_Index
  File VrResource3_Index
  File VrResource4_Index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Vcf_Name
  String Mode

  command {
    java -Xmx8G -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -R ${ref_fasta} \
    -input ${Input_Vcf} \
    -mode ${Mode} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${VrResource1} \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ${VrResource2} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${VrResource3} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${VrResource4} \
    -an QD -an MQ -an ReadPosRankSum -an FS -an SOR -tranche 100.0 -tranche 99.95 \
    -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
    -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
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
  File VrResource5
  File VrResource5_Index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String Output_Vcf_Name
  String Mode

  command {
    java -Xmx8G -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -R ${ref_fasta} \
    -input ${Input_Vcf} \
    -mode ${Mode} \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 ${VrResource5} \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -tranche 100.0 -tranche 99.95 \
    -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 \
    -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -recalFile ${Output_Vcf_Name}.recal \
    -tranchesFile ${Output_Vcf_Name}.tranches \
    -rscriptFile ${Output_Vcf_Name}.plots.R \
    -mG 4
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
    java -jar -Xmx8G \
    ${GATK3} \
    -T ApplyRecalibration \
    -input ${Input_Vcf} \
    -R ${ref_fasta} \
    -mode ${Mode} \
    --ts_filter_level 99.6 \
    -tranchesFile ${TranchesFile} \
    -recalFile ${RecalFile} \
    -o ${Output_Vcf_Name}.vcf
  }
  output {
    File output_vcf = "${Output_Vcf_Name}.vcf"
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
    java -jar -Xmx8G \
    ${GATK3} \
    -T ApplyRecalibration \
    -input ${Input_Vcf} \
    -R ${ref_fasta} \
    -mode ${Mode} \
    --ts_filter_level 95.0 \
    -tranchesFile ${TranchesFile} \
    -recalFile ${RecalFile} \
    -o ${Output_Vcf_Name}.vcf
  }
  output {
    File output_vcf = "${Output_Vcf_Name}.vcf"
  }
}
