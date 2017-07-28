workflow GermlineVarCall {
# Input files
  Array[File] known_indels_sites_indices
  Array[File] known_indels_sites_VCFs
  Array[File] scatter_interval_list
  File inputSamplesFile
#  Array[File] fastqfilesR1
#  Array[File] fastqfilesR2

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

# Jar files  
  File gatk4
  File gatk3

# String names  
  String Base_Name
  String final_gvcf_name
  String recalibrated_bam_basename = Base_Name + ".aligned.duplicates_marked.recalibrated"
  String outputfolder = "/wdl_pipeline/"
  String unmapped_basename = "unmapped_bam"

  # Preprocess sample input file to remove comments and header for proper scatter gather execution
  call FixInputSamplesFile {
    input:
      SamplesFile = inputSamplesFile,
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
  }

scatter (element in FixInputSamplesFile.FixedSamplesFile) {
  call FastqToSam {
   input:
      GATK4 = gatk4,
      ID = element[1] + "-" + element[2],
      LB = element[5],
      SM = element[0],
      PL = element[6],
      Input_Fastq1 = element[3],
      Input_Fastq2 = element[4],
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
      ID = element[1] + "-" + element[2],
      LB = element[5],
      SM = element[0],
      PL = element[6],
      Input_Fastq1 = element[3],
      Input_Fastq2 = element[4],
      Base_Name = Base_Name + ".bwa",
  }

  call MergeBamAlignment {
    input:
      GATK4 = gatk4,
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
      GATK4 = gatk4,
      Base_Name = Base_Name + ".markdup.sortsam.bwa",
      Input_File = MergeBamAlignment.output_bam,
  }
    
  # Perform Base Quality Score Recalibration (BQSR) and call variants on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        GATK4 = gatk4,
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
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      Input_Bqsr_Reports = BaseRecalibrator.Recalibration_Report,
      Output_Report_Filename = Base_Name + ".recal_data.grp",
      GATK4 = gatk4,
  }

scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
  # Apply the recalibration model by interval
  call ApplyBQSR {
    input:
      GATK4 = gatk4,
      Input_Bam = MarkDup.MarkDupOutputBam,
      Input_Bam_Index = MarkDup.MarkDupOutputBai,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      Recalibration_Report = GatherBqsrReports.output_bqsr_report,
      Output_Bam_Basename = recalibrated_bam_basename,
      Sequence_Group_Interval = subgroup,
  }
}

  # Gather bam files from ApplyBQSR and use input for HaplotypeCaller
  call GatherBamFiles {
  input:
    GATK4 = gatk4,
    Input_Bams = ApplyBQSR.recalibrated_bam,
    Output_Bam_Basename = Base_Name + ".bqsr.baserecal.markdup.sortsam.bwa",
  }
    
  # Generate GVCFs
  scatter (subInterval in scatter_interval_list) {
    call HaplotypeCaller {
      input:
        GATK3 = gatk3,
        Input_Bam = GatherBamFiles.output_bam,
        Input_Bam_Index = GatherBamFiles.output_bam_index,
        Sequence_Group_Interval = subInterval,
        Gvcf_Basename = Base_Name + ".haplotypecaller.bqsr.baserecal.markdup.sortsam.bwa",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
  }
}

  # Combine GVCFs into a single sample GVCF file
  call GatherVCFs {
    input:
      GATK4 = gatk4,
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

# Preprocess sample input file to remove comments and header for proper scatter gather execution
task FixInputSamplesFile {
  File SamplesFile

  command {
    grep -v "#" ${SamplesFile} | grep -v "SAMPLE" > /dev/stdout
  }
  output {
    Array[Array[String]] FixedSamplesFile = read_tsv(stdout())
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
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
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
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
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

task FastqToSam {
  File GATK4
  File Input_Fastq1
  File Input_Fastq2
  String ID
  String SM
  String LB
  String PL
  String Unmapped_Basename

    command {
      time java -Xmx16G -Dsnappy.disable=true -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${GATK4} \
      FastqToSam \
      --FASTQ ${Input_Fastq1} \
      --FASTQ2 ${Input_Fastq2} \
      -O ${Unmapped_Basename}.bam \
      --SAMPLE_NAME ${SM} \
      --READ_GROUP_NAME ${ID} \
      --LIBRARY_NAME ${LB} \
      --PLATFORM ${PL} \
      --SORT_ORDER coordinate

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
  File Input_Fastq2
  String ID
  String SM
  String LB
  String PL
  String Base_Name
  
    command {
      bwa mem -t 2 -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:${PL}\tPU:NotDefined" -M ${ref_fasta} ${Input_Fastq1} ${Input_Fastq2} > ${Base_Name}.sam
    }

  output {
    File outputfile = "${Base_Name}.sam"
  }
}

task MergeBamAlignment {
  File GATK4
  File ref_fasta_index
  File Unmapped_Bam
  File Aligned_Bam
  File ref_fasta
  File ref_dict
  String Output_Bam_Basename

    command {
      java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${GATK4} \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ${Aligned_Bam} \
      --UNMAPPED_BAM ${Unmapped_Bam} \
      -O ${Output_Bam_Basename}.bam \
      --reference ${ref_fasta} \
      --SORT_ORDER coordinate \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 200000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --PROGRAM_RECORD_ID "bwamem" \
      --PROGRAM_GROUP_VERSION "0.7.12-r1039" \
      --PROGRAM_GROUP_COMMAND_LINE "bwa mem -t 18 -R -M Input1 Input2 > output.sam" \
      --PROGRAM_GROUP_NAME "bwamem"
    } 
  output {
    File output_bam = "${Output_Bam_Basename}.bam"
  }
}

task MarkDup {
  Array[File] Input_File
  File GATK4
  String Base_Name
  
    command {
      java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${GATK4} \
      MarkDuplicates \
      --input ${sep=' --input=' Input_File} \
      -O ${Base_Name}.bam \
      --VALIDATION_STRINGENCY LENIENT \
      --METRICS_FILE ${Base_Name}.metrics \
      --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 200000 \
      --CREATE_INDEX true
    }
  output {
    File MarkDupOutputBam = "${Base_Name}.bam"
    File MarkDupOutputBai = "${Base_Name}.bai"
    File MetricsFile = "${Base_Name}.metrics"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File GATK4
  Array[File] Input_Bam
  Array[File] Input_Bam_Index
  String Recalibration_Report_Filename
  Array[String] Sequence_Group_Interval
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
    java -Dsnappy.disable=true -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx13G \
      -jar ${GATK4} \
      BaseRecalibrator \
      --reference ${ref_fasta} \
      --input ${sep=" --input" Input_Bam} \
      -O ${Recalibration_Report_Filename} \
      --knownSites ${dbsnp_vcf} \
      --knownSites ${v1000g_vcf} \
      --knownSites ${mills_vcf} \
      --intervals ${sep=" --intervals " Sequence_Group_Interval}
#      -cov ContextCovariate \
#      -cov CycleCovariate
  }
  output {
    File Recalibration_Report = "${Recalibration_Report_Filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File GATK4
  Array[File] Input_Bam
  Array[File] Input_Bam_Index
  File Recalibration_Report
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Array[String] Sequence_Group_Interval
  String Output_Bam_Basename

  command {
    java -Dsnappy.disable=true -Xmx13G -XX:ParallelGCThreads=4 \
      -jar ${GATK4} \
      ApplyBQSR \
      --reference ${ref_fasta} \
      --input ${sep=" --input " Input_Bam} \
      -O ${Output_Bam_Basename}.bam \
      --createOutputBamIndex true \
      -bqsr ${Recalibration_Report} \
      --intervals ${sep=" --intervals " Sequence_Group_Interval}
  }
  output {
    File recalibrated_bam = "${Output_Bam_Basename}.bam"
    File recalibrated_bam_index = "${Output_Bam_Basename}.bai"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  File GATK4
  Array[File] Input_Bams
  String Output_Bam_Basename

  command {
    java -Dsnappy.disable=true -Xmx16G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${GATK4} \
      GatherBamFiles \
      --input ${sep=" --input " Input_Bams} \
      -O ${Output_Bam_Basename}.bam \
      --CREATE_INDEX true
  }
  output {
    File output_bam = "${Output_Bam_Basename}.bam"
    File output_bam_index = "${Output_Bam_Basename}.bai"
  }
  runtime {
    continueOnReturnCode: [0, 1, 2]
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  File GATK4
  Array[File] Input_Bqsr_Reports
  String Output_Report_Filename

  command {
    java -Dsnappy.disable=true -Xmx16G -jar \
      ${GATK4} \
      GatherBQSRReports \
      --input ${sep=" --input " Input_Bqsr_Reports} \
      -O ${Output_Report_Filename} \
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
  Array[String] Sequence_Group_Interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String Gvcf_Basename

  command {
    java -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx13G \
      -jar ${GATK3} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${Gvcf_Basename}.g.vcf \
      -I ${Input_Bam} \
      -L ${sep=" -L " Sequence_Group_Interval} \
      -ERC GVCF
  }

  output {
    File output_gvcf = "${Gvcf_Basename}.g.vcf"
    File output_gvcf_index = "${Gvcf_Basename}.g.vcf.idx"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task GatherVCFs {
  File GATK4
  Array[File] Input_Vcfs
  Array[File] Input_Vcfs_Indexes
  String Output_Vcf_Name

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx16G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=`pwd`/tmp -jar \
    ${GATK4} \
    MergeVcfs \
    --input ${sep=" --input " Input_Vcfs} \
    -O ${Output_Vcf_Name}.g.vcf \
    --CREATE_INDEX true
  }
  output {
    File output_vcfs = "${Output_Vcf_Name}.g.vcf"
    File output_vcfs_index = "${Output_Vcf_Name}.g.vcf.idx"
  }
  runtime {
    continueOnReturnCode: [0, 1, 2, 127]
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
    java -Xmx16G -XX:ParallelGCThreads=18 -jar \
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
    java -Xmx16G -XX:ParallelGCThreads=18 -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -nt 2 \
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
    java -Xmx16G -XX:ParallelGCThreads=18 -jar \
    ${GATK3} \
    -T VariantRecalibrator \
    -nt 2 \
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
    java -XX:ParallelGCThreads=18 -jar -Xmx16G \
    ${GATK3} \
    -T ApplyRecalibration \
    -nt 2 \
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
    java -XX:ParallelGCThreads=18 -jar -Xmx16G \
    ${GATK3} \
    -T ApplyRecalibration \
    -nt 2 \
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
