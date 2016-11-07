workflow GermlineVarCall {
# Input files
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

# Jar files  
  File picard
  File gatk3
  File gatk4

# String names  
  String Base_Name


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
    Base_Name = Base_Name +".bwa",
    Input_File = BwaMem.outputfile,
    Ref_Fasta = ref_fasta
    }

call SetNm {
  input:
    Ref_Fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    PICARD = picard,
    Base_Name = Base_Name +".sortsam.bwa",
    Input_Bam = SortSam.SamOutputBam,
}

call MarkDup {
  input:
    PICARD = picard,
    Base_Name = Base_Name + ".bwa.sortsam",
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
      Input_Bam = sortSamBam,
      Input_Bam_Index = sortSamBamIndex,
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
      GATK4=gatk4,
      Input_Bams = ApplyBQSR.recalibrated_bam,
      Input_Unmapped_Reads_Bam = ApplyBQSRToUnmappedReads.recalibrated_bam,
      Output_Bam_Basename = Base_Name,
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
      bwa mem -t 2 -K 100000000 -M -v 3 ${Ref_Fasta} ${Input_Fastq1} ${Input_Fastq2} > ${Base_Name}.sam
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
      java -Xmx8G -jar \
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
      java -Xmx500m -jar \
      ${PICARD} \
      SetNmAndUqTags \
      INPUT=${Input_Bam} \
      OUTPUT=${Base_Name}.setnm.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      REFERENCE_SEQUENCE=${Ref_Fasta}
      }

  output {
    File SamOutputBam = "${Base_Name}.setnm.bam"
    File output_bam_index = "${Base_Name}.setnm.bai"
    File output_bam_md5 = "${Base_Name}.setnm.bam.md5"
  }
}


task MarkDup {
  File Input_File
  File PICARD
  String Base_Name
  
    command {
      java -jar -Xmx8G \
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
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx4000m \
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
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
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
    java -Xmx3000m -jar \
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
  File GATK4
  Array[File] Input_Bams
  File Input_Unmapped_Reads_Bam
  String Output_Bam_Basename

  command {
    java -Xmx2000m -jar \
      ${GATK4} \
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








