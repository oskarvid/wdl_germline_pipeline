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
#     Input_Fastq2 = input_fastq2,
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
  #    Input_Fastq2 = input_fastq2,
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

  # Outputs that will be retained when execution is complete
  
  output {
    MarkDup.*
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
#  File Input_Fastq2
  String Unmapped_Basename

    command {
      time java -Xmx3G -XX:ParallelGCThreads=16 -Djava.io.tmpdir=`pwd`/tmp -jar \
      ${PICARD} \
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
#  File Input_Fastq2
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
      ${PICARD} \
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
      ${PICARD} \
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
