#!/usr/bin/env python

import os           # Import the os module for basic path manipulation
import arvados      # Import the Arvados sdk module
import re
import subprocess

# the amount to weight each sequence contig
weight_seq = 120000

class InvalidArgumentError(Exception):
    pass

class FileAccessError(Exception):
    pass

class APIError(Exception):
    pass

def prepare_gatk_reference_collection(reference_coll):
    """
    Checks that the supplied reference_collection has the required 
    files and only the required files for GATK. 
    Returns: a portable data hash for the reference collection
    """
    # Ensure we have a .fa reference file with corresponding .fai index and .dict
    # see: http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
    rcr = arvados.CollectionReader(reference_coll)
    ref_fasta = {}
    ref_fai = {}
    ref_dict = {}
    ref_input = None
    dict_reader = None
    for rs in rcr.all_streams():
        for rf in rs.all_files():
            if re.search(r'\.fa$', rf.name()):
                ref_fasta[rs.name(), rf.name()] = rf
            elif re.search(r'\.fai$', rf.name()):
                ref_fai[rs.name(), rf.name()] = rf
            elif re.search(r'\.dict$', rf.name()):
                ref_dict[rs.name(), rf.name()] = rf
    for ((s_name, f_name), fasta_f) in ref_fasta.items():
        fai_f = ref_fai.get((s_name, re.sub(r'fa$', 'fai', f_name)), 
                            ref_fai.get((s_name, re.sub(r'fa$', 'fa.fai', f_name)), 
                                        None))
        dict_f = ref_dict.get((s_name, re.sub(r'fa$', 'dict', f_name)), 
                              ref_dict.get((s_name, re.sub(r'fa$', 'fa.dict', f_name)), 
                                           None))
        if fasta_f and fai_f and dict_f:
            # found a set of all three! 
            ref_input = fasta_f.as_manifest()
            ref_input += fai_f.as_manifest()
            ref_input += dict_f.as_manifest()
            break
    if ref_input is None:
        raise InvalidArgumentError("Expected a reference fasta with fai and dict in reference_collection. Found [%s]" % ' '.join(rf.name() for rf in rs.all_files()))
    # Create and return a portable data hash for the ref_input manifest
    try:
        r = arvados.api().collections().create(body={"manifest_text": ref_input}).execute()
        ref_input_pdh = r["portable_data_hash"]
    except:
        raise 
    return ref_input_pdh

def create_interval_lists(genome_chunks, reference_coll, skip_sq_sn_r):
    rcr = arvados.CollectionReader(reference_coll)
    ref_dict = []
    dict_reader = None
    for rs in rcr.all_streams():
        for rf in rs.all_files():
            if re.search(r'\.dict$', rf.name()):
                ref_dict.append(rf)
    if len(ref_dict) < 1:
        raise InvalidArgumentError("Reference collection does not contain any .dict files but one is required.")
    if len(ref_dict) > 1:
        raise InvalidArgumentError("Reference collection contains multiple .dict files but only one is allowed.")
    dict_reader = ref_dict[0]

    # Load the dict data
    interval_header = ""
    dict_lines = dict_reader.readlines()
    dict_header = dict_lines.pop(0)
    if re.search(r'^@HD', dict_header) is None:
        raise InvalidArgumentError("Dict file in reference collection does not have correct header: [%s]" % dict_header)
    interval_header += dict_header
    print "Dict header is %s" % dict_header
    sn_intervals = dict()
    sns = []
    total_len = 0
    for sq in dict_lines:
        if re.search(r'^@SQ', sq) is None:
            raise InvalidArgumentError("Dict file contains malformed SQ line: [%s]" % sq)
        interval_header += sq
        sn = None
        ln = None
        for tagval in sq.split("\t"):
            tv = tagval.split(":", 1)
            if tv[0] == "SN":
                sn = tv[1]
            if tv[0] == "LN":
                ln = tv[1]
            if sn and ln:
                break
        if not (sn and ln):
            raise InvalidArgumentError("Dict file SQ entry missing required SN and/or LN parameters: [%s]" % sq)
        assert(sn and ln)
        if sn_intervals.has_key(sn):
            raise InvalidArgumentError("Dict file has duplicate SQ entry for SN %s: [%s]" % (sn, sq))
        if skip_sq_sn_r.search(sn):
            next
        sn_intervals[sn] = (1, int(ln))
        sns.append(sn)
        total_len += int(ln)
    total_sequences = len(sns)

    # Chunk the genome into genome_chunks equally sized pieces and create intervals files
    print "Total sequences included: %s" % (total_sequences)
    print "Total genome length is %s" % total_len
    total_points = total_len + (total_sequences * weight_seq)
    print "Total points to split: %s" % (total_points)
    chunk_points = int(total_points / genome_chunks)
    chunks_c = arvados.collection.CollectionWriter(num_retries=3)
    print "Chunking genome into %s chunks of ~%s points" % (genome_chunks, chunk_points)
    for chunk_i in range(0, genome_chunks):
        chunk_num = chunk_i + 1
        chunk_intervals_count = 0
        chunk_input_name = dict_reader.name() + (".%s_of_%s.interval_list" % (chunk_num, genome_chunks))
        print "Creating interval file for chunk %s" % chunk_num
        chunks_c.start_new_file(newfilename=chunk_input_name)
        chunks_c.write(interval_header)
        remaining_points = chunk_points
        while len(sns) > 0:
            sn = sns.pop(0)
            remaining_points -= weight_seq
            if remaining_points <= 0:
                sns.insert(0, sn)
                break
            if not sn_intervals.has_key(sn):
                raise ValueError("sn_intervals missing entry for sn [%s]" % sn)
            start, end = sn_intervals[sn]
            if (end-start+1) > remaining_points:
                # not enough space for the whole sq, split it
                real_end = end
                end = remaining_points + start - 1
                assert((end-start+1) <= remaining_points)
                sn_intervals[sn] = (end+1, real_end)
                sns.insert(0, sn)
            interval = "%s\t%s\t%s\t+\t%s\n" % (sn, start, end, "interval_%s_of_%s_%s" % (chunk_num, genome_chunks, sn))
            remaining_points -= (end-start+1)
            chunks_c.write(interval)
            chunk_intervals_count += 1
            if remaining_points <= 0:
                break
        if chunk_intervals_count > 0:
            print "Chunk intervals file %s saved." % (chunk_input_name)
        else:
            print "WARNING: skipping empty intervals for %s" % chunk_input_name
    chunk_input_pdh = chunks_c.finish()
    print "Chunk intervals collection saved as: %s" % (chunk_input_pdh)
    return chunk_input_pdh

def main():
    current_job = arvados.current_job()
    skip_sq_sn_regex = '_decoy$'
    if 'skip_sq_sn_regex' in current_job['script_parameters']:
        skip_sq_sn_regex = current_job['script_parameters']['skip_sq_sn_regex']
    skip_sq_sn_r = re.compile(skip_sq_sn_regex)

    genome_chunks = int(current_job['script_parameters']['genome_chunks'])
    if genome_chunks < 1:
        raise InvalidArgumentError("genome_chunks must be a positive integer")

    # Limit the scope of the reference collection to only those files relevant to gatk
    ref_input_pdh = prepare_gatk_reference_collection(reference_coll=current_job['script_parameters']['reference_collection'])

    # Create an interval_list file for each chunk based on the .dict in the reference collection
    output_locator = create_interval_lists(genome_chunks, ref_input_pdh, skip_sq_sn_r)

    # Use the resulting locator as the output for this task.
    arvados.current_task().set_output(output_locator)

    # Done!


if __name__ == '__main__':
    main()

