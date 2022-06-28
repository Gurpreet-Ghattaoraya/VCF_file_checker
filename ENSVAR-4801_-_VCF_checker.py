import argparse
import gzip
import math
import os
import re
import statistics
import sys

from Bio import SeqIO
from LogWriter import LogWriter

import pandas as pd


def arg_passer():
    '''Allows user to provide required files as command line input flags

    Arguments
    ---------
    -v / --vcf
        path to VCF file

    -f / --fasta
        path to FASTA file

    Returns
    -------
    vcf_file : string
        user provided paths to VCF file

    fasta_file : string
        user provided paths to FASTA file
    '''

    arg_parse = argparse.ArgumentParser(description = 'VCF file checker')

    arg_parse.add_argument('-v', '--vcf', metavar = '<PATH_TO_VCF_FILE>', help = 'path to vcf file', required = True)
    arg_parse.add_argument('-f', '--fasta', metavar = '<PATH_TO_FASTA_FILE>', help = 'path to fasta file', required = True)
    arg_parse.add_argument('-m', '--memory-save', help = 'change this scripts behaviour to reduce memory usage at the expense of time', action = 'store_true')

    args = arg_parse.parse_args()

    vcf_file = args.vcf
    fasta_file = args.fasta
    memory_mode = args.memory_save

    return(vcf_file, fasta_file, memory_mode)


def is_file_check(file_path, file_type):
    '''Checks to see if user provided path to file is present as a file

    Parameters
    ----------
    file_path : string / file path
        path to file (preferably absolute path)

    file_type : string
        type of file - e.g. VCF file or FASTA file 
            For human readable output/print

    Returns
    -------
    is_file : boolean
        True/False - if provided path points to a file 
    '''

    is_file = os.path.isfile(file_path)

    if is_file == False:
        logger.critical(f'CRITICAL\tProvided path for {file_type} is not a file')

    return(is_file)


def determin_if_zipped_file(file_path):

    with gzip.open(file_path, 'r') as file_handle:
        try:
            file_handle.read(1)
        except OSError:  # In python3.8+ there will be a dedicated exception raised for "bad gzip files" gzip.BadGzipFile:
            compression_verdict = 'uncompressed'
        else:
            compression_verdict = 'compressed'

    return(compression_verdict)


def is_fasta_file(file_path):
    '''Checks to see if user provided path to --fasta is indeed FASTA file

    Parameters
    ----------
    file_path : string / file path
        path to file (preferably absolute path)

    Returns
    -------
    is_file : boolean
        True/False - if contents of provided FASTA file can be accessed by BioPython's FASTA file handling
    '''

    is_file_compressed = determin_if_zipped_file(file_path)

    if is_file_compressed == 'compressed':
        fasta_handle = gzip.open(file_path, 'rt')
    elif is_file_compressed == 'uncompressed':
        fasta_handle = open(file_path, 'r')

    fasta = SeqIO.parse(fasta_handle, 'fasta')
    is_fasta = any(fasta)  # False when `fasta` is empty, i.e. isn't a FASTA file

    fasta_handle.close()

    if is_fasta == False:
        logger.critical(f'Provided -f/--fasta path is not a FASTA file')
        raise()
    elif is_fasta == True:
        logger.info('FASTA file conforms to FASTA file format')

    return(is_fasta, is_file_compressed)


def has_vcf_header(file_path):
    '''Checks to see if user provided path to --vcf is indeed VCF file

    Parameters
    ----------
    file_path : string / file path
        path to file (preferably absolute path)

    Returns
    -------
    is_file : boolean
        True/False - if contents of provided VCF contains a header line with "##fileformat=VCF" format
    
    format : string
        states if file is "compressed" or "uncompressed" for file handling later in script (e.g. gzip)
    '''

    vcf_format_pattern = re.compile('##fileformat=VCFv\d+\.\d+')

    vcf_file_compressed = determin_if_zipped_file(file_path)

    # handle compressed binary file differently to ordinary file
    if vcf_file_compressed == 'compressed':
        with gzip.open(file_path, 'rb') as vcf_handle:
            header = vcf_handle.readline().decode().strip()

    elif vcf_file_compressed == 'uncompressed':
        with open(file_path, 'r') as vcf_handle:
            header = vcf_handle.readline().strip()
    else:
        logger.critical('Unrecognised file for -v/--vcf file')
        sys.exit()

    # Assumes the fileformat line is the first line of the file
    if re.match(vcf_format_pattern, header):
        verdict = True
        logger.info(f'VCF file contains header "{header}"')
    else:
        verdict = False
        logger.critical(f'Provided -v/--vcf path is not a VCF file')
        raise()

    return(verdict, vcf_file_compressed)


def check_user_input(vcf_file, fasta_file):
    '''Manages the file check functions for the user provided VCF and FASTA files and script exit if failed

    Parameters
    ----------
    vcf_file : string / file path
        path to VCF file (preferably absolute path)
        
    fasta_file : string / file path
        path to FASTA file (preferably absolute path)

    Returns
    -------
    file_compression : string
        states if file is "compressed" or "uncompressed" for file handling later in script (e.g. gzip)
    '''

    logger.info('Checking files')

    is_file_verdict = set()

    is_file_verdict.add(is_file_check(vcf_file, '-v/--vcf'))
    is_file_verdict.add(is_file_check(fasta_file, '-f/--fasta'))

    if False in is_file_verdict:
        raise(FileNotFoundError)
    else:
        logger.info('Both VCF and FASTA paths point to files')

    is_fasta_verdict, fasta_compression = is_fasta_file(fasta_file)
    is_vcf, vcf_compression = has_vcf_header(vcf_file)

    return(vcf_compression, fasta_compression)


def read_vcf(vcf_file, file_compression):
    '''Read in VCF file

    Parameters
    ----------
    vcf_file : string / file path
        path to VCF file (preferably absolute path)
        
    file_compression : string
        states if file is "compressed" or "uncompressed" for file handling in this function (i.e. gzip)

    Returns
    -------
    data_content : list
        each element in list is a line of the VCF file - ignores header lines (those starting with "#" or "##")
    
    header : list
        list containing the table headers for the file data contained in "data_content"
    '''

    print('Reading in FASTA file')
    data_content = []

    if file_compression == 'compressed':
        vcf_handle = gzip.open(vcf_file, 'rb')

        for vcf_read in vcf_handle.readlines():
            line = vcf_read.decode().strip()

            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.split('\t')
                continue

            data_content.append(line)

    elif file_compression == 'uncompressed':
        vcf_handle = open(vcf_file, 'r')

        for vcf_read in vcf_handle.readlines():
            line = vcf_read.strip()

            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.split('\t')
                continue

            data_content.append(line)

    vcf_handle.close()
    print('Data read in for VCF file')

    return(data_content, header)


def get_vcf_header(vcf_file, file_compression):

    logger.info('Obtaining VCF header information')

    if file_compression == 'compressed':
        vcf_handle = gzip.open(vcf_file, 'rb')

    elif file_compression == 'uncompressed':
        vcf_handle = open(vcf_file, 'r')

    for vcf_read in vcf_handle:

        if file_compression == 'compressed':
            line = vcf_read.decode().strip()
        elif file_compression == 'uncompressed':
            line = vcf_read.strip()

        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            header = line.split('\t')
            break
        else:
            continue

    vcf_handle.close()

    logger.debug('VCF header information obtained')

    return(header)


def get_vcf_data_generator(vcf_file, file_compression):

    if file_compression == 'compressed':
        with gzip.open(vcf_file, 'rb') as vcf_handle:
            for vcf_read in vcf_handle:
                line = vcf_read.decode().strip()
                if line.startswith('#'):
                    continue
                else:
                    yield(line)

    elif file_compression == 'uncompressed':
        with gzip.open(vcf_file, 'r') as vcf_handle:
            for vcf_read in vcf_handle:
                line = vcf_read.strip()
                if line.startswith('#'):
                    continue
                else:
                    yield(line)


def read_fasta(fasta_file, fasta_file_compression):
    '''Read in FASTA file

    Parameters
    ----------
    fasta_file : string / file path
        path to FASTA file (preferably absolute path)

    Returns
    -------
    data_content : list
        each element in list is a line of the VCF file - ignores header lines (those starting with "#" or "##")
    
    sequences_dict : dictionary
        keys: sequence name (usually chromosome or patch)

        values: sub-dictionary containing 
            chromosome: chromosome sequence is located on
            seq_start: sequence start position
            seq_end: sequence end position
            seq: BioPython Seq object containing sequence of entry
    '''

    sequences_dict = {}

    loci_pattern = re.compile('.*? dna:.*? .*?:GRCh38:(.*?):(\d*?):(\d*?):.*? .*?')

    if fasta_file_compression == 'compressed':
        fasta_handle = gzip.open(fasta_file, 'rt')
    elif fasta_file_compression == 'uncompressed':
        fasta_handle = open(fasta_file, 'r')

    seq_generator = SeqIO.parse(fasta_handle, 'fasta')

    for fasta_sequence in seq_generator:

        chromosome, seq_start, seq_end = re.findall(loci_pattern, fasta_sequence.description)[0]

        try:
            location = {'chromosome': chromosome,
                        'seq_start': seq_start,
                        'seq_end': seq_end,
                        'seq': fasta_sequence.seq}
        except(IndexError):
            print('IndexError raised')
            print(fasta_sequence.description)
            print(location_find)
            raise(IndexError)

            sequences_dict[fasta_sequence.name] = location

    fasta_handle.close()
    print('Sequences from FASTA file read in')

    return(sequences_dict)


def index_fasta(fasta_file, fasta_file_compression):

    '''if fasta_file_compression == 'compressed':
        fasta_handle = gzip.open(fasta_file, 'rt')
    elif fasta_file_compression == 'uncompressed':
        fasta_handle = open(fasta_file, 'r')'''

    logger.info('Indexing FASTA file')
    record_dict = SeqIO.index(fasta_file, 'fasta')
    logger.info('FASTA file indexed')

    logger.debug(f'Inside = {type(record_dict)}')

    print(record_dict)


def get_sequence_from_fasta(fasta_file_path, fasta_file_compression, chromosome_to_look_for):

    loci_pattern = re.compile('.*? dna:.*? .*?:GRCh38:(.*?):(\d*?):(\d*?):.*? .*?')

    if fasta_file_compression == 'compressed':
        fasta_handle = gzip.open(fasta_file_path, 'rt')
    elif fasta_file_compression == 'uncompressed':
        fasta_handle = open(fasta_file_path, 'r')

    seq_generator = SeqIO.parse(fasta_handle, 'fasta')

    for fasta_sequence in seq_generator:

        chromosome, seq_start, seq_end = re.findall(loci_pattern, fasta_sequence.description)[0]

        if chromosome == chromosome_to_look_for:
            sequence = fasta_sequence.seq
            break

    fasta_handle.close()

    return(sequence)


def rate_of_variations(vcf_data, vcf_header, fasta_file, fasta_file_compression, reduced_memory_mode):
    '''Measure the rate of variation.
    How many variations are present per megabase for each sequence 

    Parameters
    ----------
    vcf_data : list
        list of strings containing VCF file data (each element represents a line from the VCF file)
        
    vcf_header : list
        list of strings containing headings corresponding to "vcf_data"
        
    fasta_sequences : dictionary
        keys: sequence name (usually chromosome or patch)

        values: sub-dictionary containing 
            chromosome: chromosome sequence is located on
            seq_start: sequence start position
            seq_end: sequence end position
            seq: BioPython Seq object containing sequence of entry

    Output
    ------
    count : integer
        number of variations in the VCF file for each chromosome

    rate : float (3 decimal places)
        rate of variations per megabase within the FASTA sequence
    '''

    index_chrom = vcf_header.index('#CHROM')
    index_position = vcf_header.index('POS')
    index_vcf_ref_allele = vcf_header.index('REF')
    index_vcf_alt_allele = vcf_header.index('ALT')

    logger.info('Determining number and rate of variations per megabase (MB) of chromosome/patch sequence')

    rate_dict = {}
    total_count = 0
    total_len = 0
    bins_len = 1000000  # the megabase in "variations per megabase"

    fasta_reference = SeqIO.index(fasta_file, 'fasta')

    for line in vcf_data:

        line_split = line.split('\t')
        chromosome = line_split[index_chrom].replace('chr', '')
        position = int(line_split[index_position])

        if not chromosome in rate_dict:
            rate_dict[chromosome] = {'variation_count': 0}

        rate_dict[chromosome]['variation_count'] += 1

    for chromosome in rate_dict:

        '''if reduced_memory_mode == True:
            chromosome_length = len(get_sequence_from_fasta(fasta_reference, fasta_file_compression, chromosome))
        elif reduced_memory_mode == False:
            chromosome_length = len(fasta_reference[chromosome]['seq'])'''
        try:
            chromosome_length = len(fasta_reference[chromosome].seq)  # Delete if keeping above triple comment block
        except(KeyError):
            chromosome_length = len(fasta_reference[f'chr{chromosome}'].seq)  # Delete if keeping above triple comment block

        total_len += chromosome_length

        chromosome_megabase_bins_count = math.ceil(chromosome_length / bins_len)

        rate_dict[chromosome]['variation_rate'] = rate_dict[chromosome]['variation_count'] / chromosome_megabase_bins_count

        total_count += rate_dict[chromosome]['variation_count']

    total_rate = total_count / (math.ceil(total_len / bins_len))
    rate_dict['Total'] = {'variation_count': total_count, 'variation_rate': total_rate}

    logger.info('Calculated variations count and rate')

    return(rate_dict)


def check_vcf_ref_vs_fasta_ref(master_dict, vcf_data, vcf_header, fasta_file, fasta_file_compression, reduced_memory_mode):
    '''Compare the VCF file entries the FASTA file as reference.
    Count out how many variant entries (both reference and alternate) match the FASTA sequence

    Parameters
    ----------
    vcf_data : list
        list of strings containing VCF file data (each element represents a line from the VCF file)
        
    vcf_header : list
        list of strings containing headings corresponding to "vcf_data"
        
    fasta_sequences : dictionary
        keys: sequence name (usually chromosome or patch)

        values: sub-dictionary containing 
            chromosome: chromosome sequence is located on
            seq_start: sequence start position
            seq_end: sequence end position
            seq: BioPython Seq object containing sequence of entry

    Output
    ------
    mismatches : integer
        Number of instances where the VCF reference allele does not match the FASTA sequence
    
    total_vcf_refs : integer
        number of entries in VCF file
    
    percentage_mismatch : float (2 decimal places)
        percentage of "mismatches" from total number of variation entries in VCF file
    
    ref_to_alt : integer
        number of instances where VCF alt (alternate) allele matches that of the FASTA sequence
    
    percentage_switch : float (2 decimal places)
        percentage of "ref_to_alt" from total number of "mismatches"
    
    non_ref_allele : integer
        number of VCF entries that are not:
            FASTA reference matches
            VCF alt (alternate) allele matches to FASTA sequence
    '''

    logger.info('Checking VCF reference against FASTA file\'s reference alleles')

    master_dict['Total']['mismatches'] = 0
    master_dict['Total']['vcf_refs'] = 0
    master_dict['Total']['vcf_alt_match_to_fasta_ref'] = 0
    master_dict['Total']['non_ref_or_alt_variations'] = 0

    index_chrom = vcf_header.index('#CHROM')
    index_position = vcf_header.index('POS')
    index_vcf_ref_allele = vcf_header.index('REF')
    index_vcf_alt_allele = vcf_header.index('ALT')

    fasta_reference = SeqIO.index(fasta_file, 'fasta')

    previous_chromosome = None

    for line in vcf_data:

        line_split = line.split('\t')
        chromosome = line_split[index_chrom].replace('chr', '')
        position = int(line_split[index_position])
        vcf_ref_allele = line_split[index_vcf_ref_allele]
        vcf_alt_allele = line_split[index_vcf_alt_allele]

        if not 'mismatches' in master_dict[chromosome]:
            master_dict[chromosome]['mismatches'] = 0
            master_dict[chromosome]['vcf_refs'] = 0
            master_dict[chromosome]['vcf_alt_match_to_fasta_ref'] = 0
            master_dict[chromosome]['non_ref_or_alt_variations'] = 0

        if not chromosome == previous_chromosome:  # possible to reduce the need to lookup - does this reduce memory usage??
            '''if reduced_memory_mode == False:
                fasta_chromosome_sequence = fasta_reference[chromosome]['seq']
            elif reduced_memory_mode == True:
                fasta_chromosome_sequence = get_sequence_from_fasta(fasta_reference, fasta_file_compression, chromosome)'''

        previous_chromosome = chromosome

        try:
            fasta_ref_allele = fasta_reference[chromosome].seq[position - 1]  # Delete if keeping above triple comment block
        except:
            fasta_ref_allele = fasta_reference[f'chr{chromosome}'].seq[position - 1]
        # fasta_ref_allele = fasta_chromosome_sequence[position - 1]  # First index for BioPython is 0 so need to offset (negative 1)

        master_dict['Total']['vcf_refs'] += 1
        master_dict[chromosome]['vcf_refs'] += 1

        if not vcf_ref_allele == fasta_ref_allele:
            master_dict['Total']['mismatches'] += 1
            master_dict[chromosome]['mismatches'] += 1

            if vcf_alt_allele == fasta_ref_allele:
                master_dict['Total']['vcf_alt_match_to_fasta_ref'] += 1
                master_dict[chromosome]['vcf_alt_match_to_fasta_ref'] += 1
            else:
                master_dict['Total']['non_ref_or_alt_variations'] += 1
                master_dict[chromosome]['non_ref_or_alt_variations'] += 1

        fasta_ref_allele = None

    # percentage_mismatch = 100 * (master_dict['Total']['mismatches'] / master_dict['Total']['vcf_refs'])
    # percentage_switch = 100 * (master_dict['Total']['vcf_alt_match_to_fasta_ref'] / master_dict['Total']['mismatches'])
    # percentage_non_ref_allele = 100 * (master_dict['Total']['non_ref_or_alt_variations'] / master_dict['Total']['vcf_refs'])

    # mismatches = master_dict['Total']['mismatches']
    # total_vcf_refs = master_dict['Total']['vcf_refs']
    # ref_to_alt = master_dict['Total']['vcf_alt_match_to_fasta_ref']
    # non_ref_allele = master_dict['Total']['non_ref_or_alt_variations']

    # print(f'There are {mismatches:,} reference allele mismatches out of {total_vcf_refs:,} ({percentage_mismatch:.2f}%) variations in the VCF file compared to the FASTA file')
    # print(f'Of these {mismatches:,} reference allele mismatches, {ref_to_alt:,} ({percentage_switch:.2f}%) VCF alt alleles match the FASTA reference (i.e. reference switch)')
    # print(f'There are {non_ref_allele:,} variations ({percentage_non_ref_allele}) that are not VCF/FASTA reference switch alleles')

    # print('\nThe breakdown per chromosome or patch sequence is as follows:\n')

    # header = ['number_of_variations', 'number_of_mismatches', '(%)', 'VCF ref/alt switch', '(%)', 'VCF variation ref or alt does not match FASTA', '(%)']
    # print('\t'.join(header))

    for chromosome in master_dict:

        chromosome_total_vcf_entries = master_dict[chromosome]['vcf_refs']

        chromosome_mismatch = master_dict[chromosome]['mismatches']
        chromosome_percentage_mismatch = 100 * (chromosome_mismatch / chromosome_total_vcf_entries)
        master_dict[chromosome]['mismatch_percentage'] = chromosome_percentage_mismatch

        chromosome_vcf_alt_match_to_fasta_ref = master_dict[chromosome]['vcf_alt_match_to_fasta_ref']
        chromosome_percentage_vcf_ref_alt_switch = 100 * (chromosome_vcf_alt_match_to_fasta_ref / chromosome_mismatch)
        master_dict[chromosome]['vcf_ref_alt_switch_percentage'] = chromosome_percentage_vcf_ref_alt_switch

        chromosome_non_ref_variations = master_dict[chromosome]['non_ref_or_alt_variations']
        chromosome_percentage_non_ref_variations = 100 * (chromosome_non_ref_variations / chromosome_total_vcf_entries)
        master_dict[chromosome]['non_ref_variations_percentage'] = chromosome_percentage_non_ref_variations

        # print(f'{chromosome_total_vcf_entries}\t{chromosome_mismatch:,} ({chromosome_percentage_mismatch:.3f})\t{chromosome_vcf_alt_match_to_fasta_ref:,}\t({chromosome_percentage_vcf_ref_alt_switch:.3f})\t{chromosome_non_ref_variations:,}\t({chromosome_percentage_non_ref_variations:.3f})')

    data_frame = pd.DataFrame(master_dict)
    to_save = data_frame.T
    print(to_save)
    save_location = '/hps/nobackup/flicek/ensembl/variation/gurpreet/vcf_file_checker/vcf_stats.tsv'
    to_save.to_csv(save_location, sep = '\t')

    logger.info(f'Table of variation counts and rates have been saved to {save_location}')


def main():
    '''Arrange each of the other functions in this script to be called in correct order
    Only called if this script file is executing
    i.e. won't run main() if the functions are imported from another script
    '''

    vcf_file, fasta_file, reduced_memory_mode = arg_passer()
    vcf_file_compression, fasta_file_compression = check_user_input(vcf_file, fasta_file)

    # Loads VCF lines into a generator
    if reduced_memory_mode == True:
        vcf_header = get_vcf_header(vcf_file, vcf_file_compression)
        vcf_data = get_vcf_data_generator(vcf_file, vcf_file_compression)
        fasta_seqs = index_fasta(fasta_file, fasta_file_compression)
        logger.debug(f'outside type = {type(fasta_seqs)}')

    # Loads all of the data into memory - so obviously, very memory hungry
    elif reduced_memory_mode == False:
        vcf_data, vcf_header = read_vcf(vcf_file, vcf_file_compression)
        fasta_seqs = read_fasta(fasta_file, fasta_file_compression)

    master_dict = rate_of_variations(vcf_data, vcf_header, fasta_file, fasta_file_compression, reduced_memory_mode)

    # If using the generator, need to return the iterator to the beginning again
    if reduced_memory_mode == True:
        vcf_data = get_vcf_data_generator(vcf_file, vcf_file_compression)

    check_vcf_ref_vs_fasta_ref(master_dict, vcf_data, vcf_header, fasta_file, fasta_file_compression, reduced_memory_mode)


if __name__ == '__main__':

    # Start logging
    logger = LogWriter.log_writer('/homes/gurpreet/development_logging')
    LogWriter.start(logger)
    logger.info('Beginning VCF file check')

    main()

    LogWriter.end(logger)

