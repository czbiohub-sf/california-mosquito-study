#/usr/local/bin python3
"""Iterative de novo assembly with read subtraction.

contig_miner_v2
Hanna Retallack 2/17/2018
Modified from contig_miner.py (Greg Fedewag)

Requirements:
    Python 3.6.4
    magicblast 1.3.0
    PRICE Assembler v1.2
    samtools 0.1.19
    shuf

Usage:
with paired read fasta files:
python3 ~/scripts/contig_miner_v2.py -i $fasta_read1 $fasta_read2 -r $ref -t 16 &> log.cm.out

fasta_1="/data/hretallack/mosquito/fuc.CMS-032_R1_10k.fasta"
fasta_2="/data/hretallack/mosquito/fuc.CMS-032_R2_10k.fasta"
ref="/data/hretallack/databases/mos_viruses/mosquito_fullviruses_v2.fasta" #or can be an empty fasta file with single header

OR with interleaved fasta file:
python3 ~/scripts/contig_miner_v2.py -i $fasta_interleaved -r $ref -t 16 &> log.cm.out
"""

import argparse
import subprocess
import shlex
import textwrap
import sys
from glob import glob
import os
from os import system

def check_len():
    class InputAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if not 1 <= len(values) <= 2:
                raise argparse.ArgumentError('Incorrect number of inputs')
            setattr(namespace, self.dest, values)
    return InputAction


def fasta_subsample(fasta, n):
    '''Randomly subsample n reads from a FASTA'''
    #for fastq, use "awk 'NR % 2 == 0 && NR % 4 != 0 {print }' "
    subsample_command = '''awk 'NR % 2 == 0 {print }' ''' + f'{fasta} | shuf -n {n} | ' + r'''awk '{print ">"; print}' > temp_seed_reads.fasta'''
    fasta = 'temp_seed_reads.fasta'
    subprocess.run(subsample_command, shell=True)
    return fasta


def price_contigs(fasta_1, fasta_2, num_threads, loop_count):
    '''Use PRICE to build contigs '''
    print('Building contigs')
    amp_len = 400
    identity_percent = 99
    seed_reads = fasta_subsample(fasta_1, 1000)
    input_steps = 1
    second_input_cycle = 1
    constant_mult = 2
    cycle_count = 10
    sub_assembler_length = 72
    min_overlap = 30
    threshold_scaler = 20
    min_percent_id = 80
    min_len_filter = 500
    output_loop_skip = 10
    untargeted_cycles = 4
    output_name = f'loop_{loop_count}.contigs.fasta'
    return_file = f'loop_{loop_count}.contigs.cycle{cycle_count}.fasta'

    #Joe's price magic:

    # For PriceTI, my strategy is to reduce the read pool as much as possible, 
    # by removing everything you can ahead of time. Then pick 5000 seeds randomly, 
    # then run 4 cycles untargeted, followed by 10 cycles targeted. 
    # Use -lenf to limit to contigs > 500.  Then, collect contigs. Repeat ~5 times. 
    # Then, re-seed with all harvested contigs, running in targeted mode to build them out. 
    # Use -lenf to limit to 1kb.

    #PriceTI-1.2 -icf hku-23kb-seed.fasta 1 1 1 -fpp hku-hits.1.fasta hku-hits.2.fasta 300 99 
    # -mol 30 -a 12 -target 80 8 2 2 -nc 10 -lenf 350 4 -lenf 500 8 -o hku-round2.fasta

    price_command = (f'PriceTI -fpp {fasta_1} {fasta_2} {amp_len} {identity_percent} '
                     # should the below line be there? Am I using starting reads?
                     f'-icf {seed_reads} {input_steps} {second_input_cycle} {constant_mult} ' #+ \
                     f'-nc {cycle_count} '
                     f'-dbmax {sub_assembler_length} '
                     f'-mol {min_overlap} '
                     f'-tol {threshold_scaler} '
                     f'-mpi {min_percent_id} '
                     f'-lenf 350 4 -lenf 500 8 -lenf {min_len_filter} {cycle_count-1} '
                     f'-target 80 {untargeted_cycles} 2 2 ' #runs 6 cycles untargeted
                     f'-nco {output_loop_skip} '
                     f'-a {num_threads} '
                     f'-o {output_name}')
    subprocess.run(shlex.split(price_command),
                   stdout=open(f'loop_{loop_count}.price.log', 'w'))
    subprocess.run(f'rm {seed_reads}', shell=True)
    print(price_command)
    return return_file


def get_max_contig_len(contig_file):
    '''Grab max contig length from PRICE output FASTA file
    Returns 0 if fasta file is empty'''
    fasta_file = open(contig_file, 'r')
    l = 0
    for _ in fasta_file:
        if _.startswith('>contig_1'):
            l = _.split('>contig_1 (')[1].split('nt)')[0]
    return int(l)


def align_reads(fasta_1, fasta_2, reference, num_threads, loop_count):
    '''Align reads to the reference or new contigs.'''
    print('Aligning reads to reference')
    sam_file = f'loop_{loop_count}.alignment.sam'
    magicblast_command = (f'magicblast -query {fasta_1} -query_mate {fasta_2} '
                          f'-subject {reference} '
                          f'-paired '
                          #f'-num_threads {num_threads} '
                          f'-out {sam_file} ')
    subprocess.run(shlex.split(magicblast_command))
    return sam_file


def process_samfile(samfile, loop_count):
    '''Extract unmapped readpairs from samfile, write to fasta'''
    print('Extracting unmapped readpairs')
    out_R1 = f'loop_{loop_count}.unmapped.1.fasta'
    out_R2 = f'loop_{loop_count}.unmapped.2.fasta'

    #From here: https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd
    # # R1 unmapped, R2 mapped [flags 69,133...]
    # samtools view -f 4 -F 264
    # # R1 mapped, R2 unmapped [flags 89,137,153...]
    # samtools view -f 8 -F 260
    # # R1 & R2 unmapped [flags 77,141...]
    # samtools view -f 12 -F 256

    # grab all possible unmapped readpairs
    os.system(f'samtools view -f 4 -F 264 -S {samfile} > tmp.sam 2> /dev/null')
    os.system(f'samtools view -f 8 -F 260 -S {samfile} >> tmp.sam 2> /dev/null')
    os.system(f'samtools view -f 12 -F 256 -S {samfile} >> tmp.sam 2> /dev/null')
    # write to fasta file
    os.system(r'''samtools view -S tmp.sam -f 0x40 2> /dev/null | awk '{OFS="\t"; print ">"$1"\n"$10}' >''' + out_R1 )
    os.system(r'''samtools view -S tmp.sam -f 0x80 2> /dev/null | awk '{OFS="\t"; print ">"$1"\n"$10}' >''' + out_R2 )
    os.system(f'rm tmp.sam {samfile}')

    with open(out_R1, 'r') as reads:
        num_unmapped = int(sum(1 for line in reads)/2)

    return out_R1, out_R2, num_unmapped


def reconcile_contigs(fasta_file_list, read1, read2):
    '''PRICE to assemble contigs in fasta files together with original read files'''
    print('Now combining and extending contigs from all loops: ' + ' '.join(fasta_file_list))
    os.system(f'cat loop_*.contigs.cycle10.fasta > combined.contigs.fasta')
    final_price_command = (f'PriceTI -icf combined.contigs.fasta 1 1 5 '
        f'-fpp {read1} {read2} 400 99 '
        f'-target 80 4 1 1 ' #run 4 cycles untargeted, then alternate
        f'-a 12 -lenf 500 9 '
        f'-nc 10 -nco 10 -o final.contigs.fasta')
    subprocess.run(shlex.split(final_price_command),
                   stdout=open('final.price.log', 'w'))
    os.system('rm combined.contigs.fasta')
    os.system('mv final.contigs.cycle10.fasta final.contigs.fasta')


def main():
    in_fasta = args.input
    reference = args.ref
    num_threads = args.threads
    min_read_number = 100
    read_list_len = min_read_number + 1
    loop_count = 0
    unproductive_loops = 0 #counter
    max_unproductive_loops = 5 #3
    max_loops_ever = 20 #10 #temporary halt
    contig_file_list = []

    # Handle fasta input
    print(f'Input fasta file: {in_fasta}')
    print(f'Input reference file: {reference}')

    if len(in_fasta) == 2 : #assumes read pairs
        print('Assuming paired read files')
        fasta_1 = in_fasta[0]
        fasta_2 = in_fasta[1]
    elif len(in_fasta) == 1 : #assumes interleaved fasta file
        print('Assuming interleaved fasta file')
        os.system(f'cat {in_fasta[0]} | grep -A1 "/1" --no-group-separator > temp_R1.fasta')
        os.system(f'cat {in_fasta[0]} | grep -A1 "/2" --no-group-separator > temp_R2.fasta')
        fasta_1 = 'temp_R1.fasta'
        fasta_2 = 'temp_R2.fasta'

    in_fasta_1 = fasta_1 #save the original in_fasta files for final price step
    in_fasta_2 = fasta_2

    # Keep cycling until too many unproductive loops (longest contig is too short), or you don't have enough reads left
    while (read_list_len > min_read_number) and (unproductive_loops < max_unproductive_loops) and (loop_count < max_loops_ever):
        if loop_count>0:

            # ---section1--- De novo assembly of remaining unmapped reads
            try:
                built_contigs = price_contigs(fasta_1, fasta_2, num_threads, loop_count)
            except Exception as identifier:
                raise Exception('PRICE failed.')

            if os.path.exists(built_contigs) and os.path.getsize(built_contigs) > 0:
                print(f'- PRICE built at least one long contig in loop {loop_count}')
                try:
                    max_contig_len = get_max_contig_len(built_contigs)
                    print(f'Max contig length = {max_contig_len}')
                except:
                    pass
                reference = built_contigs
                contig_file_list.append(built_contigs)
            else:
                print('- PRICE output contig file does not exist or is empty')
                unproductive_loops +=1 #add to counter of unproductive 
                print(f'- Number of unproductive loops = {unproductive_loops}')
                continue           

        #---section2--- do this in first loop and after successful PRICE assemblies
        try:
            sam_mapping = align_reads(fasta_1, fasta_2, reference, num_threads, loop_count)
            unmapped_fasta_1, unmapped_fasta_2, n_unmapped_reads = process_samfile(sam_mapping, loop_count)              
        except Exception as identifier:
            raise Exception('Alignment and subtraction failed.')

        else:
            fasta_1 = unmapped_fasta_1
            fasta_2 = unmapped_fasta_2

        print(f'Loop {loop_count} complete. Unmapped reads remaining = {n_unmapped_reads}.')
        print(f'--- End of loop {loop_count} ---')
        loop_count += 1

    print('--- Loops complete ---')

    reconcile_contigs(contig_file_list, in_fasta_1, in_fasta_2) #check for redundancy in contigs and output compiled list
    os.system(f'rm loop_*.unmapped.*.fasta') #clean up unmapped reads
    os.system(f'rm temp_R*.fasta')
    os.system(f'rm final.price.log')
    #eventually, clean up other temp files too (cycle10, and price logs)

    return


if __name__ == '__main__':
    #Greg's argparser thing ... doesn't seem to output help comments?
    parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=True)
    parser_input = parser.add_mutually_exclusive_group(required=True)
    # parser_input.add_argument('-i', '--input', type=str, nargs='+', action=check_len(),
    #                           help='Inputs to start on.')
    parser_input.add_argument('-i', '--input', type=str, nargs='*',
                              help='Inputs to start on.')
    parser_input.add_argument('-f', '--folders', type=str,
                              help='The folder containing all the sample folders')
    parser.add_argument('-g', '--program', type=str, required=False,
                        choices=[], nargs='*',
                        help=textwrap.dedent('''\
                                Some help goes here.
                                '''))
    parser.add_argument('-r', '--ref', type=str, default=os.path.expanduser('~/.linuxbrew/opt/gmap-gsnap/share'),
                        help='The location, and name of the reference for alignment if you are using --input. Just the location if you are using --folders.')
    parser.add_argument('-t', '--threads', type=int, required=False, help='Number of threads to use.', default=1)
    args = parser.parse_args()
    main()















