#!/usr/bin/env python3
# Author:  Rachel Ehrlich

from subprocess import call
import os
import argparse

import pandas as pd


def make_args():
    parser = argparse.ArgumentParser(description='', add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-m', '--max_iterations',
                        default=10,
                        help='Max number of times to run pilon',
                        type=int)
    parser.add_argument('-t', '--threads',
                        help='Number of threads to use',
                        default=16,
                        type=int)

    required_flags = parser.add_argument_group('Required arguments')

    required_flags.add_argument('--r1',
                                help='Path to r1 paired reads in fastq format',
                                required=True)
    required_flags.add_argument('--r2',
                                help='Path to r2 paired reads in fastq format',
                                required=True)
    required_flags.add_argument('-u', '--unpaired',
                                help='Path to unpaired reads in fastq format',
                                required=True)
    required_flags.add_argument('-f', '--fasta_path',
                                help='Path to reference fasta',
                                required=True)
    required_flags.add_argument('-o', '--outdir',
                                help='Directory for output, must not exist',
                                required=True)
    required_flags.add_argument('-s', '--strain_name',
                                help='Nickname for the reference fasta',
                                required=True)

    return parser.parse_args()


def classify_change(source_seq, dest_seq):
    bases = set('ATCG')

    source_size = len(source_seq)
    dest_size = len(dest_seq)

    if source_seq == '.':
        variant = 'insertion'
        source_size = 0
    elif dest_seq == '.':
        variant = 'deletion'
        dest_size = 0
    elif (source_seq in bases) and (dest_seq in bases):
        variant = 'SNP'
    else:
        variant = 'complex variant'

    size_change = dest_size - source_size

    return size_change, variant


def get_ctg_and_postion(x):
    ctg, location = x.split(':')
    if '-' in location:
        start, end = location.split('-')
    else:
        start = location
        end = location
    return ctg, start, end


def parse_changes(log_path, iteration_num):
    if not os.path.isfile(log_path):
        print("Error:  Pilon did not create a changes file")
        raise FileNotFoundError

    changes = list()
    with open(log_path, 'r') as f:
        for line in f:
            source, dest, source_seq, dest_seq = line.split()
            s_ctg, s_start, s_end = get_ctg_and_postion(source)
            d_ctg, d_start, d_end = get_ctg_and_postion(dest)
            size_change, variant = classify_change(source_seq, dest_seq)
            changes.append(pd.Series({'source': source,
                                      'destination': dest,
                                      'source_ctg': s_ctg,
                                      'source_start': s_start,
                                      'source_end': s_end,
                                      'destination_ctg': d_ctg,
                                      'destination_start': d_start,
                                      'destination_end': d_end,
                                      'source_seq': source_seq,
                                      'destination_seq': dest_seq,
                                      'size_change': size_change,
                                      'variant': variant,
                                      'pilon_num': iteration_num}))
    return changes


def main():
    args = make_args()

    for input_path in [args.r1, args.r1, args.unpaired, args.fasta_path]:
        assert os.path.isfile(input_path)
     os.makedirs(args.outdir, exist_ok=True)

    last_fasta = args.fasta_path
    changes = list()
    for iteration_num in range(1, args.max_iterations + 1):
        outdir = args.outdir + '/pilon_' + str(iteration_num)

        cmd = ['bwa.sh', last_fasta, args.r1, args.r1, args.unpaired, outdir,
               str(args.threads)]
        call(cmd)
        last_fasta = outdir + '/pilon/pilon.fasta'
        if not os.path.isfile(last_fasta):
            print('\nError:  something went wrong while running pilon')
            print('If the error looks like "Exception in thread "main" '
                  'java.lang.OutOfMemoryError: Java heap space"'
                  ' you may want to edit the default_jvm_mem_opts in the file:')
            call(['which', 'pilon'])
            return

        log_path = outdir + '/pilon/pilon.changes'
        iteration_changes = parse_changes(log_path, iteration_num)

        if iteration_changes:
            changes.extend(iteration_changes)
        else:
            break

    changes = pd.concat(changes, axis=1).T.sort_values(['pilon_num', 'variant',
                                                        'size_change'])
    changes['strain'] = args.strain_name
    col_order = ['strain', 'pilon_num', 'variant', 'size_change', 'source',
                 'destination', 'source_seq', 'destination_seq', 'source_ctg',
                 'source_start', 'source_end', 'destination_ctg',
                 'destination_start', 'destination_end']
    assert sorted(col_order) == sorted(changes.columns)
    changes = changes[col_order]
    changes.to_csv(args.outdir + '/changes2.csv', index=False)


if __name__ == '__main__':
    main()
