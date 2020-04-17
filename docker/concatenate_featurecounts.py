#! /usr/bin/python3

import pandas as pd
import argparse
import os

def parse_args():
    '''
    Responsible for parsing the input args
    Require the output filename and an array of >0
    files to concatenate.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', \
        required=True, \
        dest = 'output_path',
        help='Path for the concatenated output matrix.'
    )
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()
    return args

def cat_tables(input_files):
    '''
    Concatenates the count files into a raw count matrix.
    Logic is specific to the format of the featureCounts output
    files.
    '''
    count_matrix = pd.DataFrame()
    for f in input_files:
        samplename = '.'.join(os.path.basename(f).split('.')[:-3])
        df = pd.read_csv(f, sep='\t', comment='#', index_col=0)
        s = df.iloc[:,-1]
        s.name = samplename
        count_matrix = pd.concat([count_matrix, s], axis=1, sort=True)
    return count_matrix

if __name__ == '__main__':
    args = parse_args()
    count_matrix = cat_tables(args.input_files)
    count_matrix.to_csv(args.output_path, sep='\t', index_label='Gene')
