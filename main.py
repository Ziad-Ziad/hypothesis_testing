"""
This script performs various statistical tests using data from input files.
"""

import argparse
import numpy as np
import two_sample

DESC = '''The data of ['z_test', 't_independent', 'paired_test', 'f_test'] has to be stored in
control.txt and variation.txt | The data of the chi2_conting, data has be stored in table.txt file
| and the data of proportion_test, in seperate lines in data.txt, write the values of obs1, n1, 
obs2, n2 in order.
'''
parser = argparse.ArgumentParser(description=DESC)
parser.add_argument('test', type=str, choices=['z_test', 't_independent', 'proportion_test', \
                                               'paired_test', 'f_test', 'chi2_conting'])               
parser.add_argument('--alpha', type=float, default=0.05, help="The siginificance level")
parser.add_argument('--alter', type=str, default='two-tailed', choices=['two-tailed', 'right', 'left'],
                    help="Currently doesnt support proportion_test, f_test,chi2_conting tests")

args = parser.parse_args()
test_stat = args.test
sig_level = args.alpha
alternative = args.alter


def read_data(file_name) -> np.ndarray:
    """
    Read floating-point data from a file.

    Args:
        file_name: The name of the file to read from.

    Returns:
        numpy.ndarray: An array containing the read floating-point data.
    """
    with open(file_name, 'r', encoding='utf-8') as file:
        data = np.array([int(line.strip()) for line in file], dtype=np.float32)
    return data

def read_table(file_name) -> np.ndarray:
    """
    Read floating-point data from a table in a file.

    Args:
        file_name: The name of the file to read from.

    Returns:
        numpy.ndarray: An array containing the read floating-point data.
    """
    with open(file_name, 'r', encoding='utf-8') as file:
        data = np.array([line.split() for line in file], dtype=np.float32)
    return data

dict_tests = {
    'z_test': two_sample.z_test,
    't_independent': two_sample.t_independent,
    'proportion_test': two_sample.proportion_test,
    'paired_test': two_sample.paired_test,
    'f_test': two_sample.f_test, 
    'chi2_conting': two_sample.chi2_conting
}

perform_test_stat = dict_tests[test_stat]

if test_stat == 'chi2_conting':
    table = read_table('table.txt')
    print(perform_test_stat(table, alpha=sig_level))

elif test_stat == 'proportion_test':
    obs_ata = read_data('data.txt')
    print(perform_test_stat(obs_ata[0], obs_ata[1], obs_ata[2], obs_ata[3], alpha=sig_level))
else:
    control = read_data('control.txt')
    variation = read_data('variation.txt')
    if test_stat == 'f_test':
        print(perform_test_stat(control, variation, alpha=sig_level))
    else:
        print(perform_test_stat(control, variation, alpha=sig_level, alter=alternative))