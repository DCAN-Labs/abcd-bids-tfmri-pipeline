#!/usr/bin/env python
# coding: utf-8

"""
Original Author: Anthony Juliano, acjulian@uvm.edu
Wrapper Author: Greg Conan, gconan@umn.edu
Wrapper Created: 2021-06-22
Wrapper Updated: 2021-11-12
"""

# Import standard libraries
from glob import glob
import os
import pandas as pd
import shutil
import subprocess
import sys
from datetime import datetime


def find_myself(flg):
    """
    Ensure that this script can find its dependencies if parallel processing
    :param flg: String in sys.argv immediately preceding the path to return
    :return: String, path to the directory this Python script is in
    """
    try:
        parallel_flag_pos = -2 if sys.argv[-2] == flg else sys.argv.index(flg)
        sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))
    except (IndexError, ValueError):
        sys.path.append(os.path.dirname(os.path.dirname(
                        os.path.abspath(__file__))))
    return sys.path[-1]


# Local custom imports
WRAPPER_LOC = '--wrapper-location'
SCRIPT_DIR = find_myself(WRAPPER_LOC)
from src.get_events_extract_eprime import (
    get_all_scanner_info, process_abcd_extract_files,
    read_in_scanners_info_from,
)
from src.pipeline_utilities import (
    count_lines_in_txt_file, exit_with_time_info, get_all_analysis_paths,
    get_lvl_paths, get_main_pipeline_arg_names, get_optional_cli_args, 
    get_replacements, get_pipeline_cli_argparser, get_sub_base, SCAN_ARG,
    rename_template_file_vars, validate_cli_args, wb_command
)


def main():
    # Time how long the script takes
    start_time = datetime.now()

    # Get necessary variables used throughout the script
    cli_args, paths, scanner_info = get_necessary_dicts()

    for run in cli_args['runs']:  # Create EV files
        make_EV_files_for_run(cli_args, paths, scanner_info, run)
    print('Generating EV files completed.')

    if cli_args['pipeline']:  # Run main pipeline
        print('Continuing onto pipeline.')
        args_to_run = [
            'python3', paths['wrapper'], *get_optional_cli_args(cli_args)
        ]
        print(args_to_run)
        subprocess.check_call(args_to_run)
        renaming_moving_and_cleanup(cli_args, paths)
    exit_with_time_info(start_time)


def get_necessary_dicts():
    """
    Collect and return 3 dictionaries necessary to run this script.
    :return: Tuple with the following elements in order:
        1. Dictionary containing all validated command-line arguments,
        2. Dictionary of additional paths, and
        3. Dictionary of scanner information (to get the number of timepoints).
    """
    cli_args = _cli()
    paths = get_all_analysis_paths(cli_args)

    # More necessary paths: full path to main pipeline script, string uniquely
    # identifying this subject+session+task, full path to abcd_extract_eprime
    # dir, and python wrapper template used to run abcd_extract_eprime scripts
    paths['wrapper'] = os.path.join(SCRIPT_DIR, 'pipeline_wrapper.py')
    paths['sub_base'] = get_sub_base(cli_args)
    paths['extract_eprime'] = {
        'template': os.path.join(
            cli_args['templates'], 'template_abcd_extract_eprime_wrapper.txt'
        ), 'master': os.path.join(SCRIPT_DIR, 'abcd_extract_eprime-master')
    }
    return cli_args, paths, get_all_scanner_info(
        cli_args, paths, read_in_scanners_info_from(cli_args.pop(SCAN_ARG))
    )


def _cli():
    """
    :return: Dictionary containing all command-line arguments from user
    """
    arg_names = get_main_pipeline_arg_names().union({SCAN_ARG, })
    parser = get_pipeline_cli_argparser(arg_names)
    parser.add_argument(
        '--pipeline', '-pipeline', action='store_true',
        help=('Include this flag to run the main pipeline script immediately '
              'after generating the all_trials.tsv files.')
    )
    return validate_cli_args(vars(parser.parse_args()), parser, arg_names)


def make_EV_files_for_run(cli_args, paths, scanner_info, run):
    """
    Create all event files for 1 run of 1 session of 1 subject doing 1 task.
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param scanner_info: Dictionary of scanner info including how many timepoints
    :param run: Whole number (as an int or a string) defining which run this is
    """
    sub_run_basename = get_sub_base(cli_args, run)
    paths['lvl_1'] = get_lvl_paths(
        paths['dir_lvl']['1'], # os.path.join(cli_args['output'], 'level-1'),
        sub_run_basename, paths['feat_name'],
        list(), 'abcd_extract', 'fsf', 'EV', 'intermediate'
    )
    for eachdir in paths['lvl_1'].values():
        os.makedirs(eachdir, exist_ok=True)

    eprime_path = os.path.join(cli_args['study_dir'], 'sourcedata',
                               cli_args['subject'], cli_args['ses'], 'func')  
    paths['extract_eprime']['python'] = os.path.join(
        paths['lvl_1']['abcd_extract'], 'abcd_extract_eprime_wrapper.py'
    )

    # Copy .fsf template file, changing placeholder strings to specific values
    rename_template_file_vars(
        paths['extract_eprime']['template'], paths['extract_eprime']['python'],
        get_replacements(cli_args, **{
            'EPRIME_FILE': os.path.join(eprime_path, get_sub_base(
                    cli_args, '0{}'.format(run)
                ) + '_bold_EventRelatedInformation.txt'),
            'OUTPUTS': paths['lvl_1']['abcd_extract'], 
            'MATLAB_SCRIPT': 'abcd_extract_eprime_' + cli_args['task'].lower(),
            'MMIL': os.path.join(paths['extract_eprime']['master'],
                                 'mmil_utils'),
            'AUX': os.path.join(paths['extract_eprime']['master'], 'aux'),
            'EPRIME_PATH': eprime_path, 
            'EXTRACT_EPRIME': paths['extract_eprime']['master']
        })
    )

    # Run script to get event .tsv files for this run/task/subject/session
    subprocess.check_call(('python3', paths['extract_eprime']['python']))
    clean_up_abcd_extract_files(cli_args, paths, sub_run_basename)
    process_abcd_extract_files(cli_args, paths, run, scanner_info) #, rm_EV=False)

    # Read in EV files, excluding the empty ones which just contain '0 0 0'
    basename_prefix = os.path.join(paths['lvl_1']['EV'], sub_run_basename)
    dfs_EVs = list()
    for eachtsv in glob(basename_prefix + '*.tsv'):
        is_valid_EV = True
        if count_lines_in_txt_file(eachtsv) == 1:
            with open(eachtsv) as infile:
                contents = infile.read().strip().split()
                if all(x == '0' for x in contents):
                    is_valid_EV = False
        if is_valid_EV:
            dfs_EVs.append(pd.read_csv(eachtsv, sep='\s+|\t+', engine='python',
                                       names=['onset', 'duration', 'weight']))

    # Concatenate all EV files into a single file consisting of all trials
    all_trials_file = basename_prefix + '_all_trials.tsv'
    all_trials = dfs_EVs[0].append(dfs_EVs[1:], ignore_index=True)
    all_trials = all_trials.sort_values('onset', axis=0).to_csv(
        all_trials_file, sep='\t', index=False, header=False
    )
    shutil.copy2(all_trials_file, os.path.join(
        cli_args['events_dir'], os.path.basename(all_trials_file)
    ))
    print('{}uccessfully saved to {}'.format(
        'S' if os.path.exists(all_trials_file) else 'Uns', all_trials_file
    ))


def clean_up_abcd_extract_files(cli_args, paths, sub_run_basename):
    """
    Delete abcd_extract temp scan files (but not all files there)
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    """
    to_delete = list()
    pairs_runs = {1: '2', 2: '1'}
    for run in cli_args['runs']:
        to_delete += glob(os.path.join(
            paths['lvl_1']['parent'], sub_run_basename + str(run),
            'abcd_extract_files', '*scan{}*'.format(pairs_runs[run])
        )) 
    for eachfile in to_delete:
        if os.path.exists(eachfile):
            os.remove(eachfile)


def renaming_moving_and_cleanup(cli_args, paths):
    """
    Rename and move essential gfeat files, then delete the rest
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    """
    # Make level 2 processing subdirectories 
    paths['lvl_2'] = get_lvl_paths(
        paths['output'], paths['sub_base'],
        cli_args['study_name'] + '.gfeat', cli_args['runs'],
        'cope', 'L_gii', 'R_gii', 'subcort'
    )
    for eachdir in paths['lvl_2'].values():
        os.makedirs(eachdir, exist_ok=True)

    # Find all generated output files to move/delete
    sub_cope_out = os.path.join(cli_args['output'],
                                cli_args['study_name'] + '_gfeat_outputs',
                                'cope_files')
    if os.path.exists(os.path.join(sub_cope_out, paths['sub_base']
                                   + '_contrast_1_cope.dtseries.nii')):
        contrasts = glob('%s/%s_contrast_*_cope.dtseries.nii'
                         %(sub_cope_out, paths['sub_base']))
        num_contrasts = len(contrasts) + 1
        dts = '.dtseries.nii'

        # Split every contrast into subcortical and left/right hemisphere files
        for contrast_num in range(1, num_contrasts):
            get_contrast = get_generic_contrast_function(paths['sub_base'],
                                                         contrast_num)
            os.rename(get_contrast(sub_cope_out, dts),
                      get_contrast(paths['lvl_2']['cope'], dts))
            wb_command(cli_args, ' -cifti-separate-all',
                       get_contrast(paths['lvl_2']['cope'], dts), '-volume',
                       get_contrast(paths['lvl_2']['subcort'], '.nii'), '-left',
                       get_contrast(paths['lvl_2']['L_gii'], '_L.func.gii'),
                       '-right', get_contrast(paths['lvl_2']['R_gii'], '_R.func.gii'))
            
        # Delete all remaining output files
        for d in (paths['lvl_2'][1], paths['lvl_2'][2],
                  cli_args['output'], paths['lvl_1']['parent']):
            if os.path.exists(d):
                shutil.rmtree(d)


def get_generic_contrast_function(sub_basename, contrast_num):
    """
    :param sub_basename: String identifying 1 subject, task, and session
    :param contrast_num: Int identifying the contrast/cope number
    :return: Function which accepts 2 strings (one path to a parent directory
             and one file extension) and returns a string, the full path to a 
             contrast file with its parent directory
    """
    return lambda parent, ext: os.path.join(
        parent, sub_basename + '_contrast_{}_cope'.format(contrast_num) + ext
    )


if __name__ == '__main__':
    main()
