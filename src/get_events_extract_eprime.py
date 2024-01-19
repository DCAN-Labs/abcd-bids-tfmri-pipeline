#!/usr/bin/env python3
# coding: utf-8

"""
Extract event files for ABCD-specific task-fMRI data processing
Greg Conan: gconan@umn.edu
Created: 2021-06-07
Updated: 2024-01-19 (WRAPPER_LOC fix); 2021-12-03 (the rest)
"""

# Import standard libraries
from datetime import datetime  # for seeing how long scripts take to run
from glob import glob
import matlab.engine
import os
import pandas as pd
import shutil
import sys


def find_local_imports(flg):
    """
    Ensure that this script can find its dependencies if parallel processing
    """
    try:
        parallel_flag_pos = -2 if sys.argv[-2] == flg else sys.argv.index(flg)
        sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))
    except (IndexError, ValueError):
        dirpath = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
        sys.argv.append(flg)
        sys.argv.append(dirpath)
        sys.path.append(dirpath)
    return sys.path[-1]


# Local custom imports
WRAPPER_LOC = '--wrapper-location'
SCRIPT_DIR = find_local_imports(WRAPPER_LOC)
from src.pipeline_utilities import (
    as_cli_attr, exit_with_time_info, extract_from_json, get_all_analysis_paths,
    get_pipeline_cli_argparser, get_sub_base, get_TR_and_ntpts, glob_and_copy,
    SCAN_ARG, validate_cli_args, valid_readable_dir
)
from src.run_level_1_analysis import add_lvl_1_paths_to


def main():
    start_time = datetime.now()
    cli_args = _cli()
    paths = get_all_analysis_paths(cli_args)
    for each_run in cli_args['runs']:

        # Get required paths and information about the scanner
        run_args = {'run_number': each_run, **cli_args}
        run_paths = add_lvl_1_paths_to(paths, run_args, 'abcd_extract')
        scan = get_all_scanner_info(
            run_args, run_paths,
            read_in_scanners_info_from(run_args.pop(SCAN_ARG))
        )

        # Extract event and eprime files
        extract_eprime_via_matlab(run_args, paths, each_run)
        process_abcd_extract_files(run_args, paths, each_run, scan)

    # Delete everything except event files, unless the user said to --keep-all
    if not cli_args['keep_all']:
        shutil.rmtree(paths['dir_lvl']['1'])
        shutil.rmtree(os.path.join(cli_args['output'], 'level-1'))

    # Show user how long this script took
    exit_with_time_info(start_time)  


def _cli():
    """
    :return: Dictionary of all validated command-line arguments from user
    """
    arg_names = {'events_dir', 'fd', 'keep_all', 'output', 'runs', SCAN_ARG,
                 'bids_dir', 'ses', 'spat_smooth', 'subject',  'study_dir',
                 'study_name', 'task', 'templates', 'template1', 'template2',
                 'wb_command', as_cli_attr(WRAPPER_LOC)}
    parser = get_pipeline_cli_argparser(arg_names)
    return validate_cli_args(vars(parser.parse_args()), parser, {SCAN_ARG})


def get_all_scanner_info(cli_args, paths, all_scanners):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param all_scanners: Dictionary containing scanners' info
    :return: Dictionary with 4 keys mapped to scanner info: repetition time
             ('TR'), number of timepoints ('ntpts'), 'increment' to process
             abcd extract via Matlab, and 'censor' for censoring.
    """
    start_time_QC = datetime.now()
    all_TRs = dict()
    tpts_Ns = dict()
    scanners = dict()
    msg = '{} for dtseries run {}: {}'
    for run in cli_args['runs']:
        prefix = os.path.join(paths['sub_ses']['func'],
                              get_sub_base(cli_args, run))
        dtseries = prefix + '_bold_timeseries.dtseries.nii'

        # TODO (from computing team meeting 2021-08-09) Make this file
        #   unnecessary by reading in FD and motion-regressors separately?
        motion = prefix + '_desc-filteredincludingFD_motion.tsv'
        if os.path.exists(motion):

            # Get the repetition time and number of timepoints
            all_TRs[run], tpts_Ns[run] = get_TR_and_ntpts(
                dtseries, cli_args['wb_command']
            )

            # Quality-check (QC) the ntpts and TR, then get the scanner name
            tpts_Ns = correct_type_of_value_in(tpts_Ns, run, int)
            all_TRs = correct_type_of_value_in(all_TRs, run, float)
            if not (0.799 <= all_TRs[run] <= 0.801):
                all_TRs.pop(run)
            scanners[run] = get_scanner_name(cli_args, all_scanners,
                                             run, tpts_Ns[run])

    # Tell user whether QC succeeded
        else:
            print('One of these required files does not exist:\n{}\n{}'
                  .format(dtseries, motion))
        for msg_parts in [
            ('Repetition time', run, 'OK' if run in all_TRs
             else 'MISSING OR INCOMPLETE'),
            ('Number of volumes/timepoints', run, 'OK' if run in tpts_Ns 
             else 'ERROR\nPipeline will not be completed for run {} due to '
             'error in the number of timepoints.'.format(run)),
            ('QC procedures', run, 'PASS' if run in tpts_Ns and run in all_TRs
             else 'FAIL')
        ]:
            print(msg.format(*msg_parts))
        sub_ses_task = 'subject {} session {} task {} run {}'.format(
            cli_args['subject'], cli_args['ses'], cli_args['task'], run
        )
        print('The QC process for {} took this long to run: {}'
            .format(sub_ses_task, datetime.now() - start_time_QC))
    if all(len(x.keys()) == len(cli_args['runs'])
            for x in (all_TRs, tpts_Ns)):
        print('QC procedures successfully passed')
    else:
        runs_passed_QC = set(all_TRs.keys()
                             ).intersection(set(tpts_Ns.keys()))
        if runs_passed_QC:
            print('Only including runs that passed QC: ' + str(runs_passed_QC))
            cli_args['runs'] = list(runs_passed_QC)
        else:
            sys.exit('Error: Pipeline cannot be run on {} because QC '
                     'failed.'.format(sub_ses_task))

    # Verify that all runs have the same value for scanner name, TR, and ntpts
    scan_info = all_scanners[get_dict_val_if_same(scanners, 'Scanner')]
    scan_info['TR'] = get_dict_val_if_same(all_TRs, 'Repetition time')
    scan_info['ntpts'] = get_dict_val_if_same(tpts_Ns, 'Number of timepoints')
    return scan_info


def get_dict_val_if_same(a_dict, dict_name):
    """
    Verify that all values in a dictionary are the same; if they aren't, exit
    :param a_dict: Dictionary
    :param dict_name: String naming the dictionary
    :return: Object that's the value mapped to all dictionary keys
    """
    values = set(a_dict.values())
    if len(values) != 1:
        sys.exit('Error: {} must be the same for all runs.'.format(dict_name))
    return values.pop()


def correct_type_of_value_in(a_dict, dict_key, correct_type):
    """
    :param a_dict: Dictionary
    :param dict_key: Object mapped to the value to convert to the correct type
    :param correct_type: Function returning its input, but in a specific type
    :return: a_dict, but with the value mapped to dict_key as the correct_type
    """
    try:
        a_dict[dict_key] = correct_type(a_dict[dict_key])
    except TypeError:
        a_dict.pop(dict_key)
    return a_dict


def read_in_scanners_info_from(scanners_csv):
    """
    :param scanners_csv: String, path to a real .csv file with scanners' info
    :return: Dictionary mapping scanners' name strings to dicts with their info
    """
    scanners_dict = dict()
    scanners_info = pd.read_csv(scanners_csv)

    def add_scanner_to_dict(info):
        scanners_dict[info.get('ScannerType')] = {
            'censor': info.get('NumberOfTimePointsToCensor'),
            'increment': info.get('IncrementValue'),
            'ntpts': {info.index[x].split('Timepoints_')[-1]: info.iloc[x]
                      for x in range(3, len(info.index))}
        }

    scanners_info.apply(add_scanner_to_dict, axis='columns')
    return scanners_dict


def get_scanner_name(cli_args, all_scanners, run, ntpts=None):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param all_scanners: Dictionary containing scanners' info
    :param run: Integer or string identifying which of the runs to get info for
    :param ntpts: Integer, the number of timepoints
    :return: String, the name of the scanner used for this run of this subject
    """
    # Try to get scanner name from wb_command
    scanner = None
    if ntpts: 
        for scanner_name, scanner_params in all_scanners.items():
            if scanner_params['ntpts'][cli_args['task']] == ntpts:
                scanner = scanner_name
                break

    # If wb_command fails, try to read scanner name from subject's .csv file
    if not scanner: 
        runstr = '0{}'.format(run) if int(run) < 10 else str(run)
        scan_info_json_path = os.path.join(
            cli_args['study_dir'], cli_args['subject'], cli_args['ses'],
            'func', get_sub_base(cli_args, runstr) + '_bold.json'
        )
        if os.path.exists(scan_info_json_path):
            sub_ses_task_info = extract_from_json(scan_info_json_path)
            scanner = {
                'GE': 'DV25' if 'DV25' in sub_ses_task_info['SoftwareVersions']
                else 'DV26', 'Siemens': 'S_P', 'Philips': 'S_P'
            }[sub_ses_task_info['Manufacturer']]
    return scanner


def extract_eprime_via_matlab(cli_args, paths, run):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param run: Int (1 or 2) or string ('1' or '2') defining which run this is
    """
    # Local variables: Paths to eprime file, eprime directory, and
    # abcd_extract directory
    sourcedata_dir = os.path.join(cli_args['study_dir'], 'sourcedata',
                                  cli_args['subject'], cli_args['ses'], 'func')
    eprime_file = os.path.join(sourcedata_dir,
                               get_sub_base(cli_args, '0{}'.format(run))
                               + '_bold_EventRelatedInformation.txt')
    extractdirs = {'parent': os.path.join(cli_args['wrapper_location'],
                                          'abcd_extract_eprime-master')}
    for eachdir in ('mmil_utils', 'aux'):
        extractdirs[eachdir] = os.path.join(extractdirs['parent'], eachdir)
    paths_to_add = list(extractdirs.values())
    paths_to_add.append(paths['lvl_1']['abcd_extract'])
    paths_to_add.append(sourcedata_dir)

    # Run MATLAB code to extract eprime data
    eng = matlab.engine.start_matlab()
    for eachpath in paths_to_add:
        eng.addpath(eachpath, nargout=0)
    os.chdir(extractdirs['parent'])
    extract = {'SST': eng.abcd_extract_eprime_sst,
               'MID': eng.abcd_extract_eprime_mid,
               'nback': eng.abcd_extract_eprime_nback}[cli_args['task']]
    extract(eprime_file, 'outdir', paths['lvl_1']['abcd_extract'], nargout=0)


def process_abcd_extract_files(cli_args, paths, run, scan_info): 
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param run: Whole number (as an int or a string) defining which run this is
    """
    # Find remaining abcd_extract files and process them in event files dir
    os.system('find %s -maxdepth 1 -type f -empty -execdir mv {} %s \;'
              % (paths['lvl_1']['abcd_extract'], paths['lvl_1']['EV']))
    os.system('find %s -maxdepth 1 -type f -empty -exec bash -c \'echo \"0 '
              '0 0\" > {}\' \';\' ' % (paths['lvl_1']['EV']))

    # Move FSL scan data files in abcd_extract to intermediate_files
    types_of_files = '*scan{}*_fsl.txt'.format(run)
    glob_and_copy(paths['lvl_1']['intermediate'],
                  paths['lvl_1']['abcd_extract'], types_of_files)

    # Adjust FSL data
    for fname in glob(os.path.join(paths['lvl_1']['intermediate'],
                                   types_of_files)):
        output = os.path.join(paths['lvl_1']['EV'], os.path.basename(fname))
        df = pd.read_csv(os.path.join(paths['lvl_1']['intermediate'], fname),
                         sep='\s', engine='python', names=['A', 'B', 'C'])
        df.A = df.A + scan_info['increment']
        df.to_csv(output, sep='\t', encoding='utf-8', header=None, index=False)
            
    # Save out adjusted FSL data and copy to event files directory
    os.chdir(paths['lvl_1']['EV'])
    output_file_paths = list()
    for name in glob('*_fsl.txt'):
        f = name.replace(
                    'run-01_bold_EventRelatedInformation_scan1', 'run-1'
                ).replace(
                    'run-02_bold_EventRelatedInformation_scan2', 'run-2'
                ).replace(
                    '_fsl.txt', '.tsv'
                )
        output_file_paths.append(f)
        os.rename(name, f)
    copy_event_files_to_default_dir(cli_args, output_file_paths)


def copy_event_files_to_default_dir(cli_args, all_event_files):
    """
    Copy all event files into the default event files directory
    :param cli_args: Dictionary containing all command-line arguments from user
    :param all_event_files: List of strings that are valid paths to real files
    """
    for each_EV_file in all_event_files:
        try:
            shutil.copy(each_EV_file, cli_args['events_dir'])
        except shutil.SameFileError:
            pass


if __name__ == '__main__':
    main()
