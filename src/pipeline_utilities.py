#!/usr/bin/env python3
# coding: utf-8

"""
Common source for utility functions used by ABCD-BIDS task-fmri-pipeline
Greg Conan: gconan@umn.edu
Created: 2021-01-15
Updated: 2024-01-19 (as_cli_attr); 2021-11-12 (the rest)
"""

# Import standard libraries
import argparse
from datetime import datetime  # for seeing how long scripts take to run
from glob import glob
import json
import multiprocessing as mp
import os
import pandas as pd
import random  # only used by rand_string
import shutil
import string  # only used by rand_string
import subprocess
import sys
import time

# Constants: Name of scanner-info command-line argument, directory containing
# the main pipeline script, SLURM-/SBATCH-related arguments' default names, and
# name of the argument to get the directory containing the main wrapper script
SCAN_ARG = 'scanners_info'
SCRIPT_DIR = os.path.dirname(os.path.dirname(__file__))
SLURM_ARGS = ('account', 'cpus', 'memory', 'print_progress', 'sleep', 'time')
WRAPPER_LOC = 'wrapper_location'


def add_arg_if_in_arg_names(arg_name, all_args, parser, *shortnames, **kwargs):
    """
    Wrapper for argparse.ArgumentParser.add_argument. Nearly identical, but 
    will only add the argument to the parser if arg_name is in all_args.
    :param arg_name: String naming the argument to (maybe) add to parser
    :param all_args: Set of strings; each names a command-line argument
    :param parser: argparse.ArgumentParser
    :param shortnames: Unpacked list of strings; each is arg_name shortened
    :param kwargs: Unpacked dictionary of argparse attributes to give the arg
    :return: parser, but (maybe) with the argument named arg_name added
    """
    if arg_name in all_args:
        cli_arg = as_cli_arg(arg_name)
        parser.add_argument(
            cli_arg[1:], cli_arg, *shortnames, **kwargs
        )
    return parser


def add_lvl_args_to(parser):
    """
    :param parser: argparse.ArgumentParser with all command-line arguments 
                   that the user gave to pipeline_wrapper.py
    :return: parser with all command-line arguments needed for level X analysis
    """
    # 1) Top-level directory with pipeline_wrapper.py 2) Run number 3) Path to
    # .json file which stores the 'paths' dictionary 
    # parser.add_argument('--code-dir', type=valid_readable_dir, required=True)
    parser.add_argument('--run-number', type=valid_whole_number, required=True)
    parser.add_argument('--temp-json', type=valid_readable_json, required=True)
    return parser


def add_slurm_args_to(parser):
    """
    :param parser: argparse.ArgumentParser with some command-line arguments 
    :return: parser with all CLI arguments needed to run parallel SLURM jobs
    """
    default_CPUs = 1
    default_gb_mem = 8
    default_sleep = 10
    default_time_limit = "01:00:00"
    parser.add_argument(
        '-A', '--account',
        help="Name of the account to submit the SBATCH job under."
    )
    parser.add_argument(
        '-c', '--cpus', type=valid_whole_number, default=default_CPUs,
        help=('Number of CPUs to use for each Python job. By default, this '
              'argument\'s value will be {}.'.format(default_CPUs))
    )
    parser.add_argument(
        '-mem', '--memory', type=valid_whole_number, default=default_gb_mem,
        help=("Memory in gigabytes (GB) to assign to each sbatch job. The "
              "default number is {} GB.".format(default_gb_mem))
    )
    parser.add_argument(
        '-progress', '--print-progress', action='store_true',
        help=('Include this flag for the script to print updates about its '
              'progress at intervals defined by --sleep. This will also print '
              'every command that is run to submit a pipeline batch job.')
    )
    parser.add_argument(
        '-sleep', '--sleep', type=valid_whole_number, default=default_sleep,
        help=("Number of seconds to wait between batch job submissions. The "
              "default number is {}.".format(default_sleep))
    )
    parser.add_argument(
        '-time', '--time', metavar="SLURM_JOB_TIME_LIMIT",
        type=valid_time_str, default=default_time_limit,
        help=("Time limit for each automated_subset_analysis batch job. The "
              "time limit must be formatted specifically as HH:MM:SS where HH "
              "is hours, MM is minutes, and SS is seconds. {} is the default "
              "time limit.".format(default_time_limit))
    )
    return parser


def argify(argname, argval):
    """
    :param argname: String naming a parameter for a script called from terminal
    :param argval: Object to assign in string form as the value of the argument
    :return: String, a parameter assignment for a script called from terminal
    """
    return "--{}={}".format(argname, argval)


def as_cli_arg(arg_str):
    """
    :param arg_str: String naming a stored argument taken from the command line
    :return: String which is the command-line argument form of arg_str
    """
    return "--" + arg_str.replace("_", "-")


def as_cli_attr(arg_str):
    """
    :param arg_str: String which is the command-line argument form of arg_str
    :return: String naming a stored argument taken from the command line
    """
    return arg_str.strip("--").replace("-", "_")


def copy_and_rename_file(old_file, new_file):
    """
    Rename a file and copy it to a new location
    :param old_file: String, valid path to an existing file to copy
    :param new_file: String, valid path to what will be a copy of old_file
    """
    os.rename(shutil.copy2(old_file, os.path.dirname(new_file)), new_file)


def copy_event_files_to_default_dir(cli_args, all_event_files):
    """
    Copy all event files into the default event files directory
    :param cli_args: Dictionary containing all command-line arguments from user
    :param all_event_files: List of strings that are valid paths to real files
    """
    for each_EV_file in all_event_files:
        try: shutil.copy(each_EV_file, cli_args['events_dir'])
        except shutil.SameFileError: pass


def count_lines_in_txt_file(filepath):
    """
    Quickly count how many lines are in a text file.
    Taken from pynative.com/python-count-number-of-lines-in-file
    :param filepath: String, valid path to an existing readable text file
    :return: Int, the number of lines in the file at filepath
    """
    with open(filepath, 'r') as infile:  # open file in read mode
        for count, _ in enumerate(infile):
            pass
    return count + 1


def dict_has(a_dict, a_key):
    """
    :param a_dict: Dictionary (any)
    :param a_key: Object (any)
    :return: True if and only if a_key is mapped to something truthy in a_dict
    """
    return a_key in a_dict and a_dict[a_key]


def ensure_dict_has(a_dict, a_key, new_value):
    """
    :param a_dict: Dictionary (any)
    :param a_key: Object which will be a key in a_dict
    :param new_value: Object to become the value mapped to a_key in a_dict
                      unless a_key is already mapped to a value
    :return: a_dict, but with a_key mapped to some value
    """
    if not dict_has(a_dict, a_key):
        a_dict[a_key] = new_value
    return a_dict


def exit_with_time_info(start_time, exit_code=0):
    """
    Terminate the pipeline after displaying a message showing how long it ran
    :param start_time: datetime.datetime object of when the script started
    """
    print('The pipeline for this subject took this long to run {}: {}'
          .format('successfully' if exit_code == 0 else 'and then crashed',
                  datetime.now() - start_time))
    sys.exit(exit_code)


def extract_from_json(json_path):
    """
    :param json_path: String, a valid path to a real readable .json file
    :return: Dictionary, the contents of the file at json_path
    """
    with open(json_path, 'r') as infile:
        return json.load(infile)


def get_all_analysis_paths(cli_args):
    """
    Build and save paths for various variables called throughout the pipeline
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: Dictionary containing paths to all of the following variables:
             AROI2, BIDS, dir_lvl, feat_name, final_smooth, lvl_2_paths,
             sub_paths, templates
    """
    paths = {'dir_lvl': {str(lvl): os.path.join(  # Feature dirs for all levels
                              cli_args['output'], 'Level{}_feats'.format(lvl)
                         ) for lvl in (1, 2)},
             'feat_name': '{}.feat'.format(cli_args['study_name']),
             'final_smooth': ('_smoothed_{}mm'  # Spatial smoothing variable 
                              .format(cli_args['spat_smooth']))}
    for lvl in cli_args['levels']:
        tmpllv = 'template{}'.format(lvl)
        paths[tmpllv] = os.path.join(cli_args['templates'], cli_args[tmpllv])
    paths['lvl_2'] = get_lvl_paths(
        paths['dir_lvl']['2'], get_sub_base(cli_args),
        cli_args['study_name'] + '.gfeat', cli_args['runs'], 'fsf'
    )
    paths['sub_ses'] = {f_or_a: os.path.join( # Subject anat & func directories
                            cli_args['study_dir'], 'derivatives',
                            cli_args['bids_dir'], cli_args['subject'],
                            cli_args['ses'], f_or_a
                        ) for f_or_a in ('anat', 'func')} 
    paths['AROI2'] = os.path.join(cli_args['templates'], 'Atlas_ROIs.2.nii.gz')
    return paths


def get_and_print_time_since(event_name, event_time):
    """
    Print and return a string showing how much time has passed since the
    current running script reached a certain part of its process
    :param event_name: String to print after 'Time elapsed since '
    :param event_time: datetime object representing a time in the past
    :return: String with an easily human-readable message showing how much time
             has passed since {event_time} when {event_name} happened.
    """
    timestamp = ("\nTime elapsed since {}: {}"
                 .format(event_name, datetime.now() - event_time))
    print(timestamp)
    return timestamp


def get_args_to_run_film_gls(**kwargs):
    """
    :return: List of strings which are a Bash command calling film_gls
    """
    in_arg = kwargs.pop('in_arg')
    to_call = ['film_gls', '--sa', argify('in', in_arg)] 
    for argname, argval in kwargs.items():
        to_call.append(argify(argname, argval))
    return to_call


def get_default_ext_command(cmd_name):
    """
    Try to get valid path to external software command file without user input
    :param cmd_name: String naming the executable command file
    :return: String, path to the command if the user has the command alias in
             their .bashrc / $PATH; otherwise None
    """
    try:  # If the command path is already defined, then use it
        cmd = subprocess.check_output(("which", cmd_name)
                                      ).decode('utf-8').split()[-1]
    except subprocess.CalledProcessError:
        cmd = None
    return cmd


def get_LR_functions(cli_args, paths):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :return: Dictionary mapping 'surf' to a function which returns the file
             path string to a .surf.gii file, and mapping 'shape' to a function
             which returns the file path string to a .shape.gii file
    """
    return {'surf': lambda x: os.path.join(
                paths['sub_ses']['anat'], get_subj_ses(cli_args) +
                '_hemi-{}_space-MNI_mesh-fsLR32k_midthickness.surf.gii'.format(x)
            ), 'shape': lambda y: os.path.join(
                cli_args['templates'], y + '.atlasroi.32k_fs_LR.shape.gii'
            )}


def get_lvl_paths(lvl_dir, sub_base, feat_name, runs, *extra_subdirs):
    """
    Get a dictionary of paths to analysis-level-specific files for paths dict
    :param lvl_dir: String, path to the feat directory for level 1 or 2
    :param sub_base: String identifying a subject, session, and task
    :param feat_name: String naming a feature
    :param runs: List of strings or integers, each identifying a run
    :param extra_subdirs: Unpacked list of strings naming subdirectories of 
                          the level parent directory
    :return: Dictionary mapping string keys to string paths
    """
    lvl_paths = {'parent': os.path.join(lvl_dir, sub_base + '_' + feat_name)}
    for run in runs:
        lvl_paths[run] = os.path.join(lvl_paths['parent'],
                                      'level1_run-{}'.format(run))
    for subdr in extra_subdirs:
        lvl_paths[subdr] = os.path.join(lvl_paths['parent'], subdr + '_files')
    return lvl_paths


def get_main_pipeline_arg_names():
    """
    :return: Set containing strings naming all command-line arguments included
             by default in the main script, pipeline_wrapper.py
    """
    return {'bids_dir', 'censor', 'events_dir', 'fd', 'filter', 'fsl_dir',
            'keep_all', 'levels', 'no_parallel', 'output', 'runs', 'ses',
            'spat_smooth', 'subject', 'surf_smooth', 'study_dir', 'study_name',
            'task', 'temp_dir', 'templates', 'template1', 'template2',
            'vol_smooth', 'wb_command', WRAPPER_LOC}


def get_optional_cli_args(cli_args, drop_slurm=False):
    """
    :param cli_args: Dictionary with all validated command-line arguments,
                     all of which are used by this function
    :param drop_slurm: True to exclude SLURM arguments; else False
    :return: List of most cli_args optional arguments and their values
    """
    optional_args = list()
    for arg in cli_args.keys():
        if cli_args[arg] and not (drop_slurm and arg in SLURM_ARGS):
            optional_args.append(as_cli_arg(arg))
            if isinstance(cli_args[arg], list):
                for el in cli_args[arg]:
                    optional_args.append(str(el))
            elif not isinstance(cli_args[arg], bool):
                optional_args.append(str(cli_args[arg]))
    return optional_args
    

def get_pipeline_cli_argparser(arg_names=get_main_pipeline_arg_names()):
    """
    :param arg_names: Set containing strings naming all command-line arguments
    :return: argparse.ArgumentParser with all command-line arguments 
             needed to run pipeline_wrapper.py
    """
    # Default values for user input arguments
    default_BIDS_dir = 'abcd-hcp-pipeline'
    default_censor_num = 0 # 2
    default_fd = 0.9
    default_smooth = 0
    default_study_name = 'ABCD'
    default_runs_lvls = [1, 2]
    default_temporal_filter = 100
    default_wb_command = get_default_ext_command('wb_command')
    generic_dtseries_path = os.path.join(
        '(--study-dir)', 'derivatives', '(--bids-dir)',
        '(--subject)', '(--ses)', 'func',
        'sub-(--subject)_ses-(--ses)_task-(--task)_'
        'run-(--runs)_bold_timeseries.dtseries.nii'
    )
    generic_output_dirpath = os.path.join('(--study-dir)', 'derivatives',
                                          'abcd-bids-tfmri-pipeline',
                                          '(--subject)', '(--ses)')

    # Strings used in multiple help messages
    msg_default = ' By default, this argument\'s value(s) will be {}.'
    msg_pipeline = 'Name of the {} that you are running the pipeline on.'
    msg_smooth = ('Millimeters of {} smoothing that has already been applied '
                  'in the minimal processing steps.')
    msg_template = 'Name (not full path) of the Level {} .fsf template file.'
    msg_whole_num = ' This argument must be a positive integer.'

    # Create parser with command-line arguments from user
    parser = argparse.ArgumentParser(description=(
        'ABCD fMRI Task Prep pipeline. Inputs must be in the same format '
        'as ABCD-HCP-Pipeline outputs after running filemap.'
    ))
    parser = add_arg_if_in_arg_names('bids_dir', arg_names, parser,
        metavar='NAME_OF_BIDS_DERIVATIVES_PIPELINE_DIRECTORY',
        default=default_BIDS_dir,
        help=('Name of the BIDS-standard file-mapped directory with subject '
              'data in the "derivatives" subdirectory of your --study-dir. '
              'This path should be valid: ' + generic_dtseries_path +
              msg_default.format(default_BIDS_dir))
    )
    # Specify how many initial frames/volumes to censor
    parser = add_arg_if_in_arg_names('censor', arg_names, parser,
        metavar='INITIAL_NUMER_OF_TIMEPOINTS_TO_CENSOR', 
        default=default_censor_num, type=valid_whole_number,
        help=('The number of initial frames/volumes to censor.'
              + msg_whole_num + msg_default.format(default_censor_num))
    )
    parser = add_arg_if_in_arg_names('events_dir', arg_names, parser, 
        metavar='EVENT_FILES_DIRECTORY',
        type=valid_readable_dir, 
        help='Valid path to a real directory containing event .tsv files.' 
    )
    # Specify framewise displacement threshold to censor volumes with high motion
    parser = add_arg_if_in_arg_names('fd', arg_names, parser, 
        metavar='FRAMEWISE_DISPLACEMENT_THRESHOLD',
        default=default_fd, type=valid_float_0_to_1,
        help=('The framewise displace threshold for censoring volumes with '
              'high motion. This must be a decimal between 0 and 1.{}'
              .format(msg_default.format(default_fd)))
    )
    # High pass temporal filter cutoff number value
    parser = add_arg_if_in_arg_names('filter', arg_names, parser,  
        metavar='HIGH_PASS_TEMPORAL_FILTER_CUTOFF',
        default=default_temporal_filter, type=valid_whole_number,
        help=('High pass filter cutoff (in seconds).{}{}'.format(
            msg_whole_num, msg_default.format(default_temporal_filter)
        ))
    )
    parser = add_arg_if_in_arg_names('fsl_dir', arg_names, parser, 
        '-fsl', '--fsl', dest='fsl_dir', type=valid_readable_dir,
        help=('Valid path to an existing directory containing the executable '
              'files fsl, fslmerge, fslmaths, flameo, and feat_model from '
              'the FMRIB Software Library (FSL).')
    )
    parser = add_arg_if_in_arg_names('keep_all', arg_names, parser, 
        action='store_true',
        help=('Include this flag to keep all files generated during the '
              'pipeline. By default, the pipeline will only keep dtseries, '
              'dof, log, and event files.')
    )
    # Which analysis levels to run
    parser = add_arg_if_in_arg_names('levels', arg_names, parser,  
        metavar='ANALYSIS_LEVELS_TO_RUN',
        nargs='*', choices=default_runs_lvls, type=valid_whole_number,
        help=('Levels to conduct the analysis on: {0} for one run, and/or '
              '{1} to merge multiple runs.'.format(*default_runs_lvls))
    )
    parser = add_arg_if_in_arg_names('no_parallel', arg_names, parser,
        action='store_true',
        help=('Include this flag to process level 1 analysis runs '
              'sequentially. By default, the script will process the analyses '
              'in parallel simultaneously.')
    )
    parser = add_arg_if_in_arg_names('output', arg_names, parser, 
        '-out', metavar='OUTPUT_DIRECTORY', type=valid_output_dir, # required=True, 
        help=('Directory path to save pipeline outputs into.'
              + msg_default.format(generic_output_dirpath))
    )
    # Specify the number of runs each subject has
    parser = add_arg_if_in_arg_names('runs', arg_names, parser,  
        metavar='RUN', 
        default=default_runs_lvls, type=valid_whole_number, nargs="+",
        help=('Each subject\'s number of runs. This argument must be 1 or '
              'more positive integers provided as a space-delimited list. '
              'For example: 1 2 3 4. By default, this argument\'s value(s) '
              'will be 1 2.')
    )
    parser = add_arg_if_in_arg_names(SCAN_ARG, arg_names, parser, 
        type=valid_readable_file,
        help=('Path to existing .csv file listing all scanners\' parameters. '
              + msg_default.format('scan_info/{}.csv in the code directory.'
                                   .format(SCAN_ARG)))
    )
    # Which session to run the pipeline on
    parser = add_arg_if_in_arg_names('ses', arg_names, parser,  
        metavar='SESSION', required=True, # default=default_ses,
        type=lambda x: valid_subj_ses(x, 'ses-', 'session'), #, 'ses'),
        help=msg_pipeline.format('session')
    )
    # Desired spatial smoothing number
    parser = add_arg_if_in_arg_names('spat_smooth', arg_names, parser,  
        metavar='DESIRED_SPATIAL_SMOOTHING', 
        default=default_smooth, type=valid_whole_number,
        help=('Millimeters of spatial smoothing that you want for the surface '
              'and volume data.'
              + msg_whole_num + msg_default.format(default_smooth))
    )
    parser = add_arg_if_in_arg_names('subject', arg_names, parser, 
        metavar='SUBJECT_ID', required=True,
        type=lambda x: valid_subj_ses(x, 'sub-', 'subject'), #, 'NDAR', 'INV'),
        help='ID of subject to process.'
    )
    # Surface smoothing number
    parser = add_arg_if_in_arg_names('surf_smooth', arg_names, parser,  
        metavar='CURRENT_SURFACE_SMOOTHING', 
        default=default_smooth, type=valid_whole_number,
        help=''.join((msg_smooth.format('surface'), msg_whole_num,
                      msg_default.format(default_smooth)))
    )
    # Set file path for base directory and BIDS directory
    parser = add_arg_if_in_arg_names('study_dir', arg_names, parser,  
        metavar='BIDS_BASE_STUDY_DIRECTORY',
        type=valid_readable_dir, required=True, 
        help='Valid path to existing base study directory.'
    )
    parser = add_arg_if_in_arg_names('study_name', arg_names, parser, 
        metavar='STUDY_NAME', default=default_study_name,
        help=msg_pipeline.format('study')
    )
    # Which task you are running the pipeline on
    parser = add_arg_if_in_arg_names('task', arg_names, parser,  
        metavar='TASK_NAME', required=True,
        help=msg_pipeline.format('task')  # + msg_choices(choices_tasks)
    )
    parser = add_arg_if_in_arg_names('temp_dir', arg_names, parser, 
        type=valid_readable_dir, metavar='TEMPORARY_DIRECTORY',
        help=('Valid path to existing directory to save temporary files into.')
    )
    parser = add_arg_if_in_arg_names('templates', arg_names, parser, 
        type=valid_readable_dir, 
        help='Valid path to existing directory with template .fsf files.'
    )
    for lvl in default_runs_lvls:  # Specify the .fsf template files' names
        parser = add_arg_if_in_arg_names(
            'template{}'.format(lvl), arg_names, parser, 
            metavar='LEVEL_{}_TEMPLATE_NAME'.format(lvl),
            type=valid_template_filename, help=msg_template.format(lvl)
        )
    # Volume smoothing number
    parser = add_arg_if_in_arg_names('vol_smooth', arg_names, parser,  
        metavar='CURRENT_VOLUME_SMOOTHING',
        default=default_smooth, type=valid_whole_number,
        help=''.join((msg_smooth.format('volume'), msg_whole_num,
                      msg_default.format(default_smooth)))
    )
    # Specify path to wb_command
    parser = add_arg_if_in_arg_names('wb_command', arg_names, parser,  
        default=default_wb_command, type=valid_readable_file,
        help=('Path to wb_command file to run Workbench Command. If this flag '
              'is excluded, then the script will try to guess the path to '
              'the wb_command file by checking the user\'s BASH aliases. '
              'Your default wb_command is "{}". If '
              'that says "None", then you need to include this argument.'
              .format(default_wb_command))
    )
    # Argument used to get this script's dir
    parser = add_arg_if_in_arg_names(WRAPPER_LOC, arg_names, parser,
        type=valid_readable_dir, # required=True,
        help=('Valid path to existing ABCD-BIDS-task-fmri-pipeline directory '
              'that contains pipeline_wrapper.py')
    ) 
    return parser


def get_region_path_vars(cli_args, paths, run):
    """
    Build and return paths to particular brain region images' files/dirs
    by filling in the unique parts of generic path strings
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param run: Whole number (as an int or a string) defining which run this is
    :return: Tuple of string generic paths: design, func_str, subcort, surf_str
    """
    # Paths to design file base and subcortical volume stats directory
    design = os.path.join(paths['lvl_1']['fsf'],
                          get_sub_base(cli_args, run) + '_level1')
    subcort = os.path.join(paths['lvl_1']['parent'], 'SubcorticalVolumeStats')

    # Generic strings used as templates for paths later       
    func_str = os.path.join(paths['lvl_1']['intermediate'], 
                            '{}{}_filtered.atlasroi{}.{}.32k_fs_LR.func.gii')
    surf_str = os.path.join(paths['sub_ses']['anat'], (
        '{}_hemi-{}_space-MNI_mesh-fsLR32k_midthickness.surf.gii'
    ))
    return design, func_str, subcort, surf_str


def get_replacements(cli_args, **kwargs):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: Dictionary mapping variables' generic names in template files to 
             those variables' actual values provided by the user
    """
    replacements = {'SUBID': cli_args['subject'],
                    'FEAT_NAME': cli_args['study_name'], # Not paths['feat_name']
                    'FIN_SMOOTH': str(cli_args['spat_smooth']),
                    'HP_FILTER': str(cli_args['filter']), 
                    'SESSION': cli_args['ses'], 'TASK': cli_args['task'],
                    'OUTPUT_DIR': cli_args['output'],
                    'EVENTS_DIR': cli_args['events_dir'],
                    'STUDY_DIR': cli_args['study_dir']}
    replacements.update(kwargs)
    return replacements


def get_sbatch_args(cli_args, job):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param job: String 1-8 characters long naming the SBATCH job
    :return: List of strings, SLURM-related arguments to pass to the main
             script or level 1 analysis script for parallelization
    """
    return [argify('time', cli_args['time']), '-c', str(cli_args['cpus']),
            '-J', job, argify('mem', '{}gb'.format(cli_args["memory"]))]


def get_sub_base(cli_args, run_num=None):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param run_num: Whole number as an int or string defining which run this is
    :return: String identifying a subject, session, task, and maybe run
    """
    parts = [get_subj_ses(cli_args), 'task-' + cli_args['task']]
    if run_num is not None:
        parts.append('run-{}'.format(run_num))
    return '_'.join(parts)


def get_subj_ses(cli_args):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: String which combines --subject and --ses from command line
    """
    return '_'.join((cli_args['subject'], cli_args['ses']))


def get_TR_and_ntpts(dtseries_path, wb_command_path):
    """
    :param dtseries_path: String, the full path to a .dtseries.nii file
    :param wb_command_path: String, the full path to the wb_command executable
    :return: Tuple of 2 numbers, the number of timepoints and repetition time
    """
    if not os.path.exists(dtseries_path):
        sys.exit('Error: {} does not exist'.format(dtseries_path))
    else:
        ntpts = wb_command_get_info(wb_command_path, dtseries_path,
                                    'number-of-maps')
        rep_time = wb_command_get_info(wb_command_path, dtseries_path,
                                       'step-interval')
    return rep_time, ntpts


def glob_and_copy(dest_dirpath, *path_parts_to_glob):
    """
    Collect all files matching a glob string, then copy those files
    :param dest_dirpath: String, a valid path of a directory to copy files into
    :param path_parts_to_glob: Unpacked list of strings which join to form a
                               glob string of a path to copy files from
    """
    for file_src in glob(os.path.join(*path_parts_to_glob)):
        shutil.copy(file_src, dest_dirpath)


def make_and_save_confound_matrix(cli_args, desc_tsv_file, lvl_paths,
                                  sub_run_basename):
    """
    Create the confound matrix and copy it to subjects fsf_paths for each run
    :param cli_args: Dictionary containing all command-line arguments from user
    :param desc_tsv_file: String naming a .tsv file in intermediate_files/ dir
    :param lvl_paths: Dictionary mapping keys to dir path strings
    :param sub_run_basename: String naming the subject and the run number (?)
    :return: String, the base name of the confound matrix .csv file
    """
    # Local variables: File paths, step filename, adjusted variable to censor
    # initial frames based on user-specification, and result (confounds fname)
    in_file = os.path.join(lvl_paths['intermediate'], desc_tsv_file)
    def tsv_file_for_step(stepnum):
        return os.path.join(lvl_paths['intermediate'],
                            ('{0}_desc-filteredincludingFD_motion_step{1}.tsv'
                            .format(sub_run_basename, stepnum)))
    censor_volumes = list(range(0, cli_args['censor']))
    confounds_name = str(sub_run_basename + '_confound_matrix.tsv')
    
    # Read and write framewise displacement step1 .csv file
    df = pd.read_csv(in_file, sep='\s+')
    df.framewise_displacement.iloc[[censor_volumes]] = 1
    df.framewise_displacement[df.framewise_displacement < cli_args['fd']] = 0
    df.framewise_displacement[df.framewise_displacement > 0] = 1
    df.framewise_displacement.to_csv(tsv_file_for_step(1), header=False,
                                     encoding='utf-8', sep='\t', index=False)
    
    # Read and write step2 .csv file
    df = pd.read_csv(in_file, sep='\s+')
    cols = ['trans_x_mm', 'trans_y_mm', 'trans_z_mm', 'rot_x_degrees',
            'rot_y_degrees', 'rot_z_degrees', 'trans_x_mm_dt',
            'trans_y_mm_dt', 'trans_z_mm_dt', 'rot_x_degrees_dt',
            'rot_y_degrees_dt', 'rot_z_degrees_dt']
    df = df[cols]  # the 'cols' intermediate variable is needed to avoid error
    df.to_csv(tsv_file_for_step(2), sep='\t', encoding='utf-8',
              index=False, header=False)

    # Read and write step3 .csv file
    df = pd.read_csv(tsv_file_for_step(1), names=['A'], sep='\t')
    df = pd.concat([pd.get_dummies(df[df['A'] == 1].index)
                    .transpose(), df], axis=1).fillna(0)
    del df['A']
    df.to_csv(tsv_file_for_step(3), sep='\t', encoding='utf-8',
              index=False, header=False)
    
    # Read and write confound matrix .csv file; return its name
    pd.concat([pd.read_csv(tsv_file_for_step(x), header=None, sep='\t')
               for x in (2, 3)], axis=1).to_csv(
        os.path.join(lvl_paths['fsf'], confounds_name),
        sep='\t', encoding='utf-8', header=None, index=False
    )
    return confounds_name


def individualize_subprocess_run(run_args, run, to_replace):
    """
    Cycle through every argument in run_args and replace instances of the
    to_replace string with run, then return the arguments.
    :param run_args: List of strings, all arguments to call via subprocess
    :param run: Whole number (as an int or a string) defining which run this is
    :param to_replace: String to find and replace with each run name/id
    :return: run_args, but with to_replace replaced by run in them all
    """
    for i in range(len(run_args)):
        run_args[i] = str(run_args[i]).replace(to_replace, str(run))
    return run_args


def make_fake_nifti(cli_args, generic, old_smoothed, unique_part, cmd, *args):
    """
    Create a fake nifti from the smoothed dtseries for high-pass filtering 
    :param cli_args: Dictionary containing all command-line arguments from user
    :param generic: String, new smoothed nifti file path but with a '{}' in it
    :param old_smoothed: String, the path to a real old smoothed nifti file
    :param unique_part: String/etc inserted into generic to make a valid path
    :param cmd: String which is a Bash command but with '{}'s in it to replace
    :return: String, the valid path to the now-real new smoothed nifti file
    """
    started = datetime.now()
    new_smoothed = generic.format(unique_part)
    cmd_args = cmd.format(old_smoothed, *args, new_smoothed).split()
    if cmd_args[0] == 'wb_command':
        wb_command(cli_args, *cmd_args[1:])
    else:  # if cmd_args[0] in ('fsl', 'feat_model', 'film_gls', ):
        run_fsl_tool(cli_args, *cmd_args)
    if cli_args['print_progress']:
        get_and_print_time_since('started making '
                                 + os.path.basename(new_smoothed), started)
    return new_smoothed


def merge_files_in_range(cli_args, file_names, range_to, args):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param file_names: List of strings where each is a filename
    :param range_to: Integer, the number of files to merge
    :param args: List, the rest of the arguments to call merge_to_make_dtseries
    """
    for r in range(0, range_to):
        for f in file_names:
            merge_to_make_dtseries(cli_args, str(f) + str(r + 1), *args)


def merge_to_make_dtseries(cli_args, fname, lvl_paths, substats, AROI2, shape):
    """
    :param fname: String, base name of the files to merge into a dtseries
    :param lvl_paths: Dictionary mapping keys to dir path strings
    :param substats: String, the path to the subcortical stats directory
    :param AROI2: String, path to Atlas ROI file
    :param shape: Function takes 'L' or 'R' & returns path to shape.gii file
    """
    cii_out = os.path.join(lvl_paths['GOStats'], fname + '.dtseries.nii')
    subcort_in = os.path.join(substats, fname + '.nii.gz')
    func = lambda x: os.path.join(lvl_paths['parent'], x + '_SurfaceStats',
                                  fname + '.func.gii')
    fake_nifti = os.path.join(lvl_paths['GOStats'], fname + '.nii.gz')
    wb_command(cli_args, '-cifti-create-dense-timeseries', cii_out, '-volume',
               subcort_in, AROI2, '-left-metric', func('L'), '-roi-left', 
               shape('L'), '-right-metric', func('R'), '-roi-right', shape('R'))
    wb_command(cli_args, '-cifti-convert', '-to-nifti', cii_out, fake_nifti)


def organize_lvl_paths(lvl_paths, *keys_to_remove):
    """
    :param lvl_paths: Dictionary mapping keys to dir path strings
    :param keys_to_remove: Unpacked list of strings which are lvl_paths keys
                           to exclude from the return list
    :return: List of all values in lvl_paths (except the ones mapped to
             keys_to_remove), sorted alphabetically
    """
    lvl_paths = lvl_paths.copy()
    for each_key in keys_to_remove:
        lvl_paths.pop(each_key)
    to_return = list(lvl_paths.values())
    to_return.sort(reverse=False)
    return to_return


def overwrite_dirs(dirs_to_overwrite, mkdir=False):
    """
    :param dirs_to_overwrite: List of strings which are paths to directories
                              to create or overwrite with empty directories
    :param mkdir: True to remake all the dirs after overwrite, else False
    """
    for each_dir in dirs_to_overwrite:
        if os.path.isdir(each_dir):
            shutil.rmtree(each_dir)
        elif os.path.exists(each_dir):
            os.remove(each_dir)
        if mkdir:
            os.makedirs(each_dir)


def rand_string(L):
    """
    :param L: Integer, length of the string to randomly generate
    :return: String (of the given length L) of random characters
    """
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=L))
 

def rename_template_file_vars(old_template, new_template, replacements):
    """
    :param old_template: String, path to existing template file
    :param new_template: String, path to new template file which will be 
                         written with old_template variables but renamed
    :param replacements: Dictionary mapping each string in old_template to the
                         string to replace it with in new_template
    """
    with open(old_template) as infile:  # Open the level 1 or 2 template

        # Create new .fsf file; name the output "sub-*_ses-*_task-*_level*.fsf"
        with open(new_template, 'w') as outfile:
            for line in infile:

                # Use replacements dict to replace variables in the .fsf file
                for src, target in replacements.items():
                    line = line.replace(src, target)

                # Output the new subject-, (run-,) and task-specific .fsf file
                outfile.write(line)


def run_fsl_tool(cli_args, toolname, *args):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param toolname: String naming the executable tool in --fsl-dir to run
    :param args: Unpacked list of arguments to run toolname with
    """
    subprocess.check_call([
        valid_readable_file(os.path.join(cli_args['fsl_dir'], toolname)), *args
    ])


def run_parallel_or_sequential(script_path, cli_args, runs, to_replace, 
                               extra_args, second_fn=None, second_args=None):
    """
    Run a Python script via subprocess, either sequentially or in parallel
    depending on cli_args --no-parallel
    :param script_path: String, valid path to real script to run in parallel
    :param cli_args: Dictionary containing all command-line arguments from user
    :param runs: List of unique strings identifying differences between scripts
    :param to_replace: String to find and replace with each job name/id
    """
    if cli_args['no_parallel']: # Run processes serially/sequentially
        for run in runs:
            run_python_subscript(script_path, run, to_replace, *extra_args)

    else:  # Run processes in parallel using Python multiprocessing module
        to_run = list()
        all_args = list()
        for run in runs:
            all_args.append([script_path, run, to_replace, *extra_args])
            to_run.append(mp.Process(args=all_args[-1], name=all_args[-1][0],
                                     target=run_python_subscript))
        if second_fn and second_args:
            all_args.append(second_args)
            to_run.append(mp.Process(target=second_fn, args=second_args, 
                                     name=second_args[0]))
        if dict_has(cli_args, 'print_progress'):
            print('Running parallel:\n' + '\n\n'.join(str(x) for x in all_args))
        try:
            run_parallel(os.path.basename(script_path), to_run,
                        cli_args['sleep'], cli_args['print_progress'])
        except Exception as e:
            sys.exit(e)


def run_parallel(scriptname, processes, sleep, show):
    """
    Run a script multiple times in parallel
    :param scriptname: String describing the script being run in parallel
    :param processes: List of multiprocessing.Process objects ready to run
    :param sleep_secs: Integer, how many seconds to wait between (a) process
                       submissions and (b) checking if all processes finished
    :param show: True to show the user what's running at sleep_secs intervals;
                 otherwise False
    """
    started = datetime.now()
    submitted = list()
    failed = False
    for each_process in processes:
        submitted.append(each_process.start())
        time.sleep(sleep)
    while any((p.exitcode is None) for p in processes):
        time.sleep(sleep)
        if show:
            get_and_print_time_since(scriptname + ' started', started)
        if not all(p.exitcode is None or p.exitcode == 0 for p in processes):
            failed = True
            for p in processes:
                p.terminate()
        if failed:
            sys.exit('Error: {} subprocess failed.'.format(scriptname))


def run_python_subscript(path_to_subscript, run, to_replace, *args):
    """
    Use subprocess to run a Python 3.6+ script from this code base
    :param path_to_subscript: String, valid path to real Python 3.6+ script
    :param cli_args: Dictionary containing all command-line arguments from user
    :param run: Whole number (as an int or a string) defining which run this is
    :param to_replace: String to find and replace with each run name/id
    :param args: Unpacked list of parameters to run subscript with
    """
    start_time = datetime.now()
    try:
        subprocess.check_call(individualize_subprocess_run(
            ['python3', path_to_subscript, *args], run, to_replace
        )) 
    except subprocess.CalledProcessError:
        err_type, err_msg, _ = sys.exc_info()  # TODO make this into a reusable function? See run_level_1_analysis.get_events_make_template
        sys.exit('\n\n{}: {}\n\n'.format(err_type.__name__, err_msg))
    get_and_print_time_since(os.path.basename(path_to_subscript)
                             + ' started', start_time)
    return  # Explicitly end this function so multiprocessing knows it's done


def save_to_json_and_get_path(a_dict, dict_name, out_dir):
    """
    :param a_dict: Dictionary with only string keys
    :param dict_name: String naming a_dict
    :param out_dir: String, a valid path to a real directory to save
                    the .json file containing a_dict into
    :return: String, the full path to the .json file containing a_dict
    """
    json_path = os.path.join(out_dir, 'abcd-bids-pipeline-{}_{}.json'.format(
        dict_name, datetime.now().strftime('%Y-%b-%d_%H-%M')
    ))
    with open(json_path, 'w+') as json_file:
        json_file.write(json.dumps(a_dict))
    return json_path


def valid_float_0_to_1(val):
    """
    :param val: Object to check, then throw an error if it is invalid
    :return: val if it is a float between 0 and 1 (otherwise invalid)
    """
    return validate(val, lambda x: 0 <= float(x) <= 1, float,
                    'Value must be a number between 0 and 1')


def valid_output_dir(path):
    """
    Try to make a folder for new files at path; throw exception if that fails
    :param path: String which is a valid (not necessarily real) folder path
    :return: String which is a validated absolute path to real writeable folder
    """
    return validate(path, lambda x: os.access(x, os.W_OK),
                    valid_readable_dir, 'Cannot create directory at {}', 
                    lambda y: os.makedirs(y, exist_ok=True))


def valid_readable_dir(path):
    """
    :param path: Parameter to check if it represents a valid directory path
    :return: String representing a valid directory path
    """
    return validate(path, os.path.isdir, valid_readable_file,
                    'Cannot read directory at {}')


def valid_readable_file(path):
    """
    Throw exception unless parameter is a valid readable filepath string. Use
    this, not argparse.FileType('r') which leaves an open file handle.
    :param path: Parameter to check if it represents a valid filepath
    :return: String representing a valid filepath
    """
    return validate(path, lambda x: os.access(x, os.R_OK),
                    os.path.abspath, 'Cannot read file at {}')


def valid_readable_json(path):
    """
    :param path: Parameter to check if it represents a valid .json file path
    :return: String representing a valid .json file path
    """
    return validate(path, lambda x: os.path.splitext(path)[-1] == '.json',
                    valid_readable_file, '{} is not a readable .json filepath')


def valid_subj_ses(in_arg, prefix, name): #, *keywords):
    """
    :param in_arg: Object to check if it is a valid subject ID or session name
    :param prefix: String, 'sub-' or 'ses-'
    :param name: String describing what in_arg should be (e.g. 'subject')
    :return: True if in_arg is a valid subject ID or session name; else False
    """
    return validate(in_arg, lambda _: True, # lambda x: any([key in x for key in [prefix, *keywords]]),
                    lambda y: (y if y[:len(prefix)] == prefix else prefix + y),
                    '{}' + ' is not a valid {}'.format(name))


def valid_template_filename(fname):
    """
    :param fname: Parameter to check if it represents a .fsf file name
    :return: String representing the .fsf file name
    """
    return validate(fname, lambda x: os.path.splitext(x)[-1] == '.fsf',
                    lambda y: y, '{} is not an .fsf file name')


def valid_time_str(in_arg):
    """
    :param in_arg: Object to check if it's a time string in the HH:MM:SS format
    :return: True if in_arg is a time limit string in that format; else False
    """
    try:
        split = in_arg.split(":")
        assert len(split) == 3
        for each_num in split:
            assert each_num.isdigit()
            assert int(each_num) >= 0
        return in_arg
    except (TypeError, AssertionError, ValueError):
        raise argparse.ArgumentTypeError('Invalid time string.')


def valid_whole_number(to_validate):
    """
    Throw argparse exception unless to_validate is a positive integer
    :param to_validate: Object to test whether it is a positive integer
    :return: to_validate if it is a positive integer
    """
    return validate(to_validate, lambda x: int(x) >= 0, int, 
                    '{} is not a positive integer')


def validate(to_validate, is_real, make_valid, err_msg, prepare=None):
    """
    Parent/base function used by different type validation functions. Raises an
    argparse.ArgumentTypeError if the input object is somehow invalid.
    :param to_validate: String to check if it represents a valid object 
    :param is_real: Function which returns true iff to_validate is real
    :param make_valid: Function which returns a fully validated object
    :param err_msg: String to show to user to tell them what is invalid
    :param prepare: Function to run before validation
    :return: to_validate, but fully validated
    """
    try:
        if prepare:
            prepare(to_validate)
        assert is_real(to_validate)
        return make_valid(to_validate)
    except (OSError, TypeError, AssertionError, ValueError, 
            argparse.ArgumentTypeError):
        raise argparse.ArgumentTypeError(err_msg.format(to_validate))


def validate_cli_args(cli_args, parser, arg_names=set()):
    """
    Validate types and set defaults for any arg whose validation depends on
    another arg and therefore was not possible in get_pipeline_cli_argparser
    :param cli_args: Dictionary containing all command-line arguments from user
    :param parser: argparse.ArgumentParser to raise error if anything's invalid
    :param arg_names: Set containing SCAN_ARG if that argument is needed
    :return: cli_args, but fully validated
    """
    # Default levels, template file directory, and scanner info file path
    cli_args = ensure_dict_has(cli_args, 'levels', [1, 2]
                               if len(cli_args['runs']) > 1 else [1])
    cli_args = ensure_dict_has(cli_args, 'templates',
                               os.path.join(SCRIPT_DIR, 'templates'))
    if SCAN_ARG in arg_names:
        cli_args = ensure_dict_has(cli_args, SCAN_ARG, os.path.join(
            SCRIPT_DIR, 'scan_info', SCAN_ARG + '.csv'
        ))

    for lvl in cli_args['levels']:  # Default template file names
        cli_args = ensure_dict_has(cli_args, 'template{}'.format(lvl), (
            'template_DCAN_version_{}_level{}_UPDATED_FINAL.fsf'
            .format(cli_args['task'], lvl)
        ))  
        validate_template_file(cli_args, lvl, parser)
    
    # Default paths to FSL and wb_command
    ERR_MSG = 'No {} found. Please include the {} argument.'
    if not (dict_has(cli_args, 'wb_command') and
            os.access(cli_args['wb_command'], os.X_OK)):
        parser.error(ERR_MSG.format('wb_command executable', '--wb-command'))
    if not dict_has(cli_args, 'fsl_dir'):
        fsl = get_default_ext_command('fsl')
        cli_args['fsl_dir'] = os.path.dirname(fsl) if fsl else parser.error(
            ERR_MSG.format('FSL directory', '--fsl-dir')
        )

    # Default output/temp/event files directories. Avoiding ensure_dict_has to
    if not dict_has(cli_args, 'output'):       # prevent permissions error from
        cli_args['output'] = valid_output_dir( # valid_output_dir making dirs.
            os.path.join(cli_args['study_dir'], 'derivatives', 'abcd-bids-tfm'
                         'ri-pipeline', cli_args['subject'], cli_args['ses'])
        )
    for arg in ('temp_dir', 'events_dir'):
        if not dict_has(cli_args, arg):
            cli_args[arg] = valid_output_dir(
                os.path.join(cli_args['output'], 'level-1', arg.split('_')[0])
            )
    return cli_args


def validate_template_file(cli_args, lvl, parser):
    """
    Verify that template .fsf file exists
    :param cli_args: Dictionary containing all command-line arguments from user
    :param lvl: String or int defining the analysis level, 1 or 2 or "1" or "2"
    :param parser: argparse.ArgumentParser to raise error if anything's invalid 
    """
    tmpl = 'template{}'.format(lvl)
    tmpl_fpath = os.path.join(cli_args['templates'], cli_args[tmpl])
    if not os.path.exists(tmpl_fpath):
        parser.error('{} does not exist. Please re-run with a different --{} '
                     'or --templates argument.'.format(tmpl_fpath, tmpl))


def wb_command(cli_args, *args):
    """
    Call wb_command executable with any given parameters
    :param cli_args: Dictionary mapping 'wb_command' key to wb_command filepath
    :param args: List of all parameters to call wb_command with, in order
    """
    subprocess.check_call([cli_args['wb_command'], *args])


def wb_command_get_info(wb_command, dtseries, arg_only):
    """
    Call wb_command with -file-information and -no-map-info parameters
    :param wb_command: String, path to existing workbench wb_command executable
    :param dtseries: String, the path to a .dtseries.nii file with file info
    :param arg_only: String, the last part of the name of a wb_command
                     argument starting with '-only-'
    :return: String representing a numerical value describing the dtseries
    """
    return os.popen('{} -file-information {} -no-map-info -only-{}'
                    .format(wb_command, dtseries, arg_only)).read().rstrip()


def wb_LR_pair(func_LR, arg_LR=None, after=True):
    """
    Get wb_command left- and right- arguments
    :param func_LR: Function which accepts 'L' or 'R' and returns a filepath
    :param arg_LR: String naming the left- or right- argument
    :param after: True if arg_LR goes after the word 'left'/'right'; else False
    :return: List with 4 elements, arg name and then value for left then right
    """
    arg_LR = '-' + arg_LR if arg_LR else ''
    arg_fmt = '-{}' + arg_LR if after else arg_LR + '-{}'
    return [arg_fmt.format('left'), func_LR('L'),
            arg_fmt.format('right'), func_LR('R')]
