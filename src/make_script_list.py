#!/usr/bin/env python3
# coding: utf-8

"""
Create list of Bash commands to run ABCD task prep jobs simultaneously with SBATCH
Greg Conan: gconan@umn.edu
Created: 2021-02-26
Updated: 2021-09-20
"""

# Import standard libraries
import argparse
from glob import glob
import os
import sys


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
                        os.path.abspath(sys.argv[0]))))
    return sys.path[-1]


# Constants: Paths to this dir and level 1 analysis script
WRAPPER_LOC = '--wrapper-location'
SCRIPT_DIR = find_myself(WRAPPER_LOC)
from src.pipeline_utilities import (
    dict_has, ensure_dict_has, get_main_pipeline_arg_names,
    get_optional_cli_args, get_pipeline_cli_argparser, get_sub_base,
    valid_output_dir, valid_readable_dir, valid_readable_file,  
)


def main():
    # Get all command-line arguments from user
    cli_args = ensure_dict_has(get_cli_args(), 'script', os.path.join(
        SCRIPT_DIR, 'pipeline_wrapper.py'
    ))

    # If user didn't use --sourcedata flag, use --study-dir instead
    sourcedata_given = dict_has(cli_args, 'sourcedata')
    source = (cli_args['sourcedata'] if sourcedata_given else
              cli_args['study_dir'])

    # Get commands to run every subject, session, & task in the sourcedata dir
    # thisdir = os.path.abspath(os.path.dirname(__file__))
    all_commands = list()
    for task in cli_args['tasks']:
        for subject in os.scandir(source):
            if subject.is_dir() and 'sub-NDAR' in subject.name:
                for session in os.scandir(subject.path):
                    if (session.is_dir() and 'ses-' in session.name and
                        has_events(subject, session, task, cli_args,
                                   sourcedata_given)):
                            all_commands.append(get_cmd(
                                cli_args, ses=session.name, task=task,
                                subject=subject.name
                            ))

    # Write out all commands to output file
    with open(cli_args['script_list'], 'w+') as outfile:
        outfile.write('\n'.join(' '.join(cmd) for cmd in all_commands))


def has_events(sub, ses, task, cli_args, sourcedata_given):
    """
    :param sub: os.DirEntry representing a subject's sourcedata dir
    :param ses: os.DirEntry representing a subject's session's sourcedata dir
    :param task: String, the name of the task to process for a subject session
    :param cli_args: Dictionary containing all validated command-line arguments
    :param sourcedata_given: True if the --sourcedata flag was used; else False
    :return: True if sourcedata_given or sub has task's event files; else False
    """
    cli_args.update(subject=sub.name, ses=ses.name, task=task)
    return True if not sourcedata_given else len(glob(os.path.join(
        ses.path, 'func', get_sub_base(cli_args) + '*'
    ))) > 0


def get_cli_args():
    """
    :return: argparse.Namespace with command-line arguments from user
    """
    args = get_main_pipeline_arg_names().difference({
        'output', 'ses', 'subject', 'task', WRAPPER_LOC[2:].replace('-', '_')
    })
    tasks = ('SST', 'MID', 'nback')
    parser = get_pipeline_cli_argparser(arg_names=args)
    parser.add_argument('-all-events', '--all-events', type=valid_readable_dir,
                        help=('Valid path to an existing directory which has '
                              '1 folder per subject, with the folder structure '
                              '(--all-events)/(subject ID)/(session name)/'
                              'level-1/events.'))
    parser.add_argument('-all-outputs', '--all-outputs', type=valid_output_dir,
                        help=('Valid path to your output directory root. In '
                              'other words, the "--output" argument for each '
                              'command in the --script-list file will be '
                              'subject- and session-specific subdirectories '
                              'of this --all-outputs directory.'))
    parser.add_argument('-output', '--output', type=valid_output_dir, required=False)
    parser.add_argument('-script', '--script', type=valid_readable_file)
    parser.add_argument('-script-list', '--script-list', required=True)
    parser.add_argument('-slurm', '--slurm', action='store_true')
    parser.add_argument('-sourcedata', '--sourcedata', type=valid_readable_dir)
    parser.add_argument('-tasks', '--tasks', nargs='+', default=tasks) # choices=tasks, 
    parser.add_argument(WRAPPER_LOC, type=valid_readable_dir,
                        default=SCRIPT_DIR) #, dest='loc')
    return vars(parser.parse_args())


def get_cmd(cli_args, **run_specific_args): # thisdir, subject, session, task):
    """
    :param cli_args: Dictionary containing all validated command-line arguments
    :param thisdir: String, the path to this Python script's parent directory
    :param subject: os.DirEntry representing a subject
    :param session: os.DirEntry representing a subject's session
    :param task: String, the name of the task to process for a subject session
    :return: List of all parameters in a Bash command calling the full pipeline
             wrapper Python script for one task for one session for one subject
    """
    final_args = cli_args.copy()
    if dict_has(cli_args, 'all_events'):
        final_args['events_dir'] = os.path.join(
            cli_args['all_events'], run_specific_args['subject'],
            run_specific_args['ses'], 'level-1', 'events'
        )
    if dict_has(cli_args, 'all_outputs'):
        final_args['output'] = os.path.join(cli_args['all_outputs'],
                                            run_specific_args['subject'],
                                            run_specific_args['ses'])
    for excluded_arg in ('script', 'script_list', 'slurm', 'sourcedata', 'tasks',
                         'all_events', 'all_outputs'):
        final_args.pop(excluded_arg)
    for included_arg_name, included_arg_value in run_specific_args.items():
        final_args[included_arg_name] = included_arg_value
    return [cli_args['script'], *get_optional_cli_args(final_args)]
    
    '''
    pll = (['--parallel', '--time=1:30:00']
           if dict_has(cli_args, 'parallel') else [])
    loc = ['--wrapper-location', cli_args['wrapper_location']] # if cli_args['loc'] else []
    studydir = (['--study-dir', cli_args['study_dir']]
                if dict_has(cli_args, 'study_dir') else [])
    outdirs = [cli_args['output_dir'], subject.name, session.name]
    if dict_has(cli_args, 'level'):
        outdirs.append('level-{}'.format(cli_args['level']))
    return [cli_args['script'], *studydir, '-ses', session.name, *pll,
            '-task', task, '--subject', subject.name, *loc, '-out',
            os.path.join(*outdirs)]
    '''

if __name__ == '__main__':
    main()