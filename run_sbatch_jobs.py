#!/usr/bin/env python3
# coding: utf-8

"""
Run many ABCD Task Prep Wrapper jobs in parallel via SBATCH
Greg Conan: gconan@umn.edu
Created: 2021-03-05
Updated: 2021-11-29
"""

# Import standard libraries
import argparse
from datetime import datetime
import os
import subprocess
import sys
import time


def find_local_imports(flg):
    """
    Ensure that this script can find its dependencies if parallel processing
    :param flg: String in sys.argv immediately preceding the path to return
    :return: String, path to the directory this Python script is in
    """
    try:
        parallel_flag_pos = -2 if sys.argv[-2] == flg else sys.argv.index(flg)
        sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))
    except (IndexError, ValueError):
        sys.path.append(os.path.dirname(os.path.abspath(sys.argv[0])))
    return sys.path[-1]


# Local custom imports
WRAPPER_LOC = '--wrapper-location'
find_local_imports(WRAPPER_LOC)
try:
    from src.pipeline_utilities import (
        add_slurm_args_to, argify, count_lines_in_txt_file, ensure_dict_has,
        get_sbatch_args, valid_output_dir, valid_readable_dir,
        valid_readable_file, valid_whole_number
    )
except (IndexError, ValueError):
    sys.exit('Please include the {} flag followed by the path to the '
             'directory containing pipeline_wrapper.py'.format(WRAPPER_LOC))


def main():
    # Start timing pipeline, get command line arguments from user, and move to
    # the folder where the user said to put the slurm out files
    started = datetime.now()
    cli_args = _cli()
    os.chdir(cli_args['slurm_out_dir'])

    # Submit every command listed in the --jobs-file as an SBATCH job
    all_commands = get_all_commands(cli_args)
    for command in all_commands:
        if cli_args['print_progress']:
            print(command)
        submit_job(lambda: subprocess.Popen(command), cli_args['sleep'],
                   cli_args['keep_trying'])
        time.sleep(cli_args['sleep'])
        submit_job(lambda: print(
                        'Submitted 1 more job. Jobs running: {}'
                        .format(count_jobs_running(cli_args['jobs_name']))
                   ), cli_args['sleep'], cli_args['keep_trying'])

    print('{} took this long to submit all {} pipeline jobs: {}'
          .format(os.path.basename(__file__), len(all_commands), datetime.now() - started))


def _cli():
    """
    :return: Dictionary containing all command-line arguments from user
    """
    parser = add_slurm_args_to(argparse.ArgumentParser())
    parser.add_argument(
        '--jobs-file', type=valid_readable_file, required=True,
        help=('Path to an existing readable file that has a list of exactly 1 '
              'bash command per line.')
    )
    parser.add_argument(
        '--jobs-name', type=valid_jobname, default='ABCDsbch',
        help='8-character string naming the SBATCH jobs.'
    )
    parser.add_argument(
        '--keep-trying', action='store_true',
        help=('Include this flag to re-submit any failed job. This flag will '
              'ensure that no job is skipped, but risks (1) stalling the '
              'submitter at 1 failed job or (2) thrashing the SLURM cluster.')
    )
    parser.add_argument('--partition', '-partition')
    parser.add_argument(
        '--slurm-out-dir', type=valid_output_dir, default='.',
        help='Path to a temp directory to save slurm out files in.'
    )
    parser.add_argument(
        '--start-ix', '-start', type=valid_whole_number, default=0,
        help='1 positive integer, the first line of --jobs-file to run.'
    )
    parser.add_argument(
        '--stop-ix', '-stop', type=valid_whole_number, # default=100,
        help='1 positive integer, the last line of --jobs-file to run.'
    )
    parser.add_argument(WRAPPER_LOC, type=valid_readable_dir)
    return validate_cli(vars(parser.parse_args()))


def valid_jobname(name):
    """
    :param name: Object to check if it's a valid string naming a slurm job
    :return: name if name is a 2- to 8-character-long string, else an error
    """
    strname = str(name)
    return strname if 1 < len(strname) < 9 else argparse.ArgumentTypeError(
        '{} is not a valid SBATCH job name. Please enter a job name between '
        '2 and 8 characters long instead.'.format(name)
    )


def validate_cli(cli_args):
    """
    If user didn't give a --stop-ix, then set the --stop-ix to the number of
    lines in the --jobs-file
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: cli_args, but with its --stop-ix value validated
    """
    cli_args = ensure_dict_has(
        cli_args, 'stop_ix', count_lines_in_txt_file(cli_args['jobs_file'])
    )
    valid_whole_number(cli_args['stop_ix'])
    return cli_args


def get_all_commands(cli_args):
    """
    Collect all SBATCH commands in a list
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: List of lists; each sublist is an SBATCH command ready to run
    """
    with open(cli_args['jobs_file'], 'r') as infile:
        cmd_num = 0
        ixs = None
        commands = list()
        for command in infile:
            if cli_args['start_ix'] <= cmd_num <= cli_args['stop_ix']:
                split_cmd = command.split()
                if not ixs:  # Get command structure from starting line/command
                    ixs = get_script_sub_ses_and_task_from_cmd(split_cmd)
                commands.append(get_sbatch_command(split_cmd, ixs, cli_args))
            cmd_num += 1
    return commands


def get_script_sub_ses_and_task_from_cmd(cmd_parts):
    """
    :param cmd_parts: List of string parts of complete command calling pipeline
    :return: Dictionary mapping each of several pipeline .py script flags 
             to the index of its value in split_cmd
    """
    flags_to_find = ['-subject', '-ses', '-task']
    flags_found = dict()
    for cmd_ix in range(len(cmd_parts)):
        each_flag = cmd_parts[cmd_ix]
        if 'script' not in flags_found and each_flag[-2:] == '.py':
            flags_found['script'] = cmd_ix
        for each_flag in flags_to_find:
            if cmd_parts[cmd_ix][-len(each_flag):] == each_flag:
                flags_to_find.pop(flags_to_find.index(each_flag))
                flags_found[each_flag[1:]] = cmd_ix + 1
    return flags_found


def get_sbatch_command(split_cmd, ixs, cli_args):
    """
    :param split_cmd: List of string parts of complete command calling pipeline
    :param ixs: Dictionary mapping each of several pipeline .py script flags 
                to the index of its value in split_cmd
    :param cli_args: Dictionary containing all command-line arguments from user
    :return: List of strings; an SBATCH command ready to run
    """
    sub_base = '_'.join([split_cmd[ixs['subject']], split_cmd[ixs['ses']],
                        split_cmd[ixs['task']]])
    out = argify('output', os.path.join(cli_args["slurm_out_dir"],
                                        "slurm-out-{}_{}.txt"))
    return ['sbatch', out.format(
                datetime.now().strftime("%Y-%m-%d_%H-%M"), sub_base
            ), *get_sbatch_args(cli_args, cli_args['jobs_name']),
            '--partition=' + cli_args['partition'], *split_cmd]


def submit_job(submit_the_job, sleeptime, keep_trying):
    """
    :param submit_the_job: Function to submit a job
    :param sleeptime: Integer, the number of seconds to wai
    :param keep_trying: True to keep trying to submit a job if its initial
                        submission fails, or False to skip failed jobs
    """
    # Try to submit every job, and if any of them fails, then try it again
    if keep_trying:   # (if user explicitly said to)
        submitted = False
        while not submitted:
            try:
                submit_the_job()
                submitted = True
            except subprocess.CalledProcessError:
                time.sleep(sleeptime)

    else: # Otherwise, only try to submit each job once
        submit_the_job()


def count_jobs_running(jb):
    """
    :param jb: String naming currently-running parallel slurm jobs
    :return: Integer counting how many ASA batch jobs are running right now
    """
    return subprocess.check_output("squeue", universal_newlines=True).count(jb)


if __name__ == '__main__':
    main()