#!/usr/bin/env python3
# coding: utf-8

"""
ABCD-BIDS Task fMRI Pipeline
Original: template_subject_specific_pipeline_DCAN_Public_Users_V2.py
Original Author: Anthony Juliano, PsyD, NERVE Lab, acjulian@uvm.edu
Original Created: 2020-10-27
Wrapper Author: Greg Conan, DCAN Lab, gconan@umn.edu
Wrapper Created: 2020-12-22
Wrapper Updated: 2021-12-03
"""

# Import standard libraries
from datetime import datetime
from glob import glob
import math
import os
import pandas as pd
import random
import shutil
import subprocess
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
        sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    return sys.path[-1]


# Constants: Paths to this dir and level 1 analysis script
WRAPPER_LOC = '--wrapper-location'
SCRIPT_DIR = find_myself(WRAPPER_LOC)
RUN_LVL_1_PY_FILE = os.path.join(SCRIPT_DIR, 'src', 'run_level_1_analysis.py')
COPEFEAT = 'cope{}.feat'


# Custom local imports
from src.pipeline_utilities import (
    add_slurm_args_to, argify, exit_with_time_info, get_all_analysis_paths,
    get_LR_functions, get_optional_cli_args, get_pipeline_cli_argparser,
    get_replacements, get_sub_base, organize_lvl_paths,
    rename_template_file_vars, run_fsl_tool,
    run_parallel_or_sequential, save_to_json_and_get_path, 
    validate_cli_args, wb_command
)


def main():
    # Time how long the script takes and get command-line arguments from user 
    start = datetime.now()
    cli_args = get_cli_args()

    # Calculate smoothing parameters for the dtseries files and get paths dict
    sigma = {x: add_smooth_sigma_round(cli_args, '{}_smooth'.format(x))
             for x in ('vol', 'surf')}
    paths = get_all_analysis_paths(cli_args)
    
    # Save the paths and sigma values into a .json file for analyses to use
    temp_json = save_to_json_and_get_path({'paths': paths, 'sigma': sigma},
                                          'temp', cli_args['temp_dir'])

    # Verify paths to left/right surface/shape .gii files
    run_sanity_checks(get_LR_functions(cli_args, paths), paths['AROI2'], start)

    # Skip processing if the subject/session/task has already processed
    sub_base = get_sub_base(cli_args)
    glob_existing_outputs = os.path.join(
        cli_args['output'], 'level-2', '*_files', sub_base + '*'
    )
    if not cli_args.pop('overwrite') and glob(glob_existing_outputs):
        print('\nSkipped pipeline for {} because output files already exist at'
              ' {}\n\nRerun with --overwrite to make new output files.\n'
              .format(sub_base, glob_existing_outputs))

    else:  # Main processing
        if 1 in cli_args['levels']:  # Level 1 processing
            print('\nStarting Level 1 analyses.')
            to_replace = '{replace_this}'
            run_parallel_or_sequential(
                RUN_LVL_1_PY_FILE, cli_args, cli_args['runs'], to_replace, 
                ['--temp-json', temp_json, *get_optional_cli_args(cli_args),
                 '--run-number=' + to_replace]
            )

        if 2 in cli_args['levels']:  

            # Catch error of user trying to run level-2 without level-1 files
            if not os.path.isdir(paths['dir_lvl']['1']):
                print('Error: Level 1 files must exist at {} in order to run '
                      'level 2 processing.'.format(paths['dir_lvl']['1']))
                exit_with_time_info(start, exit_code=1)

            # Level 2 processing
            os.chdir(cli_args['study_dir'])
            print('\nContinuing onto Level 2 analyses.')
            for lv2_dir in organize_lvl_paths(paths['lvl_2']): 
                os.makedirs(lv2_dir, exist_ok=True)
            make_dofs_and_masks_for_flameo_lv2(cli_args, paths)
            finish_lvl_2_analysis(cli_args, paths)
    
        if not cli_args['keep_all']:
            os.remove(temp_json)  # Cleanup by deleting temporary .json file

    # Show user how long the pipeline took and end the pipeline here
    exit_with_time_info(start)


def get_cli_args(slurm=True):
    """
    :return: Dictionary containing all command-line arguments from user
    """
    # Only include SLURM-specific command-line arguments if running in parallel
    parser = get_pipeline_cli_argparser()
    if slurm:
        parser = add_slurm_args_to(parser)

    # Extra argument for this main wrapper script, not used by other scripts
    parser.add_argument(
        '--overwrite', '-overwrite', action='store_true',
        help=('By default, the script will end if files exist at the --output '
              'directory. Include this flag to overwrite them.')
    )
    cli_args = validate_cli_args(vars(parser.parse_args()), parser)
    return cli_args


def run_sanity_checks(lr_fns, AROI2, start_time):
    """
    :param lr_fns: Dictionary with 2 local functions to get surface
                   and shape paths for both hemispheres
    :param AROI2: String, the path to the Atlas ROI file
    :param start_time: datetime.datetime object of when the script started
    """
    sanity_check((lr_fns['surf']('L'), lr_fns['surf']('R')),
                 'Surface Shape GIFTI', start_time)
    sanity_check((AROI2, lr_fns['shape']('L'), lr_fns['shape']('R')),
                 'dtseries template', start_time)
    print('All sanity checks were passed for this subject.\nThe sanity checks '
          'took this long to complete: {}'.format(datetime.now() - start_time))


def sanity_check(files_to_check, files_name, start_time):
    """
    :param files_to_check: List of strings, paths to the files to verify
    :param files_name: String naming the kind of files to check for
    :param start_time: datetime.datetime object of when the script started
    """
    if (os.path.exists(eachfile) for eachfile in files_to_check):
        print(files_name + ': OK')
    else:
        print('{0}s: MISSING\nPipeline is unable to be run for this subject '
              'due to them missing one or more {0} files.\nThe sanity checks '
              'for this subject failed.'.format(files_name))
        exit_with_time_info(datetime.now() - start_time, exit_code=1)


def add_smooth_sigma_round(cli_args, cur_smooth):
    """
    Pull user-defined smoothing variables, calculate the sigma value for
    additional smoothing, and round the sigma value to 7 decimal places
    :param cli_args: Dictionary containing all command-line arguments from user
    :param cur_smooth: String, cli_args key mapped to current smoothing value
    :return: Float; add sigma (surface or volume) round
    """
    add_smooth = math.sqrt((cli_args['spat_smooth']**2)
                           - (cli_args[cur_smooth]**2))
    return round((add_smooth / (2 * (math.sqrt(2 * math.log(2))))), 7)


def make_dofs_and_masks_for_flameo_lv2(cli_args, paths):
    """
    Create dof and mask files for input to flameo for Level 2 analysis by
    copying needed files from lvl 1 GOStats dir, reading them in as dataframes,
    making dof masks for each run based on residuals and dof values, merging
    the dof masks for each run, and then making a mask of the dof nii files
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    """
    dof_masks = list()
    fslmath_path = get_fslmath_path_fn(paths['lvl_2'])
    for run in cli_args['runs']:
        for filename in glob(os.path.join(paths['dir_lvl']['1'],
                                          (get_sub_base(cli_args, run)
                                              + '_' + paths['feat_name']),
                                          'GrayordinatesStats', '*')):
            shutil.copy(filename, paths['lvl_2'][run])   
        dof = pd.read_csv(os.path.join(paths['lvl_2'][run], 'dof'),
                          names=['A']).loc[0, 'A']
        dof_masks.append(fslmath_path(run, 'dofmask'))
        run_fsl_tool(cli_args, 'fslmaths', fslmath_path(run, 'res4d'), '-Tstd',
                    '-bin', '-mul', str(dof), dof_masks[-1])
    dof_nii = fslmath_path('parent', 'dof')
    mask = fslmath_path('parent', 'mask')
    run_fsl_tool(cli_args, 'fslmerge', '-t', dof_nii, *dof_masks)
    run_fsl_tool(cli_args, 'fslmaths', dof_nii, '-Tmin', '-bin', mask)


def get_fslmath_path_fn(lvl_2_paths):
    """
    :param lvl_2_paths: Dictionary mapping strings to lvl 2 file path strings
    :return: Function which accepts the key to a lvl_2_paths pair and a
             filename string, then returns the path to a .nii.gz file.
    """
    return lambda key, fname: os.path.join(lvl_2_paths[key], fname + '.nii.gz')


def get_cope_filenames(lvl_2_paths, copefile):
    """
    :param lvl_2_paths: Dictionary with paths to feature directory and
                        cope file directories as its values
    :param copefile: String, the basename of the files to return paths to
    :return: Tuple of 3 strings, the paths to cope files 1 and 2 plus the path
             to the cope file in the feat directory
    """
    return [os.path.join(lvl_2_dir, copefile) for lvl_2_dir in
            organize_lvl_paths(lvl_2_paths, 'fsf')]


def finish_lvl_2_analysis(cli_args, paths):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    """
    # Merge cope & varcope files for each run's flameo (Level 2 analysis) call
    lv2_base = os.path.join(paths['lvl_2']['fsf'],
                            get_sub_base(cli_args) + '_level2')
    shutil.copy(paths['template2'], paths['lvl_2']['fsf'])
    rename_template_file_vars(
        os.path.join(paths['lvl_2']['fsf'], cli_args['template2']),
        os.path.join(paths['lvl_2']['fsf'], lv2_base + '.fsf'),
        get_replacements(cli_args)
    )
    run_fsl_tool(cli_args, 'feat_model', lv2_base)

    numbers_of_copes = list()
    for eachrun in cli_args['runs']:
        numbers_of_copes.append(len(glob(
            os.path.join(paths['lvl_2'][eachrun], 'cope*dtseries.nii')
        )) + 1)
        fslmerge_for_lvl_2(cli_args, numbers_of_copes[-1], paths['lvl_2'])
    number_of_copes = max(numbers_of_copes)

    # Run 2nd level analysis using flameo for each contrast
    for r in range(1, number_of_copes):
        flameo_lvl_2_analysis(
            cli_args, paths['lvl_2']['parent'], r, runmode='fe',
            tc=lv2_base + '.con', cope='cope{}.nii.gz', ld=COPEFEAT,
            dm=lv2_base + '.mat', vc='varcope{}.nii.gz', dvc='dof.nii.gz', 
            cs=lv2_base + '.grp', mask='mask.nii.gz'
        )

    # Convert Level 2 fakeNIFTI Files back to CIFTI
    ciitmpl = os.path.join(
        paths['dir_lvl']['1'], '_'.join((
            get_sub_base(cli_args, random.choice(cli_args['runs'])),
            paths['feat_name']
        )), 'GrayordinatesStats', 'pe1.dtseries.nii'
    )
    for r in range(1, number_of_copes):
        for fname in ('zstat1', 'zflame1uppertstat1', 'zflame1lowertstat1',
                      'weights1', 'varcope1', 'tstat1', 'tdof_t1', 'res4d',
                      'pe1', 'mean_random_effects_var1', 'mask', 'cope1'):
            convert_lv2_fake_nii_to_cii(cli_args, fname, r, paths, ciitmpl)

    # If user doesn't want to keep all files, then move/rename the necessary
    # files and delete the rest, including temp files from level 1 stats dirs
    only_keep_needed_files(cli_args, paths, number_of_copes)


def fslmerge_for_lvl_2(cli_args, number_of_copes, lv2_paths):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param number_of_copes: Integer, how many cope files are in each Level 2 dir
    :param lv2_paths: Dictionary with paths to feature directory and
                      cope file directories as its values    
    """
    for r in range(1, number_of_copes):
        for each_fname in ('cope{}.nii.gz', 'varcope{}.nii.gz'):
            run_fsl_tool(cli_args, 'fslmerge', '-t',
                         *get_cope_filenames(lv2_paths, each_fname.format(r)))


def flameo_lvl_2_analysis(cli_args, feat_dir, loop_n, runmode, **kwargs):
    """
    Call the FSL command 'flameo' with a certain set of parameters
    :param feat_dir: String, the valid path to a feature directory
    :param loop_n: Integer, the loop-number to add to some of the files' names
    :param runmode: String, the --runmode value for a flameo call
    :param kwargs: Unpacked dictionary mapping each flameo option name to a
                   file path which is either complete or has '{}' in r's place
    """
    args_list = ['flameo', argify('runmode', runmode)] 
    for argname, filename in kwargs.items():
        filename = filename.format(loop_n) if '{}' in filename else filename
        args_list.append(argify(argname, os.path.join(feat_dir, filename)))
    run_fsl_tool(cli_args, *args_list)


def convert_lv2_fake_nii_to_cii(cli_args, fname, loop_n, paths, cii_template):
    """
    Run wb_command -cifti-convert to turn level 2 'fake NIFTI' files into CIFTI
    .dtseries.nii files
    :param cli_args: Dictionary containing all command-line arguments from user
    :param fname: String naming the fake nifti file to convert to cifti
    :param loop_n: Integer, the loop-number to add to the file parent dir name
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param cii_template: String, path to CIFTI template file
    """
    lv2_feat_dir = paths['lvl_2']['parent']
    fake_nii = get_cope_fpath(lv2_feat_dir, loop_n, '{}.nii.gz'.format(fname))
    cii = get_cope_fpath(lv2_feat_dir, loop_n, '{}.dtseries.nii'.format(fname))
    wb_command(cli_args, '-cifti-convert', '-from-nifti', fake_nii,
               cii_template, cii, '-reset-timepoints', '1', '1')


def get_cope_fpath(feat_dir, loop_n, fname_w_ext):
    """
    :param feat_dir: String, the path to a feature directory
    :param loop_n: Integer, loop-number to add to cope file parent dir name
    :param fname_w_ext: String, basename and file extension of the cope file
    :return: String, the full valid path to a cope file
    """
    return os.path.join(feat_dir, COPEFEAT.format(loop_n), fname_w_ext)


def only_keep_needed_files(cli_args, paths, num_dtseries):
    """
    Move necessary output files into a separate directory from temp files
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param num_dtseries: Int, how many .dtseries.nii files (+1 for 0 indexing)
    """
    # Make a new directory structure to keep essential pipeline output files in
    essentials = {'cope': 'cope1', 'log': 'logfile', 'dof': 'tdof_t1', 
                  'mask': None, 'res4d': None, 'varcope': 'varcope1'}
    ess_dir_lv = {lvl: os.path.join(cli_args['output'], 'level-{}'.format(lvl))
                  for lvl in (1, max(cli_args['levels']))}
    for ess in essentials.keys():
        os.makedirs(os.path.join(ess_dir_lv[2], ess + '_files'), exist_ok=True)

    # Move essential pipeline output files to new directory structure
    old_fdir = os.path.join(paths['lvl_2']['parent'], COPEFEAT)
    new_file = os.path.join(ess_dir_lv[2], '{}_files',
                            get_sub_base(cli_args) + '_contrast_{}_')
    move_or_copy = shutil.copy2 if cli_args['keep_all'] else shutil.move
    for r in range(1, num_dtseries):
        for eachfname, renamed in essentials.items():
            fname = renamed if renamed else eachfname
            fname = fname + '.dtseries.nii' if fname != 'logfile' else fname
            new_location = new_file.format(eachfname, r) + (
                fname[:-14] + '.dtseries.nii' if fname[:4] == 'cope' else fname
            )

            try:  # Handle the error of a missing contrast output file
                move_or_copy(os.path.join(old_fdir.format(r), fname),
                             new_location)
            except FileNotFoundError:
                print(old_fdir.format(r) + ' not found')

    # Keep the level 1 .dtseries.nii GrayordinatesStats, log, and dof files
    if not cli_args['keep_all']:
        for run in cli_args['runs']:
            dst = os.path.join(ess_dir_lv[1],
                               os.path.basename(paths['lvl_2'][run]))
            os.makedirs(dst, exist_ok=True)
            for eachfile in os.scandir(paths['lvl_2'][run]):
                prefix = eachfile.name[:4]
                if eachfile.name[-6:] != 'nii.gz' and (
                    (prefix == 'cope' and not is_blank_cope(eachfile.path,
                                                            cli_args))
                    or prefix == 'res4' or '.' not in eachfile.name
                ): 
                    # Rename the files to specify by task (and subject/session)
                    os.rename(eachfile.path, os.path.join(
                        dst, '_'.join((get_sub_base(cli_args, run),
                                       eachfile.name))
                    ))

        # Delete all files for this sub/ses/task in both levels' feat dirs
        shutil.rmtree(paths['lvl_2']['parent']) # each_lvl_path)
        for each_lvl_path in paths['dir_lvl'].values():
            for feat_dir in glob(os.path.join(
                each_lvl_path, get_sub_base(cli_args) + "*feat"
            )):
                shutil.rmtree(feat_dir, ignore_errors=True)

            # Remove any empty feat directories
            try:
                os.rmdir(each_lvl_path)
            except OSError:
                pass


def is_blank_cope(cope_path, cli_args):
    """
    :param cope_file_path: String, valid path to existing cope file
    :param cli_args: Dictionary mapping 'wb_command' to a valid path
    :return: True if the cope file at cope_file_path is blank; else False
    """
    return subprocess.check_output((
        cli_args['wb_command'], '-cifti-stats', cope_path, '-reduce', 'STDEV'
    )).decode('utf-8').rstrip() == '0'


if __name__ == '__main__':
    main()