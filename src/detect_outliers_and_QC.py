#!/usr/bin/env python3
# coding: utf-8

"""
Greg Conan: gconan@umn.edu
Created: 2021-04-16
Updated: 2021-10-01
"""

# Import standard libraries
from datetime import datetime
from glob import glob
import numpy as np
import os
import pandas as pd
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


# Constants
BLANK = '_blank_EVs'
COLS = [x + y for y in ('MIN', 'MAX', 'MEAN', 'STDEV')
        for x in ('vol_', 'L_', 'R_', 'dtseries_')]
READY1 = '_ready_for_stats'  # String in names of .tsv files made by step 3
READY2 = '_ready_for_analyses' # String in names of .tsv files made by step 5

# Local custom imports
WRAPPER_LOC = '--wrapper-location'
SCRIPT_DIR = find_myself(WRAPPER_LOC)
from src.pipeline_utilities import (
    dict_has, ensure_dict_has, exit_with_time_info, get_pipeline_cli_argparser,
    get_sub_base, valid_output_dir, valid_subj_ses, valid_whole_number
)


def main():
    cli_args = _cli()
    start_time = datetime.now()

    # If user didn't give contrast numbers, get defaults for each task
    cli_args = ensure_dict_has(cli_args, 'contrasts', [x for x in range(
        1, {'MID': 42, 'SST': 15, 'nback': 19}[cli_args['task']]
    )])

    # Run whichever step(s) the user specified, then print how long it took
    print('{} started at {}'.format(sys.argv[0], start_time))
    all_steps = {1: save_summary_stats_for_subject,
                 2: gather_organize_and_sort_EVs,
                 3: separate_all_blank_EVs,
                 4: screen_out_blank_dependent_contrasts,
                 5: screen_out_4_stdev_outliers,
                 6: make_subjects_copes_chart}
    for step in cli_args.pop('step'):
        all_steps[step](cli_args)
    exit_with_time_info(start_time, 0)


def _cli():
    """
    Get all command-line arguments from user
    - A lot of this is copied from pipeline_utilities.py instead of imported
      because some flags 'required' there are not 'required' here, due to
      variation between which parameters each 'step' requires
    :return: Dictionary containing all command-line arguments from user
    """
    # Default values for user input arguments
    default_steps = [x + 1 for x in range(6)]

    # Strings used in (usually multiple) help messages
    msg_default = ' By default, this argument\'s value(s) will be {}.'
    msg_choices = lambda x: (' This argument\'s value must be one of these: {}'
                             .format(x) + msg_default.format(x[0]))
    msg_steps = (  # Help message defining what each step does
        'Which step(s) of the QC and outlier-detection process to run.\n1. '
        'Calculate and save summary statistics for each contrast for every '
        'subject.\n2. Collect those summary statistics from each specific '
        'subject\'s .tsv files and save them all into one master .tsv file '
        'per contrast.\n3. Split the master .tsv file into 2 more: blank and '
        'non-blank.\n4. Split the non-blank .tsv file into 2 more: contrasts '
        'whose meaning depends on earlier (blank) contrasts, and contrasts '
        'whose does not.\5. Split the non-blank-dependent (or just non-blank) '
        '.tsv file into 2 more: 1 with all outliers (defined as 4 standard '
        'deviations from the mean) and 1 with the rest.\n6. Make a new .tsv '
        'chart of which subjects have valid copes for which contrast numbers.'
    ) + msg_choices(default_steps)

    # Create parser with command-line arguments from user
    parser = get_pipeline_cli_argparser(arg_names=(
        'contrasts', 'output', 'ses', 'study_dir', 'task', 'wb_command'
    ))
    parser.add_argument(
        '-step', '--step', type=valid_whole_number, choices=default_steps,
        default=default_steps, nargs='+', help=msg_steps
    )
    parser.add_argument(
        '-subject', '--subject', metavar='SUBJECT_ID',
        type=lambda x: valid_subj_ses(x, 'sub-', 'subject', 'NDAR', 'INV'),
        help='ID of subject to process. Not needed for most steps.'
    )
    parser.add_argument(WRAPPER_LOC)
    return vars(parser.parse_args())


def save_summary_stats_for_subject(cli_args):
    """
    Calculate summary statistics for each contrast for 1 subject, then save
    those summary stats into 3 new files (left hemi, right hemi, and subcort)
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    # Get file paths for directories containing subject data
    cli_args = ensure_dict_has(cli_args, 'output', (
        get_QC_dirpath(cli_args, subject=cli_args['subject'])
    ))
    sub_basename = get_sub_base(cli_args)
    outfilepath = get_generic_contrast_path(cli_args['output'], sub_basename)
    get_cope_path = get_cope_path_function(cli_args['output'], sub_basename)
    cope_check = get_cope_path(1, '.dtseries.n')

    if os.path.exists(cope_check):
        for contrast in cli_args['contrasts']:
            
            # Get all cope file paths
            dtseries = get_cope_path(contrast, '.dtseries.n')
            volume = get_cope_path(contrast, '.n')
            hemi_L = get_cope_path(contrast, '_L.func.g')
            hemi_R = get_cope_path(contrast, '_R.func.g')

            # Split each cope file into 3: volume, left hemi, and right hemi
            os.system('wb_command -cifti-separate-all {} -volume {} -left {} '
                      '-right {}'.format(dtseries, volume, hemi_L, hemi_R))

            # Calculate summary statistics for all copes and save them to .csv
            d = {'sub': sub_basename}
            for reduct in ('MIN', 'MAX', 'MEAN', 'STDEV'):
                d['vol_' + reduct] = wb_stats('volume', volume, reduct)
                d['L_' + reduct] = wb_stats('metric', hemi_L, reduct)
                d['R_' + reduct] = wb_stats('metric', hemi_R, reduct)
                d['dtseries_' + reduct] = wb_stats('cifti', dtseries, reduct)
            save_to_csv(d, outfilepath.format(contrast))

    else:  # If a cope does not exist, save out a blank stats file for it
        d = {'sub': sub_basename}
        for reduct in ('MIN', 'MAX', 'MEAN', 'STDEV'):
            d['vol_' + reduct] = 'n/a'
            d['L_' + reduct] = 'n/a'
            d['R_' + reduct] = 'n/a'
            d['dtseries_' + reduct] = 'n/a'
        for contrast in cli_args['contrasts']:   
            save_to_csv(d, outfilepath.format(contrast))


def get_cope_path_function(outdir, sub_base):
    """
    :param outdir: String, valid path to 'level-2' output directory
    :param sub_base: String identifying subject, session, and task
    :return: Function accepting a contrast number and a string defining a valid
             NIfTI file extension and returning a valid path to a contrast file
    """
    return lambda contrast, nii_ext: os.path.join(
        outdir, 'cope_files', sub_base + ('_contrast_{}_cope1{}ii'
                                          .format(contrast, nii_ext))
    )


def wb_stats(stat, filepath, reduction):
    """
    Call wb_command to get statistical information about a contrast file
    :param stat: String naming the wb_command statistic argument, namely
                 'volume', 'metric', or 'cifti'
    :param filepath: String, valid path to NIfTI file
    :param reduction: String defining the statistical output to return, namely
                      'MAX', 'MEAN, 'MIN', or 'STDEV'
    :return: String representing a number (the max/mean/min/stdev) for the file
    """
    return os.popen('wb_command -{}-stats {} -reduce {}'
                    .format(stat, filepath, reduction)).read().rstrip()


def save_to_csv(data, outpath):
    """
    :param data: Dictionary of lists, to be made into a pandas.DataFrame
    :param outpath: String, valid path to .csv file to save data into
    """
    df = pd.DataFrame(data=data, index=[0])
    df.to_csv(outpath, sep='\t', encoding='utf-8', header=None, index=False)


def gather_organize_and_sort_EVs(cli_args):
    """
    Given that each subject already has its own directory containing QC .tsv
    files, collect the data from all of those individual subject files and
    save it into one .tsv file per contrast containing all subjects' data
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    cli_args = ensure_dict_has(cli_args, 'output', get_QC_dirpath(cli_args))
    path_allsubjs = valid_output_dir(cli_args['output'])
    for contrast in cli_args['contrasts']:  # Prepare .tsv files
        newtsv = get_generic_contrast_path(path_allsubjs, cli_args['task']
                                           ).format(contrast)
        if not os.path.exists(newtsv):
            with open(newtsv, 'w+') as tsv_file:
                tsv_file.write('\t'.join(['sub_basename', *COLS]) + '\n')
    if dict_has(cli_args, 'subject'):
        sub_basename = get_sub_base(cli_args)
        path_per_subj = get_QC_dirpath(cli_args, cli_args['subject'])
        organize_EVs_for_subj(sub_basename, path_per_subj,
                              path_allsubjs, cli_args)
    else:
        for subj_dir in os.scandir(cli_args['study_dir']):
            path_per_subj = get_QC_dirpath(cli_args, subj_dir.name)
            sub_basename = get_sub_base({'subject': subj_dir.name,
                                         'ses': cli_args['ses'],
                                         'task': cli_args['task']})
            organize_EVs_for_subj(sub_basename, path_per_subj,
                                  path_allsubjs, cli_args)


def get_QC_dirpath(cli_args, subject=None):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param subject: String, ID of subject to get the QC .tsv file path for
    :return: String, the full path to a directory with QC .tsv files, either
             for all subjects (if subject=None) or otherwise for 1 subject
    """
    subj_dirs = ([cli_args['study_dir'], subject, cli_args['ses'], 'level-2']
                 if subject else [os.path.dirname(cli_args['study_dir'])])
    return os.path.join(*subj_dirs, 'QC_files', 'Intensity_TSVs')


def organize_EVs_for_subj(sub_basename, path_in, path_out, cli_args):
    """
    Get the QC data from one subject's (session's task's) .tsv files for each
    contrast, then append that data as a row in a master .tsv file containing
    all subjects' data for a given contrast number
    :param sub_basename: String naming the subject ID, session, and task
    :param path_in: Valid path to existing .tsv file to read QC data in from
    :param path_out: Valid path to directory with .tsv file to append QC data into
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    contrast_generic = get_generic_contrast_path(path_in, sub_basename)
    new_ix_col = 'subject'
    for r in cli_args['contrasts']:
        file_in = contrast_generic.format(r)
        if os.path.exists(file_in):
            outfile = get_generic_contrast_path(path_out, cli_args['task']
                                                ).format(r)
            df = pd.read_csv(file_in, sep='\t')
            df[sub_basename.split('_')[0]] = new_ix_col
            with open(outfile, 'a') as outfile:
                df.to_csv(outfile, sep='\t', header=True, index=False)


def get_generic_contrast_path(parent_dir, unique_identifier):
    """
    :param parent_dir: Valid path to an existing directory
    :param unique_identifier: String identifying 1 subject, task, and session
    :return: String, the full path to a .tsv file containing QC data for a
             contrast, but with the contrast number replaced by '{}' so the 
             resulting generic path can be used for any contrast 
    """
    return os.path.join(parent_dir, unique_identifier +
                        '_contrast_{}_Intensity_QC.tsv')


def separate_all_blank_EVs(cli_args):
    """
    Separate blank copes from non-blank ones for each contrast number
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    cli_args = ensure_dict_has(cli_args, 'output', get_QC_dirpath(cli_args))
    contrast_generic = get_generic_contrast_path(cli_args['output'],
                                                 cli_args['task'])
    for contrast in cli_args['contrasts']:
        df = pd.read_csv(contrast_generic.format(contrast), sep='\t')
        df.set_index('sub_basename', inplace=True, drop=True)
        df = df.apply(lambda x: (x.apply(lambda _: np.nan)
                                 if is_messed_up(x) else x), axis=1)
        separate_out_blank_EVs_from(
            df, get_generic_intensity_QC_TSV_path(cli_args, contrast) 
        )


def separate_out_blank_EVs_from(a_df, outfilepath_generic):
    """
    Split a .tsv file containing all subjects' data from 1 contrast into 2 new
    .tsv files containing blank and non-blank contrast files' data respectively
    :param a_df: pandas.DataFrame with a 'vol_STDEV' column containing mostly
                 numerical data and some NaNs
    :param outfilepath_generic: String, the path to a .tsv file with 1 unique
                                part in its name replaced by '{}'
    """
    df2 = a_df[pd.to_numeric(a_df['vol_STDEV'], errors='coerce').isnull()
               ].applymap(lambda _: np.nan)
    a_df = a_df[~a_df.isin(df2)].dropna().apply(floatify, axis='columns')
    to_tsv(a_df, outfilepath_generic.format(READY1))
    to_tsv(df2, outfilepath_generic.format(BLANK))


def to_tsv(a_df, tsv_path):
    """
    Save a pandas.DataFrame into a .tsv file at a given path
    :param a_df: pandas.DataFrame (any as long as it has a header and index)
    :param tsv_path: String, valid path to a .tsv file to write a_df into
    """
    a_df.to_csv(tsv_path, sep='\t', header=True, index=True)


def get_generic_intensity_QC_TSV_path(cli_args, contrast):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param contrast: Integer, the contrast number to get the path for
    :return: String, the path to a .tsv file for one contrast and one path,
             but with a unique part of that .tsv filename replaced by '{}'
    """
    return os.path.join(cli_args['output'], cli_args['task'] + (
        '_contrast_{}'.format(contrast) + '_Intensity_QC{}.tsv'
    ))


def is_messed_up(series, printit=False):
    """
    :param series: pandas.Series filled with either numerical data or strings
                   representing numerical data, some of which may have 1 too
                   many decimal points (let's call them "messed up strings")
    :param printit: True to print that series has messed up strings; else False
    :return: True if the series has messed up strings; else False
    """
    messups = list()
    for i in range(2, 16):  # TODO Is this hardcoded to a specific task? If so, has that caused problems running the script on other tasks? Or is <16 a minimum applying to all tasks?
        sr_val = str(series[i]
                     ).rsplit('.', maxsplit=1)[-1].rsplit(': ', maxsplit=1)[-1]
        if printit:
            print('{0}=={1} or {0}=={2}'.format(sr_val, i, i+1))
        messups.append(sr_val == str(i) or sr_val == str(i + 1))
    return all(messups)
    

def floatify(series):
    """
    :param series: pandas.Series filled with either numerical data or strings
                   representing numerical data which may have 1 too many 
                   decimal points
    :return: series, but filled completely with valid numeric data (stripped of
             of the extra decimal point and any digits after it)
    """
    try:
        result = pd.to_numeric(series)
    except ValueError:
        result = pd.to_numeric(series.apply(lambda x: x if isinstance(x, float)
                                            else x.rsplit('.', maxsplit=1)[0]))
    return result


def screen_out_blank_dependent_contrasts(cli_args):
    """
    Given a .tsv file with data for all non-blank contrasts for all subjects,
    split it into two .tsv files, with one containing all of the subjects whose
    contrasts depended on earlier blank contrasts, and the other containing
    the rest; save both into new .tsv files
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    cli_args = ensure_dict_has(cli_args, 'output', get_QC_dirpath(cli_args))

    # Get the matrix mapping each contrast number to the number(s) of the 
    # contrast(s) it depends on (if the dependency is blank, so is the other)
    dependency_matrix = pd.read_csv(os.path.join(
            SCRIPT_DIR, 'scan_info',
            'contrast_dependency_matrix_{}.csv'.format(cli_args['task'])
        ), index_col=[0])
        
    # This dict maps each contrast number to a generic name of its .tsv file
    fname_of = {x: get_generic_intensity_QC_TSV_path(cli_args, x)
                for x in cli_args['contrasts']}

    # This dict maps a contrast number to a pd.Series containing subject IDs 
    # of the subjects whose contrast for that number is blank (read in from the
    # intensity_QC_blank_EVs.tsv file for that contrast), stored for reuse
    subjects_blank_for_contrast = dict()

    # All subject-session-task indices to drop from a ready_for_stats.tsv file
    to_drop = set()

    for contrast in cli_args['contrasts']:
        if '{}_COPE'.format(contrast) not in dependency_matrix:
            dependencies = dependency_matrix.iloc[contrast - 1]
            dependencies = [int(x[0]) for x in
                            dependencies[dependencies != 0.0
                                         ].index.values.tolist()]

            ready_matrix = pd.read_csv(fname_of[contrast].format(READY1),
                                       index_col=[0], sep='\t')
            for con_dep in dependencies:  # con_dep means contrast_dependency
                subjects_blank_for_contrast = ensure_dict_has_tsv_mx_ix(
                    subjects_blank_for_contrast, con_dep, fname_of, BLANK
                )
                to_drop.update(set(subjects_blank_for_contrast[con_dep])
                                .intersection(ready_matrix.index))
                                
            print('Dropping from contrast {}: {}'.format(contrast, len(to_drop)))
            to_tsv(ready_matrix[~ready_matrix.index.isin(to_drop)],
                   fname_of[contrast].format('_ready_no_blank_dependents'))
            to_drop = set()


def ensure_dict_has_tsv_mx_ix(a_dict, a_key, fnames_dict, identifier):
    """
    Ensure that a_dict maps a_key to the first/pd.Index column of the matrix in
    the file at the path mapped to a_key in fnames_dict.
    :param a_dict: Dictionary (any)
    :param a_key: Object (any, but in this case, contrast number integers)
    :param fnames_dict: Dictionary mapping contrast numbers to generic .tsv
                        filepaths w/ unique/identifying parts replaced by '{}'
    :param identifier: String, the unique part of the generic filename
    :return: a_dict, but with a_key mapped to a pandas.Series of subject IDs
             with blank contrasts for this number, if a_key wasn't in a_dict 
    """
    if a_key not in a_dict:
        a_dict[a_key] = pd.read_csv(
                fnames_dict[a_key].format(identifier), index_col=[0], sep='\t'
            ).index
    return a_dict


def screen_out_4_stdev_outliers(cli_args):
    """
    Given a .tsv file with data for all valid contrasts for all subjects, split 
    that .tsv into two more, where one contains all outliers of the data 
    (defined as +/- 4-standard-deviations away from the mean) and the other
    contains the rest of the data; then save both into new .tsv files.
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    cli_args = ensure_dict_has(cli_args, 'output', get_QC_dirpath(cli_args))
    dt_mean_col = 'dtseries_MEAN'
    dt_4std_col = dt_mean_col + '{}_4_std'
    for contrast in cli_args['contrasts']:

        # Get input and output .tsv file names
        file_base = get_generic_intensity_QC_TSV_path(cli_args, contrast)
        file_in = file_base.format('_ready_no_blank_dependents')
        if not os.path.exists(file_in):
            file_in = file_base.format(READY1)
        outliers_out = file_base.format('_4SD_outliers')
        reduced_out = file_base.format(READY2)

        df = pd.read_csv(file_in, sep='\t', header=0, index_col=0,
                         names=['sub_basename', *COLS])

        # Calculate 4 standard deviation threshold to identify outliers        
        dtseries_MEAN_mean = df[dt_mean_col].mean()
        dtseries_MEAN_4_std = df[dt_mean_col].std() * 4

        # Add columns to df specifying whether each row is an outlier
        df.loc[df[dt_mean_col] > dtseries_MEAN_mean + dtseries_MEAN_4_std,
               dt_4std_col.format('plus')] = 1
        df.loc[df[dt_mean_col] < dtseries_MEAN_mean - dtseries_MEAN_4_std,
               dt_4std_col.format('minus')] = 1
        df = df.fillna(0)
        df[dt_4std_col.format('')] = (df[dt_4std_col.format('minus')]
                                      + df[dt_4std_col.format('plus')])

        # Save outliers df and non-outliers df into separate .tsv files
        df2 = df[df.dtseries_MEAN_4_std == 1]
        df2 = df2.drop([dt_4std_col.format(x) for x in ('plus', 'minus', '')],
                       axis='columns')
        to_tsv(df2, outliers_out)
        df = df[df.dtseries_MEAN_4_std != 1]
        df = df.drop([dt_4std_col.format(x) for x in ('plus', 'minus', '')],
                     axis='columns')
        to_tsv(df, reduced_out)


def make_subjects_copes_chart(cli_args):
    """
    Given that each subject has a .tsv file containing data for each of its
    contrasts, read in all contrast .tsv files for all subjects and create
    a chart showing which subjects have valid data for which contrasts and
    which do not; then save that chart into a new .tsv file.
    :param cli_args: Dictionary containing all command-line arguments from user
    """
    cli_args = ensure_dict_has(cli_args, 'output', get_QC_dirpath(cli_args))
    all_subjects = pd.DataFrame(index=[os.path.basename(x) for x in glob(
            os.path.join(cli_args['study_dir'], 'sub-NDARINV*')
        )], columns=cli_args['contrasts'])
    for contrast in cli_args['contrasts']:
        subjects_with_valid_contrast = list()
        file_base = get_generic_intensity_QC_TSV_path(cli_args, contrast)
        contrast_subjs = pd.read_csv(file_base.format(READY2), sep='\t')

        contrast_subjs['subjectID'] = contrast_subjs.apply(
            lambda row: row['sub_basename'].split('_', 1)[0], axis=1
        )
        subjects_with_valid_contrast = contrast_subjs['subjectID'
                                                      ].values.tolist()
        all_subjects[contrast].loc[subjects_with_valid_contrast] = 1
    all_subjects.fillna(0).to_csv(os.path.join(
        cli_args['output'],
        cli_args['task'] + '_all_subjects_which_contrasts.tsv'
    ), index=True, header=True, sep='\t')


if __name__ == '__main__':
    main()