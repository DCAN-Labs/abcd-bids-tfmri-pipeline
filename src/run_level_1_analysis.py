#!/usr/bin/env python3
# coding: utf-8

"""
Run ABCD Task Prep Level 1 analysis for pipeline_wrapper.py
Greg Conan: gconan@umn.edu
Created: 2021-01-15
Updated: 2021-12-03
"""

# Import standard libraries
from glob import glob
import os
import shutil
import subprocess
import sys


def find_local_imports(flg):
    """
    Ensure that this script can find its dependencies if parallel processing
    """
    try:
        parallel_flag_pos = -2 if sys.argv[-2] == flg else sys.argv.index(flg)
        sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))
    except (IndexError, ValueError):
        sys.exit('Please include the {} flag followed by the path to the '
                 'directory containing pipeline_utilities.py'.format(flg))


# Local custom imports
WRAPPER_LOC = '--wrapper-location'
find_local_imports(WRAPPER_LOC)
from src.pipeline_utilities import (
    add_lvl_args_to, add_slurm_args_to, copy_and_rename_file,
    extract_from_json, get_args_to_run_film_gls, get_LR_functions,
    get_lvl_paths, get_optional_cli_args, get_pipeline_cli_argparser,
    get_region_path_vars, get_replacements, get_sub_base, get_TR_and_ntpts,
    glob_and_copy, make_and_save_confound_matrix, make_fake_nifti,
    merge_files_in_range, merge_to_make_dtseries, organize_lvl_paths,
    overwrite_dirs, rename_template_file_vars, run_fsl_tool,
    run_parallel_or_sequential, validate_cli_args, wb_command, wb_LR_pair
)


def main():
    # Get command line arguments and variables from main pipeline wrapper
    cli_args = _cli()
    json_dict = extract_from_json(cli_args['temp_json'])
    sigma = json_dict['sigma']
    lr_fns = get_LR_functions(cli_args, json_dict['paths'])

    # Get/build Level-1-specific paths and then prepare for Level 1 analysis
    hemi_py_script = os.path.join(cli_args['wrapper_location'], 'src',
                                  'process_hemisphere.py')
    paths = add_lvl_1_paths_to(json_dict['paths'], cli_args)
    
    dtseries, scan = get_dtseries_paths(cli_args, paths, lr_fns, sigma)
    (subcort_vol, confound_matrix
     ) = get_vol_and_confound_mx_paths_after_wb_prep(cli_args, paths, dtseries)

    # Level 1 analysis
    get_events_make_template(cli_args, paths, confound_matrix, dtseries, scan)
    subcort_stats = run_fsl_for_hemisphere_and_subcort(
        cli_args, paths, confound_matrix, subcort_vol, hemi_py_script
    )
    merge_files_to_make_all_dtseries(cli_args, paths, lr_fns, subcort_stats)


def _cli():
    """
    :return: argparse.ArgumentParser with all needed command-line arguments 
    """
    parser = add_slurm_args_to(add_lvl_args_to(get_pipeline_cli_argparser()))
    return validate_cli_args(vars(parser.parse_args()), parser)


def add_lvl_1_paths_to(paths, cli_args, *extra_lvl_1_dirs):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :return: paths, but with level 1 paths for this task added
    """
    paths['sub_run'] = get_sub_base(cli_args, cli_args['run_number'])
    paths['lvl_1'] = get_lvl_paths(paths['dir_lvl']['1'], paths['sub_run'],
                                   paths['feat_name'], list(), 'fsf', 'EV',
                                   'intermediate', *extra_lvl_1_dirs)
    paths['lvl_1']['GOStats'] = os.path.join(paths['lvl_1']['parent'],
                                             'GrayordinatesStats')
    overwrite_dirs(organize_lvl_paths(paths['lvl_1']), mkdir=True)
    return paths


def get_dtseries_paths(cli_args, paths, lr_fns, sigma_smooth):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries thereof
    :param lr_fns: Dictionary with 2 local functions to get surface
                   and shape paths for both hemispheres
    :param sigma_smooth: Dictionary with sigma value floats for more smoothing
    :return: Dictionary containing .dtseries.nii paths
    """
    # Create smoothed dtseries file, or copy regular dtseries into its place
    dtseries = {0: get_dt_str(paths['sub_ses']['func'], paths['sub_run'])}
    dtseries['smoothed'] = get_dt_str(paths['lvl_1']['intermediate'],
                                      paths['sub_run'], paths['final_smooth'])
    smooths = ('surf', 'vol')
    if cli_args['spat_smooth'] > 0 and (cli_args['spat_smooth'] - cli_args[x]
                                        != 0 for x in smooths):
        sigmas = [str(sigma_smooth[x]) for x in smooths]
        wb_command(cli_args, '-cifti-smoothing', dtseries[0],
                   *sigmas, 'COLUMN', dtseries['smoothed'],
                   *wb_LR_pair(lr_fns['surf'], 'surface'))
        print('Smoothing completed.')
    else: 
        copy_and_rename_file(dtseries[0], dtseries['smoothed'])
        print('Continuing on without smoothing being applied...')

    # Get repetition time and number of timepoints
    scan = dict()
    scan['TR'], scan['ntpts'] = get_TR_and_ntpts(dtseries[0],
                                                 cli_args['wb_command'])
    # Create fake NIFTI files and get dtseries paths used by other functions
    dtseries['genericfakenii'] = os.path.join(paths['lvl_1']['intermediate'], (
        '{}_bold_timeseries{}_{}fake_nifti.nii.gz'
        .format(paths['sub_run'], paths['final_smooth'], '{}')
    ))
    dtseries['smoothed_fake_nifti'] = make_fake_nifti(
        cli_args, dtseries['genericfakenii'], dtseries['smoothed'], '',
        cli_args['wb_command'] + ' -cifti-convert -to-nifti {} {}'
    )
    dtseries['smoothed_mean_fake_nifti'] = make_fake_nifti(
        cli_args, dtseries['genericfakenii'], dtseries['smoothed_fake_nifti'],
        'mean_', 'fslmaths {} -Tmean {}'
    )
    dtseries['smooth_generic'] = os.path.join(
        paths['lvl_1']['intermediate'], '{}{}_bold_timeseries{}_filtered{}'
    )
    dtseries['smooth_filtered_fake_nifti'] = make_fake_nifti(
        cli_args, dtseries['smooth_generic'].format(
            paths['sub_run'], '', paths['final_smooth'], '{}'
        ) + '.nii.gz', dtseries['smoothed_fake_nifti'], '_fake_nifti',

        # Apply high pass temporal filter via fslmaths & add mean back to image
        'fslmaths {} -bptf {} -1 -add {} {}',

        # Compute smoothing kernel sigma 
        ((0.5 * cli_args['filter']) / float(scan['TR'])),
        dtseries['smoothed_mean_fake_nifti']

    # Skip temporal filtering if --filter is 0
    ) if cli_args['filter'] > 0 else dtseries['smoothed_fake_nifti']
    dtseries['smoothed_filtered'] = dtseries['smooth_generic'].format(
        paths['sub_run'], '', paths['final_smooth'], ''
    ) + '.dtseries.nii'
    return dtseries, scan


def get_dt_str(dirname, sub_ID, smoothing=''):
    """
    Build .dtseries.nii file path by filling in a generic string's unique parts
    :param dirname: String, a valid path to the real parent directory
    :param sub_ID: String identifying the subject whose dtseries to get
    :param smoothing: String showing (if present) that the dtseries is smoothed
    :return: String, the full path to a .dtseries.nii file
    """
    return os.path.join(dirname, ('{}_bold_timeseries{}.dtseries.nii'
                                  .format(sub_ID, smoothing)))


def get_vol_and_confound_mx_paths_after_wb_prep(cli_args, paths, dtseries):
    """
    Make, convert, and copy Level 1 feat input files to where they belong, and
    then create confound matrix
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param dtseries: Dictionary containing .dtseries.nii file path strings
    :return: Tuple of 2 strings; paths to subcortical volume & confound matrix
    """
    # Create the input files for Level 1 feats for each run
    wb_command(cli_args, '-cifti-convert', '-from-nifti',
               dtseries['smooth_filtered_fake_nifti'],
               dtseries['smoothed'], dtseries['smoothed_filtered'])
    subcort_vol = os.path.join(paths['lvl_1']['intermediate'],
                              '{}_AtlasSubcortical{}_filtered.nii.gz'
                              .format(paths['sub_run'], paths['final_smooth']))
    hemi = os.path.join(paths['lvl_1']['intermediate'],
                        paths['sub_run'] + paths['final_smooth']
                        + '_filtered.atlasroi.{}.32k_fs_LR.func.gii')
    wb_command(cli_args, '-cifti-separate-all', dtseries['smoothed_filtered'],
               '-volume', subcort_vol, *wb_LR_pair(lambda lr: hemi.format(lr)))
        
    # Copy fsf template to fsf dir and motion files to intermediate dir
    shutil.copy(paths['template1'], paths['lvl_1']['fsf'])
    desc_tsv_file = paths['sub_run'] + '_desc-filteredincludingFD_motion.tsv'
    glob_and_copy(paths['lvl_1']['intermediate'], paths['sub_ses']['func'],
                  desc_tsv_file)

    # Make confound matrix, copy it subject's fsf_paths dir, and return it
    confounds_name = make_and_save_confound_matrix(
        cli_args, desc_tsv_file, paths['lvl_1'], paths['sub_run']
    )
    return subcort_vol, os.path.join(paths['lvl_1']['fsf'], confounds_name)


def get_events_make_template(cli_args, paths, confound_mx, dtseries, scan):
    """
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param confound_mx: String, path to confounds matrix file
    :param dtseries: Dictionary containing .dtseries.nii file path strings
    :param scan: Dictionary with two values, the repetition time and number of
                 timepoints for this subject, session, task, and run
    """
    # If user-given event files dir has events for this subj/ses/run, use those
    event_file_base = (cli_args['events_dir'], paths['sub_run'] + '*.tsv')
    event_fpaths = glob(os.path.join(*event_file_base))
    if any(event_fpaths):
        glob_and_copy(paths['lvl_1']['EV'], *event_file_base)

    # Create the subject- and run-specific fsf file for Level 1 analyses
    fsf_generic = os.path.join(paths['lvl_1']['fsf'], '{}')
    fsf_template = fsf_generic.format(os.path.basename(paths['template1']))
    run_N = 'run-{}'.format(cli_args['run_number'])

    # Create new .fsf file by copying level 1 template & renaming its variables
    new_fsf_template = fsf_generic.format("{}_{}_{}_{}_level1.fsf".format(
        cli_args['subject'], cli_args['ses'], 'task-' + cli_args['task'], run_N
    ))
    rename_template_file_vars(fsf_template, new_fsf_template, get_replacements(
        cli_args, CONFOUND_MATRIX=confound_mx, RUNNUM=run_N,
        NTPTS=str(scan['ntpts']), REP_TIME=str(scan['TR']),
        DT_SF_FAKE_NIFTI=dtseries['smooth_filtered_fake_nifti']
    ))
    # Verify that all necessary event files have been created
    event_fpaths.append(confound_mx)  # <-- Prevents verify_event_files from 
    try:                              #   calling confound_mx a missing event
        verify_event_files(cli_args['task'], new_fsf_template, event_fpaths)
    except (FileNotFoundError, ValueError):
        err_type, err_msg, _ = sys.exc_info()
        sys.exit('\n{}: {}\n'.format(err_type.__name__, err_msg))


def verify_event_files(task_name, template_fsf_path, event_fpaths):
    """
    :param task_name: String naming the --task
    :param template_fsf_path: String, path to .fsf file with no placeholders 
    :param event_fpaths: List of strings, each a path to an extant event file
    """
    # Check that all event files listed in the template .fsf file exist in
    # the event files directory supplied by the user
    event_fnames = set(os.path.basename(ev_fpath) for ev_fpath in event_fpaths)
    ev_paths = list()
    with open(template_fsf_path) as infile:
        for eachline in infile:
            eachline = eachline.strip()
            if '.tsv' in eachline[-5:]:
                ev_paths.append(eachline.rsplit(' ', 1)[-1].strip('"'))
                basename = os.path.basename(ev_paths[-1])
                if basename not in event_fnames:
                    raise FileNotFoundError('\n{} event file not found at {}'
                                            .format(basename, ev_paths[-1]))
    # Warn user if they are using a template .fsf file for the wrong task
    if ev_paths and not any(task_name in ev_name for ev_name in ev_paths):
        raise ValueError('None of the event files in the template .fsf file '
                         '{} are named with your --task {}. You may have '
                         'chosen the wrong --template1 file. Please use '
                         'another with event file names including your --task.'
                         .format(template_fsf_path, task_name))


def run_fsl_for_hemisphere_and_subcort(cli_args, paths, confound_matrix,
                                       subcortical_volume, hemi_py_script):
    """
    Generate subject's .con and .mat files, run feat_model for confound matrix,
    then run hemisphere and subcortical processing--in parallel or sequentially
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param hemi_py_script: String, valid path to existing Python script to run
                           for hemispheric processing (parallel or not)
    :return: String, the path to the subcortical statistics directory
    """
    design, _, subcort_stats, _ = get_region_path_vars(cli_args, paths, 
                                                       cli_args['run_number'])
    run_fsl_tool(cli_args, 'feat_model', design, confound_matrix)
    to_replace = '{replace_here}'
    film_args_subcort = {'in_arg': subcortical_volume, 'rn': subcort_stats,
                         'con': design + '.con', 'thr': 1, 
                         'pd': design + '.mat', 'mode': 'volumetric'}
    if os.path.exists(design + '.fts'):
        film_args_subcort['fcon'] = design + '.fts'
    film_args_subcort = get_args_to_run_film_gls(**film_args_subcort)
    if not cli_args['no_parallel']:
        film_args_subcort = [film_args_subcort]
    run_parallel_or_sequential(
        hemi_py_script, cli_args, ('L', 'R'), to_replace,
        [*get_optional_cli_args(cli_args, True), '--hemisphere', to_replace],
        second_fn=subprocess.check_call, second_args=film_args_subcort
    )
    return subcort_stats


def merge_files_to_make_all_dtseries(cli_args, paths, lr_fns, subcort_stats):
    """
    Merge files to create their .dtseries.nii files by running wb_command
    -cifti-create-dense-timeseries on files occurring (#1) only once or (#2) as
    often as feat model contrasts; (#3) parameter estimates, & (#4) f-testing
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param subcort_stats: String, the path to the subcortical stats directory
    """
    for fname in ('sigmasquareds', 'threshac1', 'res4d'):
        merge_to_make_dtseries(cli_args, fname, paths['lvl_1'],              #1
                               subcort_stats, paths['AROI2'], lr_fns['shape'])
    files_count = {f: f + '*.nii.gz' for f in ('pe', 'cope', 'fstat')}
    for label, glob_base in files_count.items():
        files_count[label] = len(glob(os.path.join(subcort_stats, glob_base)))
    others = (paths['lvl_1'], subcort_stats, paths['AROI2'], lr_fns['shape'])
    merge_files_in_range(cli_args, ['cope', 'tstat', 'zstat', 'varcope'],    #2
                         files_count['cope'], others)
    merge_files_in_range(cli_args, ['pe'], files_count['pe'], others)        #3
    if files_count['fstat'] > 0:
        merge_files_in_range(cli_args, ['fstat', 'zfstat'],                  #4
                             files_count['fstat'], others)


if __name__ == '__main__':
    main()