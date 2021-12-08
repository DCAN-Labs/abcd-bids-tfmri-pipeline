#!/usr/bin/env python3
# coding: utf-8

"""
Hemisphere analysis for run_level_1_analysis.py parallel processing 
Greg Conan: gconan@umn.edu
Created: 2021-01-29
Updated: 2021-09-20
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
    add_lvl_args_to, extract_from_json, get_args_to_run_film_gls, 
    get_lvl_paths, get_pipeline_cli_argparser, get_region_path_vars,
    get_sub_base, get_subj_ses, wb_command
)


def main():
    cli_args = _cli()  # Get command-line input arguments

    # Get all necessary paths from JSON
    paths = extract_from_json(cli_args['temp_json'])['paths']
    paths['lvl_1'] = get_lvl_paths(
        paths['dir_lvl']['1'], get_sub_base(cli_args,
                                            cli_args['run_number']),
        paths['feat_name'], [cli_args['run_number']],
        'abcd_extract', 'EV', 'fsf', 'intermediate'
    )
    paths['lvl_1']['GOStats'] = os.path.join(paths['lvl_1']['parent'],
                                             'GrayordinatesStats')

    # Run level 1 hemispheric processing
    process_hemisphere(cli_args, paths, cli_args['run_number'])


def _cli():
    """
    :return: argparse.ArgumentParser with all command-line arguments 
             needed to run this script
    """
    parser = add_lvl_args_to(get_pipeline_cli_argparser())
    parser.add_argument('--hemisphere', choices=('L', 'R'), required=True)
    return vars(parser.parse_args())


def process_hemisphere(cli_args, paths, run):
    """
    Run and complete all level 1 processing for one hemisphere
    :param cli_args: Dictionary containing all command-line arguments from user
    :param paths: Dictionary of path strings, and of dictionaries of path
                  strings, used throughout processing in both levels
    :param run: Whole number (as an int or a string) defining which run this is
    """
    # Define brain region dirpath variables from run_level_1_analysis.py
    (design, func_str, subcort_stats, surf_str
     ) = get_region_path_vars(cli_args, paths, run)

    # Dilate the smoothed cortical func.gii's for each hemisphere and each run
    func = func_str.format(get_sub_base(cli_args, run),
                           paths['final_smooth'], '', cli_args['hemisphere'])
    surf = surf_str.format(get_subj_ses(cli_args), cli_args['hemisphere'])
    dil_f = func_str.format(get_sub_base(cli_args, run), 
                            paths['final_smooth'], '_dil',
                            cli_args['hemisphere'])
    wb_command(cli_args, '-metric-dilate', func, surf, '50', dil_f, '-nearest')

    # Run FSL's film on cortical surface data for each hemisphere for each run
    cort_stats = os.path.join(paths['lvl_1']['parent'],
                              '{}_SurfaceStats'.format(cli_args['hemisphere']))
    film_args = {'in_arg': dil_f, 'rn': cort_stats, 'pd': design + '.mat',
                 'ms': 15, 'con': design + '.con', 'mode': 'surface',
                 'in2': surf, 'epith': 5}
    if os.path.exists(design + '.fts'):
        film_args['fcon'] = design + '.fts'
    subprocess.check_call(get_args_to_run_film_gls(**film_args))

    # Copy log and dof files to GrayordinatesStats directory
    with open(os.path.join(paths['lvl_1']['GOStats'], 'logfile'), 'w+') as out:
        for fname in [os.path.join(stats_dir, 'logfile') for stats_dir
                      in (subcort_stats, cort_stats)]:
            with open(fname) as infile:
                out.write(infile.read())
    for dof_file in glob(os.path.join(subcort_stats, 'dof')):
        shutil.copy(dof_file, paths['lvl_1']['GOStats'])
    

if __name__ == '__main__':
    main()