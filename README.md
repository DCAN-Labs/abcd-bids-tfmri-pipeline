# abcd-bids-tfmri-pipeline

## Purpose
This tool is a modified version of the TaskfMRIAnalysis stage of the HCP-pipeline (Glasser et al., 2013) developed at University of Vermont by Anthony Juliano and can be used to complete Level 1 and Level 2 analyses of task fMRI dtseries data. Specifically, this tool was designed to work with data that were minimally processed using the [DCAN Labs' ABCD-HCP-Pipeline](https://github.com/DCAN-Labs/abcd-hcp-pipeline), including their [Collection 3165 release through NDA](https://collection3165.readthedocs.io/en/stable/). Inputs must be in the same format as ABCD-HCP-Pipeline outputs after running filemap. The files output from this pipeline were fully processed and prepared for higher-level statistics.

## Overview 
The abcd-bids-tfmri pipeline is optimized for HCP-style CIFTI files and thus heavily relies on [HCP workbench commands](https://www.humanconnectome.org/software/workbench-command). This includes user-specified spatial smoothing (wb_command -cifti-smoothing), converting the smoothed data to and from a format that FSL (Jenkinson et al. 2012) can interpret (wb_command -cifti-convert), separating the dtseries data into its comprised components (wb_command -cifti-separate-all), and reading in pertinent information from the dtseries data (wb_command -file-information), among others. Based on the user-specified parameters for censoring volumes (i.e. initial and/or high-motion frames), the pipeline will read in the filtered motion file (Fair et al., 2020) produced by the ABCD-BIDS processing pipeline and create a matrix for nuisance regression. Finally, high-pass filtering, with a cutoff of 0.005 Hz (200 seconds), is completed before running FSL's FILM (Woolrich et al. 2001).

For FILM to run, users must supply their own subject-, task-, and run-specific event timing files that are in the FSL standard three column format (i.e. onset, duration, weight/magnitude). Additionally, users need to supply a task-specific fsf template file per task that they will be processing using the abcd-bids-tfmri pipeline. As the abcd-bids-tfmri pipeline modifies this template to make it subject- and run-specific, certain values need to be replaced with specific variables that the abcd-bids-tfmri pipeline will be able to recognize. An example fsf file template for ABCD’s MID task is made available for users to review on ABCC (https://osf.io/psv5m/).

Users can specify which task data they would like to process by providing a list of task names within the abcd-bids-tfmri pipeline’s command line interface. If the user specifies multiple runs of the task, the pipeline will complete higher-level analyses (i.e. fixed effects modeling) to combine a given subject's run-level data. Therefore if a study has three different fMRI tasks that consist of two runs, all six level 1 analyses and all three level 2 analyses can be completed for a subject with a single run of the abcd-bids-tfmri pipeline.

## Requirements and Dependencies

1. [Python 3.6.13](https://www.python.org/downloads/release/python-3613/) or greater
1. Washington University [Workbench Command (wb_command)](https://github.com/Washington-University/workbench)
1. [FMRIB Software Library (FSL) v6.0.4](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
1. Python [`numpy` package](https://numpy.org)
1. Python [`pandas` package](https://pandas.pydata.org)
 
## Usage

### Installation

1. Clone this repository into a new directory on your local filesystem.
1. Install the required Python modules using `pip`.
1. To verify that the code is set up, run `python3 pipeline_wrapper.py --help` within that same directory.

Once this repository is cloned into a local directory, and all of its dependencies are also downloaded onto the local filesystem, no further installation or setup is required. Running this pipeline on multiple subjects at once will require a SLURM cluster.

### Required Arguments

- `--subject` takes one string, the ID of the subject to process. If it does not start with the prefix `sub-`, that prefix will be added.

- `--ses` takes one string, the name of the session that you are running the pipeline on. If it does not start with the prefix `ses-`, that prefix will be added.

- `--task` takes one string, the name of the task that you are running the pipeline on.

- `--study-dir` takes a valid path to the existing BIDS base study directory. The `--study-dir` should have a subdirectory called `derivatives`, which should have a subdirectory called `abcd-hcp-pipeline` (or called something else, which the user can provide with the `--bids-dir` argument).

- `--events-dir` takes a valid path to a real directory containing event `.tsv` files. `--events-dir` is required unless running a certain kind of level 1 analysis. This flag is necessary unless the level 1 `.fsf` file specifies that the basic waveform shape of each EV model is square. Include this flag unless the level 1 `.fsf` file contains the line `set fmri(shape1) 0`.

- `--output` takes a valid directory path to save pipeline outputs into. If the provided `--output` directory does not already exist, then this script will create a new directory at that path.

- `--wrapper-location` takes a valid path to this directory containing `pipeline_wrapper.py`.
  - Currently, an error may occur if this argument is excluded. In future versions of the pipeline, including this argument will not be necessary unless running the pipeline (or its sub-scripts) in parallel using `SLURM`.

### Example Call

To run this pipeline on one subject from the Linux/Unix command-line terminal, use a command that looks like this:

```sh
fsl_dir=/home/software/fsl/bin/
subject_ID=sub-ABCDEFGH
session_ID=ses-QRSTUV
task_name=my_task
study_dir=/home/users/shared/data/my_study
wb_command=/home/software/workbench/bin/wb_command
events_dir=/home/users/shared/data/task_events/
wrapper_dir=.

python3 pipeline_wrapper.py --subject ${subject_ID} --ses ${session_ID} --study-dir ${study_dir} --task ${task_name} --events-dir ${events_dir} --fsl-dir ${fsl_dir} --wb-command ${wb_command} --wrapper-location ${wrapper_dir}
```

## Other Arguments/Options

### Arguments With Default Values
    
Often, the following arguments are unnecessary because they have default values that the pipeline will use. But the user can change them by including these arguments explicitly.<br />

### Files for GLMs

  - `--template1` accepts a string, the base name (but not the full path) of the Level 1 `.fsf` template file. By default, the pipeline will assume that the `--template1` filename is `template_DCAN_version_`(`--task`)`_level1_UPDATED_FINAL.fsf`.
  - `--template2` accepts a string, the base name (but not the full path) of the Level 2 `.fsf` template file. By default, the pipeline will assume that the `--template2` filename is `template_DCAN_version_`(`--task`)`_level2_UPDATED_FINAL.fsf`.

### Information for GLMs

#### Runs 

  - `--runs` takes 1 or more positive integers, each subject's number of runs. The integers must be provided as a space-delimited list, such as `1 2 3 4` for runs 1 through 4. By default, this argument's value(s) will be `1 2`.
  - `--levels` takes up to 2 positive integers, the levels to conduct the analysis on: `1` for one run, and/or `2` to merge multiple runs. By default, both levels will be run.

#### Censoring

  - `--censor` takes a positive integer, the number of initial frames/volumes to censor. By default, this argument's value will be `2`.
  - `--fd` takes one decimal number between 0 and 1, the framewise displacement threshold for censoring volumes with high motion. By default, this argument's value will be `0.9`.

#### Smoothing

  - `--spat-smooth` takes one positive integer, the number of millimeters of spatial smoothing that you want for the surface and volume data. By default, this argument's value will be `0`.
  - `--surf-smooth` takes one positive integer, the number of millimeters of surface smoothing that has already been applied in the minimal processing steps. By default, this argument's value will be `0`.
  - `--vol-smooth` takes one positive integer, the number of millimeters of volume smoothing that has already been applied in the minimal processing steps. By default, this argument's value will be `0`.

#### Highpass Filtering

  - `--filter` takes one positive integer, the high pass temporal filter cutoff (in seconds). By default, this argument's value will be `100`.

#### Paths

  - `--bids-dir` takes one string, the name of the BIDS-standard file-mapped directory with subject data in the "derivatives" subdirectory of your --study-dir. This argument is only needed if your BIDS directory is called something other than `abcd-hcp-pipeline`.
    - This path should be valid: `(--study-dir)/derivatives/(--bids-dir)/(--subject)/(--ses)/func/sub-(--subject)_ses-(--ses)_task-(--task)_run-(--runs)_bold_timeseries.dtseries.nii`
  - `--temp-dir` takes a valid path to a directory to save temporary files into. By default, the pipeline will look for a folder called `temp` in the `level-1` subdirectory of the `--output` folder. If no folder exists at the `--temp-dir` path, then the pipeline will create an empty folder there.
  - `--templates` takes a valid path to an existing directory which contains all of the required template `.fsf` files. By default, the pipeline will use a folder called `templates/` in the same file location as the main pipeline script (`pipeline_wrapper.py`), the top level of this repository.

#### Naming Conventions

  - `--study-name` takes one string, the name of the study that you are running the pipeline on. If no `--study-name` is given by the user, then the pipeline will give the outputs the default study name of `ABCD`.
    
### Optional Arguments for Different Run Modes

Each of these flags is unnecessary, but if included, will change the mode in which the pipeline runs.<br />

  - `--help` takes no parameters. Include this flag to show this help message and exit.
  - `--keep-all` takes no parameters. Include this flag to keep all files generated during the pipeline. By default, the pipeline will only keep `.dtseries.nii`, `dof`, `log`, and event files.
  - `--overwrite` takes no parameters. Include this flag to overwrite files at the `--output` directory. By default, the script will end if files exist there.
  - `--no-parallel` takes no parameters. By default, the pipeline script will use Python's `multiprocessing` module to process level 1 analysis runs in parallel simultaneously. Include this flag for the script to process the analyses sequentially instead. Using this flag is not recommended when the user can run the pipeline on multiple cores, because running the analyses sequentially will significantly increase the pipeline's operating time.

### Quasi-Required Arguments
    
These arguments are usually required, except in special circumstances.<br />

These arguments are paths to software libraries required to run the task pipeline. By default, the pipeline will try to infer the paths. But if that fails, then the user needs to give them explicitly. So, `--fsl-dir` and `--wb-command` are required unless the user has certain paths already saved as BASH environment variables.

- `--fsl-dir` takes a valid path to the existing FreeSurfer `bin` directory containing the executable files `fsl`, `fslmerge`, `fslmaths`, `flameo`, and `feat_model` from the FMRIB Software Library (FSL). If this flag is excluded, then the script will try to guess the path to the FSL by checking the user's BASH aliases.
- `--wb-command` takes a valid path to the existing `wb_command` executable file to run Workbench Command. If this flag is excluded, then the script will try to guess the path to the `wb_command` file by checking the user's BASH aliases.


#### SLURM/SBATCH Arguments for run_sbatch_jobs.py
  
`pipeline_wrapper.py` is the main script to run the pipeline on one task of one session of one subject. Use `run_sbatch_jobs.py` to run multiple subjects/sessions/tasks at once using a SLURM cluster. The arguments below are not needed to run 1 subject using `pipeline_wrapper.py` but are needed to run multiple subjects in parallel using `run_sbatch_jobs.py`

  - `--account` takes one string, the name of the account to submit each `SBATCH` job under.
  - `--cpus` takes one integer, the number of CPUs to use for each `SBATCH`/`SLURM` job. By default, this argument's value will be `1`.
  - `--memory` takes one integer, the memory in gigabytes (GB) to assign to each sbatch job. The default number is `8`.
  - `--print-progress` takes no parameters. Include this flag for the script to print updates about its progress at intervals defined by `--sleep`. This will also print every command that is run to submit a pipeline batch job.
  - `--sleep` takes one integer, the number of seconds to wait between batch job submissions. The default number is `10`.
  - `--time` takes a string formatted specifically as `HH:MM:SS` where `HH` is hours, `MM` is minutes, and `SS` is seconds. This argument is the time limit for each ABCD-BIDS-task-fmri-pipeline `SBATCH`/`SLURM` job. `01:00:00` is the default time limit.

</details>

## Current Issues

The current version of the pipeline might crash if the user is not aware of how to avoid certain pitfalls. Future versions of the pipeline will aim to fix all of these issues.

### An error may occur unless...

- The `--output` directory matches the location of the `.fsf` template file FEAT directory containing the subdirectories `level-1` and `level-2`.
- The `--events-dir` location matches the `level-1` subdirectory of the `.fsf` template file FEAT directory.
- The `--template[x]` flags are filenames only, and not full file paths. 
- The `--wrapper-location` argument is explicitly included.
    - Although the script is supposed to infer its own location (except when run using SBATCH), there may still be cases where the script needs its location explicitly passed in.
- Only one task per subject session is run at a time, or `--keep-all` is explicitly included.
    - Without the `--keep-all` flag, running the pipeline on one task may overwrite/delete the outputs of a simultaneous pipeline run on another task.

### Future Features

Future versions of the pipeline will hopefully

  - be able to run multiple tasks in one script call

  - not need `--wrapper-location` unless running `run_sbatch_jobs.py`

## Outputs

The pipeline will produce its outputs in the following [BIDS](http://bids.neuroimaging.io/)-valid directory structure:

```
output_dir
├── level-1
│   ├── events
│   │   └── sub-*_ses-*_task-*_run-*.tsv
│   ├── level1_run-*
│   │   ├── sub-*_ses-*_task-*_run-*_cope*.dtseries.nii
│   │   ├── sub-*_ses-*_task-*_run-*_dof
│   │   ├── sub-*_ses-*_task-*_run-*_logfile
│   │   └── sub-*_ses-*_task-*_run-*_res4d.dtseries.nii
│   └── temp
└── level-2
    ├── cope_files
    │   └── sub-*_ses-*_task-*_contrast_*_cope*.dtseries.nii
    ├── dof_files
    │   └── sub-*_ses-*_task-*_contrast_*_tdof_t1.dtseries.nii
    ├── log_files
    │   └── sub-*_ses-*_task-*_contrast_*_logfile
    ├── mask_files
    │   └── sub-*_ses-*_task-*_contrast_*_mask.dtseries.nii
    └── res4d_files
        └── sub-*_ses-*_task-*_contrast_*_res4d.dtseries.nii
```
 
## Metadata

Information about this `README` file:

- Created on 2020-10-22 by Eric Feczko ([@EricFeczko](https://github.com/EricFeczko)) 
- Updated on 2021-12-03 by Greg Conan ([@GregConan](https://github.com/GregConan))
