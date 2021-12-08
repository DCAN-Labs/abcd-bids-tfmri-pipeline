# ABCD-BIDS Task fMRI Pipeline: Other Scripts

The files in this `src` directory are scripts which supplement the main pipeline script (`pipeline_wrapper.py`), using much of the same functionality:

## Overview

### Brief Summary of Each Script in This Directory

- `concatenate_EV_files.py` takes a list of many subjects' event `.tsv` files and concatenates all event `.tsv` files (per run per task per session per subject) into one `all_events.tsv` file which can be used by the main pipeline script to produce main effect contrast files.
- `detect_outliers_and_QC.py` filters out invalid and outlier event `.tsv` files from those used by the main pipeline script, and shows which COPE (contrast parameter estimate) files produced by the main pipeline script are valid.
- `get_events_extract_eprime.py` is an [ABCD](https://abcdstudy.org/)-specific, [MIDB](https://midb.umn.edu/)-/[DCAN-Lab](https://innovation.umn.edu/developmental-cognition-and-neuroimaging-lab/)-specific script used to extract the event/EV `.tsv` files required for `pipeline_wrapper.py` from the `sourcedata` outputs of the [ABCD-HCP-Pipeline](https://github.com/DCAN-Labs/ABCD-HCP-pipeline). It also verifies that the event files match the (ABCD-specific) fMRI scanner information.
- `make_script_list.py` creates a text file where each line is a command to run `pipeline_wrapper.py` on a different subject and/or a different task. To run many instances of `pipeline_wrapper.py` in parallel, the user must (1) run this script to create that text file, and then (2) run the `run_sbatch_jobs.py` script to submit each line of the text file as its own `SBATCH` job.
- `pipeline_utilities.py` does nothing on its own, but it contains all of the functionality used by multiple scripts in this repository. It defines most of this codebase's lower-level tools and utility functions.
- `process_hemisphere.py` does the level 1 analysis on 1 brain hemisphere (from 1 run of 1 session of 1 subject).
- `run_level_1_analysis.py` does the entire level 1 analysis on 1 run (of 1 session of 1 subject). Level 1 analysis uses FSL and Workbench to process both hemispheres and subcortical data. This script runs 2 instances of `process_hemisphere.py`, 1 per hemisphere. It also does the level 1 analysis on subcortical data. Although a user can run this script on its own if the user provides a `--temp-json` file, it is usually easier to run `pipeline_wrapper.py`, which automatically runs `run_level_1_analysis.py` on each run for 1 session of 1 subject. 

## Explanatory Charts

<details><summary style="text-align: center;"><strong style="display: inline; font-size: 1.17em;">Table: Which Scripts Accept Which Arguments</strong></summary>

| Argument Name | ...wrapper.py | run_level... | process... | concatenate... | detect... | get_events... |
|----------------------------|:-:|:-:|:-:|:-:|:-:|:-:|
| `--bids-dir`         | Optional | Optional | Optional | Optional | | |
| `--censor`           | Optional | Optional | Optional | Optional | | |
| `--contrasts`        | | | | | Optional | |
| `--events-dir`       | Usually **Required**<sup> 1</sup> | Usually **Required**<sup> 1</sup> | Usually **Required**<sup> 1</sup> | Usually **Required**<sup> 1</sup> | | Usually **Required**<sup> 1</sup> |
| `--fd`               | Optional | Optional | Optional | Optional | | Optional |
| `--filter`           | Optional | Optional | Optional | Optional | | |
| `--fsl-dir`          | Usually **Required**<sup> 3</sup> | Usually **Required**<sup> 3</sup> | Usually **Required**<sup> 3</sup> | Usually **Required**<sup> 3</sup> | |
| `--hemisphere`       | | | **Required** | | | | |
| `--keep-all`         | Optional | Optional | Optional | Optional | | Optional |
| `--levels`           | Optional | Optional | Optional | Optional | | |
| `--output`           | **Required** | **Required** | **Required** | **Required** | **Required** | **Required** |
| `--parallel`         | Optional | Optional | Optional | Optional | | |
| `--pipeline`         | | | | Optional | | |
| `--print-progress`   | Optional | Optional | | | |
| `--runs`             | Optional | Optional | Optional | Optional | | Optional |
| `--run-number`       | | **Required** | **Required** | | | |
| `--scanners-info`    | | | | Optional | | Optional |
| `--ses`              | **Required** | **Required** | **Required** | **Required** | **Required** | **Required** |
| `--sleep`            | Optional | Optional | | | | |
| `--spat-smooth`      | Optional | Optional | Optional | Optional | | Optional |
| `--step`             | | | | | Optional | | |
| `--study-dir`        | **Required** | **Required** | **Required** | **Required** | **Required** | **Required** |
| `--study-name`       | Optional | Optional | Optional | Optional | | Optional |
| `--subject`          | **Required** | **Required** | **Required** | **Required** | Optional | **Required** |
| `--surf-smooth`      | Optional | Optional | Optional | Optional | | |
| `--task`             | **Required** | **Required** | **Required** | **Required** | **Required** | **Required** |
| `--temp-dir`         | Optional | Optional | Optional | Optional | | |
| `--temp-json`        | | **Required** | **Required** | | | | 
| `--templates`        | Optional | Optional | Optional | Optional | | Optional |
| `--template1`        | Optional | Optional | Optional | Optional | | Optional |
| `--template2`        | Optional | Optional | Optional | Optional | | Optional |
| `--vol-smooth`       | Optional | Optional | Optional | Optional | | |
| `--wb-command`       | Usually **Required**<sup> 2</sup> | Usually **Required**<sup> 2</sup> | Usually **Required**<sup> 2</sup> | Usually **Required**<sup> 2</sup> | Usually **Required**<sup> 2</sup> | Usually **Required**<sup> 2</sup> |
| `--wrapper-location` | Rarely **Required<sup> 3</sup>** | **Required** | **Required** | Rarely **Required<sup> 3</sup>** | Rarely **Required<sup> 3</sup>** | Rarely **Required<sup> 3</sup>** |

#### Footnotes

<sup>1 </sup>Required unless running a specific kind of level 1 analysis that does not require event files (unless the level 1 `.fsf` file specifies that the basic waveform shape of each EV model is square; in other words, the level 1 `.fsf` file contains the line `set fmri(shape1) 0`).

<sup>2 </sup>Required unless the user already has an environment variable defining the location of this external software. For example, if the user has `WB_COMMAND=/home/users/shared/tools/workbench/bin/wb_command` in their `.bashrc` file, then the `--wb-command` argument may not be necessary.

<sup>3 </sup>Usually not required, unless running the script as a SLURM/SBATCH job.

</details>

<details><summary style="text-align: center;"><span style="display: inline; font-size: 1.17em;">Tree Diagram: Main Pipeline Script and Subscripts</span></summary>

The diagram below is just a visual explanation of the fact that `run_sbatch_jobs.py` runs many instances of `pipeline_wrapper.py`, which runs several instances of `run_level_1_analysis.py`, which runs 2 instances of `process_hemisphere.py`.

```
                                                      run_sbatch_jobs.py
                                                     /                  \
                                 pipeline_wrapper.py                     pipeline_wrapper.py
                                /                   \                    
         run_level_1_analysis.py                     run_level_1_analysis.py      
        /                       \                   /                       \
process_hemisphere.py   process_hemisphere.py   process_hemisphere.py   process_hemisphere.py
```

</details>

## Detailed Explanations of detect_outliers_and_QC.py Steps

This Python script will complete any of the following 6 steps to filter out invalid data from the event files:

1. Calculate and save summary statistics for each contrast for every subject.
2. Collect those summary statistics from each specific subject’s event `.tsv` files and save them all into one all-subject `.tsv` file per contrast.
3. Split the all-subject `.tsv` file into 2 more, 1 with all blank event `.tsv` files and 1 with non-blank `.tsv` files.
4. Split the non-blank all-subject `.tsv` file into 2 more, 1 with contrasts whose meaning depends on earlier (blank) contrasts, and 1 with contrasts whose meaning does not.
5. Split the non-blank-dependent (or just non-blank) all-subjects `.tsv` file into 2 more, 1 with all outliers (defined as 4 standard deviations from the mean) and 1 with the rest.
6. Make a new `.tsv` chart of which subjects have valid COPE (contrast parameter estimate) files for which contrast numbers.

## Metadata

Information about this `README` file:

- Created on 2021-10-01 by Greg Conan ([@GregConan](https://github.com/GregConan))
- Updated on 2021-10-01 by Greg Conan ([@GregConan](https://github.com/GregConan))
