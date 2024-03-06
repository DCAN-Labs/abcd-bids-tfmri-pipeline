# abcd_extract_eprime-master

## Purpose

Matlab scripts abcd_extract_eprime_* process the e-prime timing files for the ABCD Study.
The outputs of these scripts are the task timing ready for the AFNI's deconvolution and also behavioral metrics.

## Example

An example e-prime file (`NDAR_INV********_WM.txt`) will be included in future updates, with all of its personally identifying data de-identified and replaced with arbitrary or randomly generated values. It will be included only to show the format/structure of the e-prime file(s) for get_events_extract_eprime.py to read.

Below is an example MATLAB command to process the e-prime file for that (nonexistent) subject:

```matlab
%% MMPS Auxiliary functions are called, so  add 'aux' and 'mmil_utils' directories to the MATLAB path.
addpath mmil_utils
addpath aux

%%
abcd_extract_eprime_nback('NDAR_INV********_WM.txt')
```

