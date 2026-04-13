# anatprep Usage Guide

`anatprep` is an anatomical preprocessing [snakebids](https://snakebids.readthedocs.io) app designed for TRIDENT (15T) mouse MRI data. It takes a [BIDS](https://bids.neuroimaging.io)-organized dataset, applies bias field correction, optional multi-scan super-resolution averaging, and registers the MRI to a mouse brain template to produce a brain-extracted, preprocessed image.

---

## Table of Contents

- [Installation](#installation)
- [Basic Usage](#basic-usage)
- [Positional Arguments](#positional-arguments)
- [Flags and Options](#flags-and-options)
  - [App-Specific Flags](#app-specific-flags)
  - [Snakebids / BIDS App Standard Flags](#snakebids--bids-app-standard-flags)
- [Input Wildcards](#input-wildcards)
- [Input Suffix and Filters](#input-suffix-and-filters)
- [Templates](#templates)
  - [SPIM Templates (`--template`)](#spim-templates---template)
  - [MRI Brain-Masking Templates (`--template_mri`)](#mri-brain-masking-templates---template_mri)
  - [Atlas Segmentations (`--atlas_segs`)](#atlas-segmentations---atlas_segs)
- [Workflow Overview](#workflow-overview)
- [Output Files](#output-files)
- [Examples](#examples)

---

## Installation

`anatprep` is managed with [pixi](https://prefix.dev/docs/pixi/). Clone the repository and install:

```bash
git clone https://github.com/farhad-99/anatprep
cd anatprep
pixi install          # default (CPU) environment
# or
pixi install -e gpu   # GPU environment (requires CUDA 12.0+)
```

Run the app through pixi:

```bash
pixi run anatprep <bids_dir> <output_dir> <analysis_level> [options]
```

Or activate the environment and call it directly:

```bash
pixi shell
anatprep <bids_dir> <output_dir> <analysis_level> [options]
```

---

## Basic Usage

```
anatprep <bids_dir> <output_dir> participant [options]
```

| Argument         | Description                                          |
|------------------|------------------------------------------------------|
| `bids_dir`       | Path to the root of the input BIDS dataset           |
| `output_dir`     | Path where outputs will be written                   |
| `analysis_level` | Must be `participant` (the only supported level)     |

---

## Flags and Options

### App-Specific Flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--work_dir PATH` | `str` | `null` | Local path to use for temporary/intermediate files (e.g. `$SLURM_TMPDIR`). Useful on HPC systems to avoid writing temp files to slow network storage. |
| `--template TEMPLATE` | choice | `ABAv3` | Mouse brain template to use for registration. See [SPIM Templates](#spim-templates---template). |
| `--atlas_segs SEG [SEG ...]` | list | template default | One or more atlas segmentation labels to apply with the chosen template. Pass multiple values space-separated. See [Atlas Segmentations](#atlas-segmentations---atlas_segs). |
| `--template_crop {left,right}` | choice | `null` | Crop the template along the X-axis to retain only one hemisphere before registration. Useful for hemispheric analyses. |
| `--template_mri TEMPLATE` | choice | `MouseIn` | Template to use for MRI-to-template registration in order to propagate a brain mask. See [MRI Brain-Masking Templates](#mri-brain-masking-templates---template_mri). |
| `--mri_resample_percent PCT` | `int` | `100` | Resampling factor (as a percentage) applied when averaging multiple MRI scans for super-resolution. Use `200` to upsample to twice the native resolution. |
| `--sloppy` | flag | `False` | Use fast, low-quality registration parameters. **For testing only — do not use for production analyses.** |
| `--skip-bids-validation` | flag | `False` | Skip BIDS dataset validation. By default the dataset is validated with `bids-validator` (if installed) or `pybids`. |

### Snakebids / BIDS App Standard Flags

These are provided by snakebids and are available in all snakebids apps:

| Flag | Description |
|------|-------------|
| `--participant-label LABEL [LABEL ...]` | Process only the listed subject(s) (e.g. `--participant-label 001 002`). |
| `--exclude-participant-label LABEL [LABEL ...]` | Exclude specific subject(s) from processing. |
| `--derivatives` | Allow the app to search for inputs inside `derivatives/` subdirectories of the BIDS dataset. |
| `--pybidsdb-dir PATH` | Path for caching the pybids database (speeds up repeated runs on large datasets). |
| `--pybidsdb-reset` | Force regeneration of the pybids database cache. |
| `-n` / `--dry-run` | Snakemake dry-run: show what would be executed without running anything. |
| `-c N` / `--cores N` | Number of CPU cores to use. |

> **Tip:** Any extra arguments after `--` are passed directly to Snakemake (e.g. `-- --profile slurm`).

---

## Input Wildcards

`anatprep` discovers input files using [pybids](https://bids-standard.github.io/pybids/). The following BIDS entities are used as wildcards and are automatically extracted from matching filenames:

| Wildcard | BIDS Entity | Example value |
|----------|-------------|---------------|
| `subject` | `sub-` | `001` |
| `session` | `ses-` | `baseline` |
| `acquisition` | `acq-` | `highres` |
| `run` | `run-` | `01` |
| `reconstruction` | `rec-` | `avgecho` |
| `suffix` | file suffix | `T2w` |
| `extension` | file extension | `nii.gz` |

These wildcards are resolved at runtime and propagate through all output file names, keeping outputs organized and traceable back to their inputs.

---

## Input Suffix and Filters

By default `anatprep` looks for files matching:

```
sub-*/[ses-*/]anat/sub-*[_ses-*][_acq-*][_run-*][_rec-*]_T2w.nii.gz
```

The filters applied are:

| Filter | Value |
|--------|-------|
| `suffix` | `T2w` |
| `extension` | `nii.gz` |
| `datatype` | `anat` |

> **Note:** If your dataset uses a different suffix (e.g. `T2starw`) you can override the input filter via `--filter-mri suffix=T2starw`. If more than one suffix is detected across the dataset, the app will raise an error and ask you to narrow the selection with `--filter-mri`.

---

## Templates

### SPIM Templates (`--template`)

These templates are used as registration targets. All are automatically downloaded from Zenodo or the mouse imaging repository on first use.

| Value | Description |
|-------|-------------|
| `ABAv3` *(default)* | Allen Brain Atlas v3. Available atlas segmentations: `all`, `allsym`, `coarse`, `mid`, `fine`, `roi22`, `roi82`, `roi198`. Default segmentations: `all`, `coarse`, `mid`, `fine`. |
| `DSURQE` | Dorr-Steadman-Ullmann-Richards-Qiu-Egan atlas (40 µm). Atlas: `all`. |
| `gubra` | Gubra mouse brain template. |
| `MBMv3` | McGill Brain Minerva v3. Atlas: `paxinos`. |
| `turone` | Turone mouse brain template. Atlas: `all`. |

### MRI Brain-Masking Templates (`--template_mri`)

A second template is used specifically to propagate a brain mask onto the subject MRI.

| Value | Description |
|-------|-------------|
| `MouseIn` *(default)* | MouseIn template (MP2RAGE inv-1 contrast). |
| `DSURQE` | DSURQE 40 µm template. |
| `MBMv3` | McGill Brain Minerva v3. |
| `turone` | Turone mouse brain template. |

### Atlas Segmentations (`--atlas_segs`)

Select one or more segmentation granularities from those defined for the chosen template. Example for `ABAv3`:

| Value | Description |
|-------|-------------|
| `all` | Full parcellation (all regions) |
| `allsym` | Symmetric full parcellation |
| `coarse` | Coarse-grained regions |
| `mid` | Mid-level parcellation |
| `fine` | Fine-grained parcellation |
| `roi22` | 22-ROI parcellation |
| `roi82` | 82-ROI parcellation |
| `roi198` | 198-ROI parcellation |

Pass multiple values:

```bash
anatprep /data/bids /data/out participant \
  --template ABAv3 \
  --atlas_segs coarse mid fine
```

---

## Workflow Overview

The workflow runs the following stages in order:

1. **N4 Bias Field Correction** — Each individual MRI is corrected for intensity inhomogeneity using ANTs `N4BiasFieldCorrection`.

2. **Resampling** — Each bias-corrected image is resampled at the percentage specified by `--mri_resample_percent` (default 100 %; use 200 for 2× super-resolution).

3. **Multi-MRI Rigid Registration & Averaging** *(when multiple scans exist for the same subject)* — All scans are rigidly aligned to the first scan and averaged to produce a super-resolved image (`desc-preproc`).

4. **Rigid + Deformable MRI → MRI Template Registration** — The preprocessed MRI is registered to the chosen `--template_mri` brain template using greedy (rigid 6 DOF followed by deformable) so that the template brain mask can be transferred to subject space.

5. **Brain Mask Propagation** — The template brain mask is warped back to subject MRI space using the inverse transformations computed in step 4.

6. **Brain Extraction** — The propagated mask is applied to the preprocessed MRI to yield a brain-extracted image (`desc-brain`).

---

## Output Files

All outputs follow BIDS derivative naming conventions and are written under `<output_dir>/`:

| File | Description |
|------|-------------|
| `sub-*/[ses-*/]anat/sub-*_desc-preproc_T2w.nii.gz` | Bias-field-corrected, optionally super-resolved (averaged) MRI. |
| `sub-*/[ses-*/]anat/sub-*_desc-brain_T2w.nii.gz` | Brain-extracted MRI (brain mask applied). |

The `T2w` suffix in the output filenames reflects the detected input suffix and will match whatever suffix was found in your dataset (e.g. `T2starw`).

---

## Examples

### Minimal run (single subject)

```bash
anatprep /data/bids /data/out participant \
  --participant-label 001
```

### Using a different brain-masking template

```bash
anatprep /data/bids /data/out participant \
  --template_mri DSURQE
```

### Super-resolution averaging of multiple scans

```bash
anatprep /data/bids /data/out participant \
  --mri_resample_percent 200
```

### Select specific atlas segmentations

```bash
anatprep /data/bids /data/out participant \
  --template ABAv3 \
  --atlas_segs coarse fine roi82
```

### Hemispheric analysis (right hemisphere only)

```bash
anatprep /data/bids /data/out participant \
  --template ABAv3 \
  --template_crop right
```

### Dry run to check what will execute

```bash
anatprep /data/bids /data/out participant \
  --participant-label 001 \
  -- -n
```

### HPC run with SLURM temporary directory

```bash
anatprep /data/bids /data/out participant \
  --work_dir $SLURM_TMPDIR \
  -- --profile slurm --jobs 50
```

### Quick test with low-quality parameters

```bash
anatprep /data/bids /data/out participant \
  --sloppy \
  --participant-label 001
```
