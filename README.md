# anatprep

Anatomical preprocessing snakebids app for TRIDENT (15T) mouse data.

Registers in-vivo MRI (T2starw) to a standard template using ANTs (rigid + affine + SyN), propagates a brain mask back to subject space, and outputs skull-stripped MRI ready for downstream analysis.

---

## Requirements

- [Pixi](https://pixi.sh) (manages all dependencies via `pixi.lock`)
- Linux x86-64

Install the environment:

```bash
pixi install
```

---

## Input data

Your dataset must follow [BIDS](https://bids-specification.readthedocs.io). The pipeline expects T2starw anatomical images:

```
bids_dataset/
├── dataset_description.json
├── sub-01/
│   └── anat/
│       └── sub-01_T2starw.nii.gz
├── sub-02/
│   └── ses-01/
│       └── anat/
│           └── sub-02_ses-01_T2starw.nii.gz
└── ...
```

---

## Basic usage

```bash
pixi run anatprep <bids_dir> <output_dir> participant [options]
```

| Argument | Description |
|---|---|
| `bids_dir` | Path to the input BIDS dataset |
| `output_dir` | Path for pipeline outputs |
| `participant` | Analysis level (only `participant` is supported) |

### Minimal example

```bash
pixi run anatprep /data/myproject/bids /data/myproject/anatprep participant \
    --use-conda --cores all
```

> **`--use-conda` is required.** Each rule runs inside its own conda environment
> (`ants`, `c3d`). Without this flag Snakemake ignores those environments and
> tools like `N4BiasFieldCorrection` will not be found.

### Specify subjects or sessions

```bash
# Single subject
pixi run anatprep /data/bids /data/out participant \
    --use-conda --participant-label 01 --cores all

# Multiple subjects
pixi run anatprep /data/bids /data/out participant \
    --use-conda --participant-label 01 02 03 --cores all
```

### Multi-echo datasets

If your dataset contains individual echo images (`echo-1`, `echo-2`, …) alongside a reconstructed magnitude image (`rec-*`), the pipeline uses the reconstructed files by default (`echo: false` filter in the config). If you need to override this and select a specific echo instead, pass `--filter_mri` at runtime:

```bash
# Use only echo-1 images
pixi run anatprep /data/bids /data/out participant \
    --use-conda --filter_mri echo=1 --cores all

# Use only a specific reconstruction
pixi run anatprep /data/bids /data/out participant \
    --use-conda --filter_mri reconstruction=preproc --cores all
```

---

## Key options

### Template selection

`--template` chooses the SPIM reference space (default: `ABAv3`):

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --template ABAv3 --cores all
```

Available templates: `ABAv3`, `DSURQE`, `gubra`, `MBMv3`, `turone`

### MRI brain-masking template

`--template_mri` selects the MRI template used to propagate the brain mask (default: `MouseIn`):

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --template_mri MouseIn --cores all
```

Available: `MouseIn`, `DSURQE`, `MBMv3`, `turone`

### Multi-MRI super-resolution averaging

If multiple MRI runs exist per subject, they are rigidly aligned and averaged. Use `--mri_resample_percent` to upsample before averaging (e.g. `200` = 2× resolution):

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --mri_resample_percent 200 --cores all
```

### Hemisphere cropping

Crop the template along the X-axis to a single hemisphere:

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --template_crop left --cores all
```

### Working directory

Redirect temporary files to a fast local disk:

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --work_dir /scratch/tmp --cores all
```

### Skip BIDS validation

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --skip-bids-validation --cores all
```

---

## Cluster / HPC execution

The pipeline uses Snakemake and supports SLURM out of the box.

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda \
    --executor slurm \
    --default-resources slurm_account=myaccount mem_mb=8000 runtime=60 \
    --cores all \
    --jobs 50
```

---

## Dry run

Preview the jobs Snakemake would execute without running anything:

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --cores all --dry-run
```

---

## Testing with low-quality parameters

Use `--sloppy` to run faster with reduced registration quality (for pipeline testing only):

```bash
pixi run anatprep /data/bids /data/out participant \
    --use-conda --sloppy --cores all
```

---

## Outputs

All outputs are written to `<output_dir>` in BIDS derivative format:

| File | Description |
|---|---|
| `sub-*/anat/sub-*_desc-preproc_T2starw.nii.gz` | N4-corrected, motion-averaged MRI |
| `sub-*/anat/sub-*_desc-brain_mask.nii.gz` | Brain mask in subject space |
| `sub-*/anat/sub-*_desc-brain_T2starw.nii.gz` | Skull-stripped MRI |
| `sub-*/xfm/sub-*_from-T2starw_to-<template>_*0GenericAffine.mat` | ANTs affine transform (subject → template) |
| `sub-*/xfm/sub-*_from-T2starw_to-<template>_*1Warp.nii.gz` | ANTs SyN warp (subject → template) |
| `sub-*/xfm/sub-*_from-T2starw_to-<template>_*1InverseWarp.nii.gz` | ANTs SyN inverse warp (template → subject) |

---

## Development

```bash
# Install dev environment
pixi install -e dev

# Format and lint
pixi run quality_fix

# Check formatting without modifying
pixi run quality_check
```
