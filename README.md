# anatprep — mri-space-atlas workflow

Builds a group MRI atlas from per-subject `desc-preproc_T2starw` images and
registers every subject into a target template space (default: ABAv3).

## Workflow overview

```
register_to_mean          (per subject/session)
    ↓
build_mri_atlas           (all subjects → mri-atlas.nii.gz)
    ↓
register_atlas_to_template
    ↓
compose_subject_to_template  (per subject/session)
```

| Rule | Tool | Output |
|---|---|---|
| `register_to_mean` | `antsRegistration` (Rigid+Affine+SyN) | per-subject affine, warp, warped image |
| `build_mri_atlas` | `AverageImages` | `mri-atlas.nii.gz` |
| `register_atlas_to_template` | `antsRegistration` (Rigid+Affine+SyN) | atlas→template transforms |
| `compose_subject_to_template` | `antsApplyTransforms` | composite `.h5` + warped image in template space |

## Requirements

- [ANTs](https://github.com/ANTsX/ANTs) ≥ 2.5 (installed via the `ants` conda env)
- [snakebids](https://github.com/akhanf/snakebids) ≥ 0.15
- Input BIDS dataset with `desc-preproc_T2starw.nii.gz` images:

```
{bids_dir}/sub-{subject}/ses-{session}/anat/sub-{subject}_ses-{session}_desc-preproc_T2starw.nii.gz
```

## Usage

### 1. Edit the config (optional)

Open `anatprep/config/snakebids.yml` and set defaults, or pass everything on
the command line.

### 2. Dry run

```bash
snakemake \
  --snakefile anatprep/workflow/Snakefile \
  --config \
    bids_dir=/path/to/bids \
    output_dir=/path/to/output \
    template_path=/path/to/initial_mean.nii.gz \
    template_path_target=/path/to/ABAv3_anat.nii.gz \
  --cores all \
  --use-conda \
  -n
```

### 3. Full run

Remove the `-n` flag:

```bash
snakemake \
  --snakefile anatprep/workflow/Snakefile \
  --config \
    bids_dir=/path/to/bids \
    output_dir=/path/to/output \
    template_path=/path/to/initial_mean.nii.gz \
    template_path_target=/path/to/ABAv3_anat.nii.gz \
  --cores all \
  --use-conda
```

### 4. Cluster (SLURM)

```bash
snakemake \
  --snakefile anatprep/workflow/Snakefile \
  --config \
    bids_dir=/path/to/bids \
    output_dir=/path/to/output \
    template_path=/path/to/initial_mean.nii.gz \
    template_path_target=/path/to/ABAv3_anat.nii.gz \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch --mem={resources.mem_mb}M --time={resources.runtime} --cpus-per-task={threads}" \
  --use-conda \
  --jobs 32
```

## Config keys

| Key | Required | Default | Description |
|---|---|---|---|
| `bids_dir` | yes | — | Path to the BIDS dataset |
| `output_dir` | yes | — | Root output directory |
| `template_path` | yes | — | Initial mean image for first-pass registration |
| `template_path_target` | yes | — | Target template NIfTI (e.g. ABAv3 anat) |
| `target_template` | no | `ABAv3` | Label used in output filenames |

## Output layout

```
{output_dir}/
  mri-atlas/
    mri-atlas.nii.gz
    from-mriatlas_to-ABAv3_0GenericAffine.mat
    from-mriatlas_to-ABAv3_1Warp.nii.gz
    from-mriatlas_to-ABAv3_1InverseWarp.nii.gz
  sub-{subject}/ses-{session}/
    xfm/
      sub-{subject}_ses-{session}_from-T2starw_to-mriatlas_0GenericAffine.mat
      sub-{subject}_ses-{session}_from-T2starw_to-mriatlas_1Warp.nii.gz
      sub-{subject}_ses-{session}_from-T2starw_to-mriatlas_1InverseWarp.nii.gz
      sub-{subject}_ses-{session}_from-T2starw_to-ABAv3_composite.h5
    anat/
      sub-{subject}_ses-{session}_space-ABAv3_desc-deformwarped_T2starw.nii.gz
```

## Selecting subjects

Use standard snakebids/pybids filters:

```bash
--config participant_label='["sub-01","sub-02"]'
```
