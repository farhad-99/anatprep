"""
Rule: trim_and_center

Per session, per modality:
  - Load image with nibabel
  - Clip intensities at 10th–90th percentile of foreground voxels
  - Reset the image origin to the intensity-weighted centre of mass in world coordinates
  - Save the result
"""


rule trim_and_center:
    """Clip hot/cold voxels and re-centre the image origin to intensity-weighted COM."""
    input:
        nii=os.path.join(
            bids_dir,
            "sub-{subject}",
            "ses-{session}",
            "anat",
            "sub-{subject}_ses-{session}_{modality}.nii.gz",
        ),
    output:
        nii=os.path.join(
            output_dir,
            "sub-{subject}",
            "ses-{session}",
            "anat",
            "sub-{subject}_ses-{session}_{modality}_desc-trimmed_anat.nii.gz",
        ),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=15,
    log:
        os.path.join("logs", "trim_and_center", "{subject}_{session}_{modality}.log"),
    run:
        import logging
        import nibabel as nib
        import numpy as np

        logging.basicConfig(filename=log[0], level=logging.INFO)
        log_msg = logging.info

        log_msg(f"Loading {input.nii}")
        img = nib.load(input.nii)
        data = img.get_fdata(dtype=np.float32)
        affine = img.affine.copy()

        fg = data > 0
        if fg.sum() == 0:
            raise ValueError(f"No foreground voxels found in {input.nii}")

        fg_vals = data[fg]
        p10, p90 = np.percentile(fg_vals, [10, 90])
        log_msg(f"Clipping intensities to [{p10:.2f}, {p90:.2f}]")
        data = np.clip(data, p10, p90)
        # re-zero background after clipping
        data[~fg] = 0.0

        # intensity-weighted centre of mass in voxel space then world space
        weights = data.copy()
        weights[~fg] = 0.0
        total_weight = weights.sum()
        coords = np.indices(data.shape).astype(np.float32)
        com_vox = np.array(
            [(coords[i] * weights).sum() / total_weight for i in range(3)]
        )
        com_world = affine[:3, :3] @ com_vox + affine[:3, 3]
        log_msg(f"COM world coords: {com_world}")

        affine[:3, 3] = -com_world
        os.makedirs(os.path.dirname(output.nii), exist_ok=True)
        out_img = nib.Nifti1Image(data, affine, img.header)
        nib.save(out_img, output.nii)
        log_msg(f"Saved trimmed image to {output.nii}")
