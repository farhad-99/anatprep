"""
MRI preprocessing and cross-modality registration workflow for SPIMquant.

This module handles co-registration of in-vivo MRI (T2w, or other) with ex-vivo SPIM data,
enabling multi-modal analysis and assessment of tissue changes due to perfusion
fixation and optical clearing.

Key workflow stages:
1. N4 bias field correction of MRI
2. MRI to MRI-template rigid+nlin registration (for brain masking)
3. Template brain mask to MRI transformation
4. MRI brain extraction
5. MRI to SPIM affine+nlin registration
6. Parameter tuning rules for optimization
7. Concatenated transformations (MRI -> SPIM -> Template)
8. Jacobian determinant calculation (tissue deformation quantification)
9. Registration QC report generation

This workflow is optional and used when both MRI and SPIM data are available
for the same subject. 
"""


def get_all_mri(wildcards):
    """Get all MRI files for a subject (for multi-MRI averaging)."""
    return sorted(
        inputs["mri"]
        .filter(**wildcards) #subject=wildcards.subject, session=wildcards.session)
        .expand(bids(root=root, datatype="anat", desc="N4", **inputs["mri"].wildcards))
    )


def get_ref_mri(wildcards):
    """gets a N4resampled reference"""
    return sorted(
        inputs["mri"]
        .filter(**wildcards)
        .expand(
            bids(
                root=root,
                datatype="anat",
                desc="N4resampled",
                **inputs["mri"].wildcards,
            )
        )
    )[0]


def get_flo_mri(wildcards):
    return get_all_mri(wildcards)[1:]


def get_mri_indices(wildcards):
    """Get list of indices for MRI files for a subject."""
    return list(range(len(get_all_mri(wildcards))))


def get_mri_by_index(wildcards):
    return get_all_mri(wildcards)[int(wildcards.mrindex)]


rule n4_mri_individual:
    """Apply N4 bias field correction to individual T2w MRI images.
    
    This rule only runs when multiple MRI images are present for a subject.
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities in
    each MRI image before registration and averaging.
    """
    input:
        nii=inputs["mri"].path,
    output:
        nii=temp(bids(root=root, datatype="anat", desc="N4", **inputs["mri"].wildcards)),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "N4BiasFieldCorrection -i {input.nii}"
        " -o {output.nii}"
        " -d 3 -v "


rule resample_mri_ref:
    """apply resampling by mri_resample_percent"""
    input:
        nii=bids(root=root, datatype="anat", desc="N4", **inputs["mri"].wildcards),
    params:
        resample_percent=config["mri_resample_percent"],
    output:
        nii=temp(
            bids(
                root=root,
                datatype="anat",
                desc="N4resampled",
                **inputs["mri"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input} -interpolation Cubic -resample {params.resample_percent}% -o {output}"


rule register_mri_to_first:
    """Rigidly register each MRI to the first MRI using greedy.
    
    This aligns multiple MRI acquisitions to enable super-resolved averaging.
    Only runs when multiple MRI images are present.
    """
    input:
        fixed=get_ref_mri,
        moving=get_mri_by_index,
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=f"{mri_suffix}",
                to=f"{mri_suffix}ref",
                type_="ras",
                desc="rigid",
                mrindex="{mrindex}",
                suffix="xfm.txt",
                **inputs.subj_wildcards,
            )
        ),
    threads: workflow.cores
    resources:
        mem_mb=8000,
        runtime=10,
    conda:
        "../envs/c3d.yaml"
    shell:
        "if [ {wildcards.mrindex} -eq 0 ]; then "
        "  echo '1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1' > {output.xfm_ras}; "
        "else "
        "  greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        "  -a -dof 6 -ia-image-centers -m SSD -o {output.xfm_ras}; "
        "fi"


rule resample_mri_to_first:
    """Apply transformation to resample mri to first"""
    input:
        fixed=get_ref_mri,
        moving=get_mri_by_index,
        xfm_ras=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=f"{mri_suffix}ref",
            type_="ras",
            desc="rigid",
            mrindex="{mrindex}",
            suffix="xfm.txt",
            **inputs.subj_wildcards,
        ),
    output:
        resampled=temp(
            bids(
                root=root,
                datatype="anat",
                desc="resampled",
                mrindex="{mrindex}",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            )
        ),
    threads: workflow.cores
    conda:
        "../envs/c3d.yaml"
    resources:
        mem_mb=8000,
        runtime=10,
    shell:
        "c3d -interpolation Cubic {input.fixed} {input.moving} -reslice-matrix {input.xfm_ras} -o {output.resampled}"


rule average_mri:
    """Average all registered and N4-corrected MRI images.
    
    Optionally upsamples images during transformation for super-resolution.
    Produces a final preprocessed MRI with desc-preproc.
    """
    input:
        resampled_images=lambda wildcards: expand(
            bids(
                root=root,
                datatype="anat",
                desc="resampled",
                mrindex="{mrindex}",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            ),
            mrindex=get_mri_indices(wildcards),
            allow_missing=True,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    threads: workflow.cores
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.resampled_images} -accum -add -endaccum -o {output.nii}"


# 


rule rigid_nlin_reg_mri_to_template:
    """Initial unmasked MRI to unmasked template MRI using ANTs registration.

    Performs rigid, affine, and SyN deformable registration to align
    non-brain-masked MRI with the template so that the template brain
    mask can be propagated to native space.
    """
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    params:
        prefix=lambda wildcards, output: output.affine[: -len("0GenericAffine.mat")],
    output:
        affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="{template}",
            suffix="0GenericAffine.mat",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="{template}",
            suffix="1Warp.nii.gz",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="{template}",
            suffix="1InverseWarp.nii.gz",
            **inputs.subj_wildcards,
        ),
        warped=bids(
            root=root,
            datatype="xfm",
            space="{template}",
            desc="deformwarped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    threads: workflow.cores
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/ants.yaml"
    shell:
        "antsRegistration"
        " --dimensionality 3"
        " --float 0"
        ' --output ["{params.prefix}","{output.warped}"]'
        " --interpolation Linear"
        " --use-histogram-matching 0"
        " --winsorize-image-intensities [0.005,0.995]"
        ' --initial-moving-transform ["{input.template}","{input.subject}",1]'
        " --transform Rigid[0.1]"
        ' --metric MI["{input.template}","{input.subject}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform Affine[0.1]"
        ' --metric MI["{input.template}","{input.subject}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform SyN[0.1,3,0]"
        ' --metric CC["{input.template}","{input.subject}",1,4]'
        " --convergence [100x70x50x20,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " -v 1"


rule all_tune_mri_mask:
    input:
        inputs["mri"].expand(
            bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix="mask.nii.gz",
                **inputs.subj_wildcards,
            ),
        ),


rule transform_template_mask_to_mri:
    input:
        mask=bids(
            root=root,
            template=config["template_mri"],
            desc="brain",
            suffix="mask.nii.gz",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=config["template_mri"],
            suffix="0GenericAffine.mat",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=config["template_mri"],
            suffix="1InverseWarp.nii.gz",
            **inputs.subj_wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs.subj_wildcards,
        ),
    threads: workflow.cores
    resources:
        mem_mb=4000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "antsApplyTransforms"
        " -d 3"
        " -i {input.mask}"
        " -r {input.ref}"
        " -o {output.mask}"
        " -n NearestNeighbor"
        " -t [{input.affine},1]"
        " -t {input.invwarp}"


rule apply_mri_brain_mask:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.nii} {input.mask} -multiply -o {output.nii}"


