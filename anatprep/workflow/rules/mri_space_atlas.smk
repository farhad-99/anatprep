"""
MRI space atlas construction and subject-to-template registration.

Workflow:
1. compute_initial_mean    — Average all desc-preproc images into a first-pass mean.
2. register_to_mean        — Register each subject/session T2starw to the initial mean.
3. build_mri_atlas         — Average warped images into an mri-atlas template.
4. register_atlas_to_template — Register mri-atlas to the target template (default: ABAv3).
5. compose_subject_to_template — Compose subject→mri-atlas and mri-atlas→template transforms.
"""


ruleorder: register_to_mean > rigid_nlin_reg_mri_to_template


def get_all_preproc(wildcards=None):
    """Get all desc-preproc T2starw images across subjects/sessions."""
    return inputs["mri"].expand(
        bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        allow_missing=False,
    )


def get_all_warped(wildcards=None):
    """Get all per-subject/session warped images for atlas construction."""
    return inputs["mri"].expand(
        bids(
            root=root,
            datatype="anat",
            space="mriatlas",
            desc="warped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        allow_missing=False,
    )


rule compute_initial_mean:
    """Average all desc-preproc T2starw images to create the first-pass mean target."""
    input:
        images=get_all_preproc,
    output:
        mean=os.path.join(root, "mri-atlas", "initial_mean.nii.gz"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/ants.yaml"
    shell:
        "mkdir -p $(dirname {output.mean}) && "
        "AverageImages 3 {output.mean} 0 {input.images}"


rule register_to_mean:
    """Register each subject/session T2starw to the computed initial mean (first pass).

    Produces per-subject affine and warp transforms, plus the warped image in
    mri-atlas space.
    """
    input:
        fixed=os.path.join(root, "mri-atlas", "initial_mean.nii.gz"),
        moving=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    params:
        prefix=lambda wildcards, output: output.affine.removesuffix("0GenericAffine.mat"),
    output:
        affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="0GenericAffine.mat",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1Warp.nii.gz",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1InverseWarp.nii.gz",
            **inputs.subj_wildcards,
        ),
        warped=bids(
            root=root,
            datatype="anat",
            space="mriatlas",
            desc="warped",
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
        ' --initial-moving-transform ["{input.fixed}","{input.moving}",1]'
        " --transform Rigid[0.1]"
        ' --metric MI["{input.fixed}","{input.moving}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform Affine[0.1]"
        ' --metric MI["{input.fixed}","{input.moving}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform SyN[0.1,3,0]"
        ' --metric CC["{input.fixed}","{input.moving}",1,4]'
        " --convergence [100x70x50x20,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --number-of-threads {threads}"
        " -v 1"


rule build_mri_atlas:
    """Average all warped images into an mri-atlas template.

    Uses ANTs AverageImages to produce a mean image across all warped
    subject/session images.
    """
    input:
        warped=get_all_warped,
    output:
        atlas=os.path.join(root, "mri-atlas", "mri-atlas.nii.gz"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    conda:
        "../envs/ants.yaml"
    shell:
        "AverageImages 3 {output.atlas} 1 {input.warped}"


rule register_atlas_to_template:
    """Register mri-atlas to the target template using ANTs.

    Produces affine and SyN warp transforms from mri-atlas to the target
    template space (default: ABAv3).
    """
    input:
        fixed=config["template_path_target"],
        moving=os.path.join(root, "mri-atlas", "mri-atlas.nii.gz"),
    params:
        prefix=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_"),
    output:
        affine=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_0GenericAffine.mat"),
        warp=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_1Warp.nii.gz"),
        invwarp=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_1InverseWarp.nii.gz"),
        warped=temp(os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_warped.nii.gz")),
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
        ' --initial-moving-transform ["{input.fixed}","{input.moving}",1]'
        " --transform Rigid[0.1]"
        ' --metric MI["{input.fixed}","{input.moving}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform Affine[0.1]"
        ' --metric MI["{input.fixed}","{input.moving}",1,32,Regular,0.25]'
        " --convergence [1000x500x250x100,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --transform SyN[0.1,3,0]"
        ' --metric CC["{input.fixed}","{input.moving}",1,4]'
        " --convergence [100x70x50x20,1e-6,10]"
        " --shrink-factors 8x4x2x1"
        " --smoothing-sigmas 3x2x1x0vox"
        " --number-of-threads {threads}"
        " -v 1"


rule compose_subject_to_template:
    """Compose subject→mri-atlas and mri-atlas→template transforms.

    Concatenates per-subject transforms with the atlas→template transforms
    using antsApplyTransforms to produce a direct subject→template composite
    warp and the subject image warped into template space.
    """
    input:
        moving=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        fixed=config["template_path_target"],
        sub_affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="0GenericAffine.mat",
            **inputs.subj_wildcards,
        ),
        sub_warp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1Warp.nii.gz",
            **inputs.subj_wildcards,
        ),
        atlas_affine=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_0GenericAffine.mat"),
        atlas_warp=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_1Warp.nii.gz"),
    output:
        composite=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=target_template,
            suffix="composite.nii.gz",
            **inputs.subj_wildcards,
        ),
        warped=bids(
            root=root,
            datatype="anat",
            space=target_template,
            desc="deformwarped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    threads: workflow.cores
    resources:
        mem_mb=8000,
        runtime=30,
    conda:
        "../envs/ants.yaml"
    shell:
        "antsApplyTransforms"
        " -d 3"
        " -r {input.fixed}"
        " -o [{output.composite},1]"
        " -t {input.atlas_warp}"
        " -t {input.atlas_affine}"
        " -t {input.sub_warp}"
        " -t {input.sub_affine}"
        " --number-of-threads {threads}"
        " -v 1"
        " && "
        "antsApplyTransforms"
        " -d 3"
        " -i {input.moving}"
        " -r {input.fixed}"
        " -o {output.warped}"
        " -t {output.composite}"
        " --interpolation Linear"
        " --number-of-threads {threads}"
        " -v 1"
