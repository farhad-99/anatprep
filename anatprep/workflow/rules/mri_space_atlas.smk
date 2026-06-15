"""
MRI space atlas construction and subject-to-template registration.

Workflow:
1. register_to_mean        — Register each subject/session T2starw to an initial mean image.
2. build_mri_atlas         — Average warped images into an mri-atlas template.
3. register_atlas_to_template — Register mri-atlas to the target template (default: ABAv3).
4. compose_subject_to_template — Compose subject→mri-atlas and mri-atlas→template transforms.
"""


def get_all_warped(wildcards):
    """Get all per-subject warped images for atlas construction."""
    return inputs["mri"].expand(
        bids(
            root=root,
            datatype="anat",
            space="mriatlas",
            desc="warped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        )
    )


rule register_to_mean:
    """Register each subject/session T2starw to an initial mean image (first pass).

    Produces per-subject affine and warp transforms, plus the warped image in
    mri-atlas space.  The initial mean is supplied via config['template_path'].
    """
    input:
        fixed=config["template_path"],
        moving=inputs["mri"].path,
    params:
        prefix=lambda wildcards, output: output.affine.removesuffix("0GenericAffine.mat"),
    output:
        affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="0GenericAffine.mat",
            **inputs["mri"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1Warp.nii.gz",
            **inputs["mri"].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1InverseWarp.nii.gz",
            **inputs["mri"].wildcards,
        ),
        warped=bids(
            root=root,
            datatype="anat",
            space="mriatlas",
            desc="warped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs["mri"].wildcards,
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
        moving=inputs["mri"].path,
        fixed=config["template_path_target"],
        sub_affine=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="0GenericAffine.mat",
            **inputs["mri"].wildcards,
        ),
        sub_warp=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to="mriatlas",
            suffix="1Warp.nii.gz",
            **inputs["mri"].wildcards,
        ),
        atlas_affine=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_0GenericAffine.mat"),
        atlas_warp=os.path.join(root, "mri-atlas", f"from-mriatlas_to-{target_template}_1Warp.nii.gz"),
    output:
        composite=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=target_template,
            suffix="composite.h5",
            **inputs["mri"].wildcards,
        ),
        warped=bids(
            root=root,
            datatype="anat",
            space=target_template,
            desc="deformwarped",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs["mri"].wildcards,
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
        " -i {input.moving}"
        " -r {input.fixed}"
        " -o Linear[{output.composite}]"
        " -t {input.atlas_warp}"
        " -t {input.atlas_affine}"
        " -t {input.sub_warp}"
        " -t {input.sub_affine}"
        " -v 1"
        " && "
        "antsApplyTransforms"
        " -d 3"
        " -i {input.moving}"
        " -r {input.fixed}"
        " -o {output.warped}"
        " -t {output.composite}"
        " --interpolation Linear"
        " -v 1"
