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
    threads: 8
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
    threads: 8
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
    threads: 8
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.resampled_images} -accum -add -endaccum -o {output.nii}"


rule rigid_nlin_reg_mri_to_template:
    """Initial unmasked MRI to unmasked template MRI using rigid + 
    deformable registration.
    
    Performs initial rigid (6 DOF) alignment followed by deformable
    registration to align non-brain-masked  MRI with similarly
    non-brain-masked MRI, so that the template brain mask can be
    propagated. Note: this step could be replaced by a ML-based 
    brain-masking if a suitable model is found/trained.
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
        iters="{iters}",  #"100x50x50",
        metric="NCC {radius}",  #3x3x3",
        sigma1="{gradsigma}vox",  #3
        sigma2="{warpsigma}vox",  #2
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=f"{mri_suffix}",
                to="{template}",
                type_="ras",
                desc="rigid",
                suffix="xfm.txt",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        warp=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=f"{mri_suffix}",
                to="{template}",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="xfm",
                from_="{template}",
                to=f"{mri_suffix}",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="xfm",
                space="{template}",
                desc="deformwarped",
                suffix=f"{mri_suffix}.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
    threads: 32
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -a -dof 6 -ia-image-centers -m {params.metric} -o {output.xfm_ras} && "
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -it {output.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.warp} {output.xfm_ras} "


rule all_tune_mri_mask:
    input:
        inputs["mri"].expand(
            bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix="mask.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            ),
            iters="100x50x50",
            gradsigma=range(3, 6),
            warpsigma=range(3, 6),
            radius=[f"{i}x{i}x{i}" for i in range(2, 5)],
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
        xfm_ras=bids(
            root=root,
            datatype="xfm",
            from_=f"{mri_suffix}",
            to=config["template_mri"],
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_=config["template_mri"],
            to=f"{mri_suffix}",
            suffix="warp.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs.subj_wildcards,
        ),
    output:
        mask=temp(
            bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix="mask.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
    shadow:
        "minimal"
    threads: 32
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        " c3d_affine_tool {input.xfm_ras} -inv -o inv_rigid.txt && "
        " greedy -threads {threads} -d 3 -ri NN -rf {input.ref} "
        "  -rm {input.mask} {output.mask} "
        "  -r inv_rigid.txt {input.invwarp}"


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
            iters=config["reg_mri"]["rigidgreedy"]["iters"],
            radius=config["reg_mri"]["rigidgreedy"]["radius"],
            gradsigma=config["reg_mri"]["rigidgreedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["rigidgreedy"]["warpsigma"],
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


