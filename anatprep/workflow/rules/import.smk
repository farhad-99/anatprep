"""Import template anatomical image and brain mask for MRI registration."""


localrules:
    import_template_anat,
    import_mask,


rule import_template_anat:
    """Import template anatomical image."""
    input:
        anat=lambda wildcards: storage(
            ancient(resources_path(config["templates"][wildcards.template]["anat"]))
        ),
    output:
        anat=bids(root=root, template="{template}", suffix="anat.nii.gz"),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    log:
        bids(
            root="logs",
            datatype="import_anat",
            template="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/copy_nii.py"


rule import_mask:
    input:
        mask=lambda wildcards: storage(
            ancient(resources_path(config["templates"][wildcards.template]["mask"]))
        ),
    output:
        mask=bids(root=root, template="{template}", desc="brain", suffix="mask.nii.gz"),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    log:
        bids(
            root="logs",
            datatype="import_mask",
            template="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/copy_nii.py"
