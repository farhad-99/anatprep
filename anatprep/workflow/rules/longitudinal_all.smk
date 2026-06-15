"""
Terminal rule: all_within_subject

Requests all subject template files and all coregistered session images
for both modalities across all subjects and sessions.
"""


def _all_templates():
    return [
        os.path.join(
            output_dir,
            f"sub-{sub}",
            "anat",
            f"sub-{sub}_{mod}_desc-subjectTemplate_anat.nii.gz",
        )
        for sub in subjects
        for mod in modalities
    ]


def _all_coregistered():
    return [
        os.path.join(
            output_dir,
            f"sub-{sub}",
            f"ses-{ses}",
            "anat",
            f"sub-{sub}_ses-{ses}_{mod}_space-subjectTemplate_desc-coregistered_anat.nii.gz",
        )
        for sub in subjects
        for ses in get_sessions(sub)
        for mod in modalities
    ]


rule all_within_subject:
    """Collect all subject templates and coregistered session images."""
    input:
        templates=_all_templates(),
        coregistered=_all_coregistered(),
