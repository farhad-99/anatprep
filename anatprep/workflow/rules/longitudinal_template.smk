"""
Rule: build_subject_template

Per subject, per modality — aggregates all sessions:
  - Collect all trimmed session images for that subject and modality
  - Run ants.build_template() to produce a bias-corrected unbiased midpoint
    average across sessions
  - Save the subject template
"""


def get_trimmed_sessions(wildcards):
    """Return all trimmed images for a given subject and modality."""
    sessions = get_sessions(wildcards.subject)
    return [
        os.path.join(
            output_dir,
            f"sub-{wildcards.subject}",
            f"ses-{ses}",
            "anat",
            f"sub-{wildcards.subject}_ses-{ses}_{wildcards.modality}_desc-trimmed_anat.nii.gz",
        )
        for ses in sessions
    ]


rule build_subject_template:
    """Build an unbiased subject template from all trimmed longitudinal sessions."""
    input:
        trimmed=get_trimmed_sessions,
    output:
        template=os.path.join(
            output_dir,
            "sub-{subject}",
            "anat",
            "sub-{subject}_{modality}_desc-subjectTemplate_anat.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        os.path.join("logs", "build_subject_template", "{subject}_{modality}.log"),
    run:
        import logging
        import shutil
        import tempfile
        import ants

        logging.basicConfig(filename=log[0], level=logging.INFO)
        log_msg = logging.info

        log_msg(f"Building template for sub-{wildcards.subject} {wildcards.modality}")
        log_msg(f"Input images: {input.trimmed}")

        image_list = [ants.image_read(f) for f in input.trimmed]

        with tempfile.TemporaryDirectory() as tmpdir:
            template = ants.build_template(
                image_list=image_list,
                iterations=3,
                gradient_step=0.2,
                blending_weight=0.75,
                verbose=True,
                outprefix=os.path.join(tmpdir, "template_"),
            )

        os.makedirs(os.path.dirname(output.template), exist_ok=True)
        ants.image_write(template, output.template)
        log_msg(f"Saved subject template to {output.template}")
