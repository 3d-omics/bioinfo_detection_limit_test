include: "functions.smk"
include: "step.smk"
include: "library.smk"


rule report:
    input:
        rules.report_step.input,
        rules.report_library.input,
