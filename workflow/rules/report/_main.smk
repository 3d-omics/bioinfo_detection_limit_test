include: "_functions.smk"
include: "step.smk"
include: "library.smk"


rule report:
    input:
        rules.report__step.input,
        rules.report__library.input,
