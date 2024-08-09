include: "__functions__.smk"
include: "step.smk"
include: "library.smk"


rule report:
    """Generate all the per step and per library reports"""
    input:
        rules.report__step.input,
        rules.report__library.input,
