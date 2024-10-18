include: "__functions__.smk"
include: "coverm.smk"


rule quantify:
    """Quantify MAG abundances"""
    input:
        rules.quantify__coverm.input,
