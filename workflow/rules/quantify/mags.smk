use rule preprocess__hosts as quantify__mags with:
    input:
        fa_gz=lambda w: features["mag_catalogues"][w.catalogue],
    output:
        fa_gz=QUANT_MAGS / "{catalogue}.fa.gz",
    log:
        QUANT_MAGS / "{catalogue}.log",
    conda:
        "../../environments/mags.yml"


rule quantify__mags__all:
    """Recompress all MAG catalogues"""
    input:
        [QUANT_MAGS / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],
