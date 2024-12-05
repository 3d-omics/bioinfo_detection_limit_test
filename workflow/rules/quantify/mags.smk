use rule preprocess__hosts as quantify__mags with:
    input:
        fa_gz=lambda w: features["mag_catalogues"][w.mag_catalogue],
    output:
        fa_gz=QUANT_MAGS / "{mag_catalogue}.fa.gz",
    log:
        QUANT_MAGS / "{mag_catalogue}.log",
    conda:
        "../../environments/mags.yml"


rule quantify__mags__all:
    """Recompress all MAG catalogues"""
    input:
        [QUANT_MAGS / f"{mag_catalogue}.fa.gz" for mag_catalogue in MAG_CATALOGUES],
