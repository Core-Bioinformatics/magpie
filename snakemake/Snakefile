include_barcode_matrix = False

if os.path.isfile('input/selected.txt'):
    file = open("input/selected.txt", "r")
    samples = [line.rstrip() for line in file]
else:
    samples = glob_wildcards("input/{sample}/msi").sample
if os.path.isfile('input/exclude.txt'):
    file = open("input/exclude.txt", "r")
    samples = list(set(samples).difference(set([line.rstrip() for line in file])))

if include_barcode_matrix:
    output_file = "spaceranger_meanIntensity"
else:
    output_file = "spaceranger"
    
rule all:
    input:
        expand("output/{sample}/"+output_file+"/filtered_feature_bc_matrix.h5", sample=samples)
        
rule perform_coreg:
    input:
        "input/{sample}/msi/MSI_metadata.csv"
    output:
        "output/{sample}/transformed.csv",
        "output/{sample}/transformed_withCoords.png",
        "output/{sample}/transformed.png"
    params:
        transform_type = 'TPS',
        sample = "{sample}"
    script:
        "scripts/alter_data.py"

rule make_spaceranger:
    input:
        "output/{sample}/transformed.csv",
        "output/{sample}/transformed.png"
    output:
        "output/{sample}/spaceranger/filtered_feature_bc_matrix.h5"
    params:
        sample = "{sample}"
    script:
        "scripts/create_mock_spaceranger.py"

rule create_barcode_matrix:
    input:
        "output/{sample}/spaceranger/filtered_feature_bc_matrix.h5"
    output:
        "output/{sample}/spaceranger_meanIntensity/filtered_feature_bc_matrix.h5"
    params:
        sample = "{sample}"
    script:
        "scripts/create_perbarcode_matrix.py"