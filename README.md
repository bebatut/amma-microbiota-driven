Atlas of Microglia-Microbiota in Aging (AMMA)
=============================================

Microglia, the brain resident macrophages, display high plasticity in response to their environment. Aging of the central nervous system (CNS), where microglial physiology is especially disrupted, is a major risk factor for a myriad of neurodegenerative diseases. Therefore, it is crucial to decipher intrinsic and extrinsic factors, like sex and the microbiome, that potentially modulate this process.

We found that **microglia follow sex-dependent dynamics in aging**. This repository stores the transcriptomics data analyses and the sources for the website explaining the analysis.

# Transcriptomics data analyses

## Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html)



- Creation of the `conda` environment with all the requirements

    ```
    $ make create-env
    ```

## Preparation of files from the sequencing facility

1. Rename the files from the sequencing facility to follow the naming convention

    ```
    $ python src/copy_rename_raw_files.py \
        --input_dir <path to input directory> \
        --file_name_description <path to csv file with the correspondance between directory structure and sample name (from the Google drive)> \
        --output_dir <path to output directory>\
    ```

The naming convention for each sample is `microbiota_age_sex_replicate`

## From sequences to gene counts (in Galaxy)

1. Upload the data on Galaxy (e.g. [https://usegalaxy.eu/](https://usegalaxy.eu/)) inside a data library
2. Update the details in [`config.yaml`](config.yaml), specially the API key
3. Prepare the history in Galaxy (import the files from the data library, merge the files sequenced on 2 different lanes (for Project_S178 and Project_S225) and move the input files into collections)

    ```
    $ python src/prepare_data.py
    ```

5. Launch Galaxy workflow to extract gene counts

    ```
    $ python src/extract_gene_counts.py
    ```

    The worklow do:
    1. Quality control and trimming using FastQC and Trim Galore!
    2. Preliminary mapping and experiment inference using STAR and RSeQC
    3. Mapping using STAR
    4. Gene counting using FeatureCounts

The workflow is applied on each dataset (organized into data collection). It can take a while.

Once it is finished:

1. Download the generated count table and the gene length file
2. Put these files in the `data` folder

## Differentially Expression Analysis (locally using Jupyter Notebooks)

1. Launch Jupyter

    ```
    $ make launch-jupyter
    ```

3. Open [http://localhost:8888/tree/](http://localhost:8888/tree/)
2. Move to `src` in Jupyter
4. Prepare the data: Open `0-prepare_data.ipynb` and execute all cells

### Full analysis

Data: SPF/GF, young/middle-aged/old, female/male

1. Move to `full-data`

2. Prepare the differential expression analysis
    1. Open `1-run_dge_analysis.ipynb` and execute all cells
    2. Open `2-previsualize_data.ipynb` and execute all cells

2. Analyze the differentially expressed genes given different comparisons: Open the related notebook and execute all cells
    
    Analysis | Notebook
    --- | ---
    Effect of microbiota (GF vs SPF) for the different ages and sexes | `3-analyze_microbiota_effect_given_ages_sexes.ipynb`
    Effect of sexes (Male vs Female) for the different ages and microbiota | `4-analyze_sex_effect_given_ages_microbiota.ipynb`
    Effect of ages (Middle-aged vs Young, Old vs Young and Old vs Middle-aged) for the different microbiota and sexes | `6-analyze_age_effect_given_microbiota_sexes.ipynb`

4. Run extra analyses

    Analysis | Noteobook
    --- | ---
    Study of CML effect | `7-analyze_cml_effect.ipynb`
    Postvisualize | `8-postvisualize.ipynb`

## Sex-driven aging analysis

Data: SPF, young/middle-aged/old, female/male

1. Move to `sex-driven-aging`

2. Prepare the differential expression analysis
    1. Open `1-extract_samples` and execute all cells
    1. Open `2-run_dge_analysis.ipynb` and execute all cells
    2. Open `3-previsualize_data.ipynb` and execute all cells

3. Analyze the differentially expressed genes given different comparisons: open the related notebook and execute all cells

    Analysis | Notebook
    --- | ---
    Effect of the sex (Male vs Female) for the 3 ages | `4-analyze-sex-effect-given-ages.ipynb`
    Effect of the ages (Young, Middle-aged, Old) for the sexes | `4-analyze-age-effect-given-sex.ipynb`


## Microbiota driven analysis

Data: SPF/GF, young/old, male

1. Move to `microbiota-driven`

2. Prepare the differential expression analysis
    1. Open `1-extract_samples` and execute all cells
    1. Open `2-run_dge_analysis.ipynb` and execute all cells
    2. Open `3-previsualize_data.ipynb` and execute all cells

3. Analyze the differentially expressed genes given different comparisons: open the related notebook and execute all cells

    Analysis | Notebook
    --- | ---
    Effect of the sex (Male vs Female) for the 3 ages | `4-analyze-sex-effect-given-ages.ipynb`
    Effect of the ages (Young, Middle-aged, Old) for the sexes | `4-analyze-age-effect-given-sex.ipynb`


# Website

This folder stores the sources of the website describing the analyses in `docs` folder.
Reports from the Jupyter notebooks are available there to show the different steps and images. 

## Generate HTML reports from the Jupyter Notebooks

```
$ make generate-reports
```

These reports are stored in the `docs` folder and are linked on the website.

## Generate the website locally

1. Install Install the website's dependencies:

    ```
    $ make install
    ```

2. Serve the website locally

    ```
    $ make serve
    ```

3. Open [http://127.0.0.1:4000/amma/](http://127.0.0.1:4000/amma/)