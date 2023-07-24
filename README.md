# Platereader screen analysis pipeline

Nextflow pipeline to process files from Biotek platereaders, normalise output, perform QC, make plots, and do statistical tests. 

## Processing steps

For each experiment in the `sample_sheet`:

1. Convert the compound source map to columns using `hts pivot`.
2. Convert platereader files to columnar format using `hts parse` and concatenate into one big table. 
3. Annotate compounds the platereader data with `sample_sheet` and pivoted compound source map using `hts join`.
4. Normalize the data using positive and negative controls using `hts normalize`.
5. Do statstical tests compared to negative controls using `hts sumamrize`.
6. (Optionally) if `counterscreen` is set, do statistical tests between another set of groups.

### Other steps

1. Carry out quality control and generate plots of Z\'-factor using `hts qc`.
2. Plot plate heatmaps using `hts plot-hm`.
3. Plot replicate scatter plots using `hts plot-rep`.
4. Plot histogram of data using `hts plot-hist`.

## Requirements

### Software

You need to have Nextflow and either `conda` or `mamba` installed on your system. If possible, use `mamba` because it will be faster.

### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow
```

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-platereader-screen
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which will not be automatically updated. If you want to ensure that you're using the very latest version of the pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-platereader-screen -latest
```
If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash 
nextflow run scbirlab/nf-platereader-screen -r v0.0.1
```

For help, use `nextflow run scbirlab/nf-platereader-screen --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: Path to a CSV with information about the samples and platereader files to be processed
- `data_dir`: Path to directory where files in `sample_sheet` can be found
- `output_dir`: Path to directory where outputs should be saved.
- `export_layout`: "plate" or "row", indicating the layout of the exported platereader data
- `positive`: Positive control label contained in the `control_column`
- `negative`: Negative control label contained in the `control_column`
- `grouping`: Columns with labels indicating groups (i.e batches) within which to normalize data. Often this will be a column of plate labels.
- `hit_grouping`: Columns with labels indicating repeated measurements of experimental conditions. Often this will be a column indicating different compounds.

The following parameters have default values can be overridden if necessary.

- `layout_content_name = "compound_name"`: A name to give the data from the wells of the compound source plates
- `control_column = "compound_name"`: Name of column containing labels indicating positive and negative controls

The following parameters are optional.

- `counterscreen`: A colon-separated column name and negative control label to do a second wave of hit calling. For example `"strain_name:WT"`.

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    sample_sheet = "path/to/sample-sheet.csv"
    data_dir = "path/to/data"
    output_dir = "path/to/outputs"

    export_layout = 'plate'

    positive = 'RIF'
    negative = 'DMSO'
    grouping = 'strain_name plate_id'

    hit_grouping = 'strain_name compound_name'
    counterscreen = 'strain_name:WT'

}
```

Alternatively, you can provide these on the command line:

```bash
nextflow run scbirlab/nf-platereader-screen \
    --sample_sheet path/to/sample_sheet.csv \
    --data_dir path/to/data \
    --output_dir path/to/outputs \
    --export_layout plate \
    --positive RIF \
    --negative DMSO \
    --grouping strain_name plate_id \
    --hit_grouping strain_name compound_name \
    --counterscreen strain_name:WT
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which platereader files belong to which named sample, which compound map belongs to each, and further optional metadata, such as timepoints, which could be used in downstream analysis.

The file must have a header with the column names below, and one line per sample to be processed.

- `experiment_id`: name of the experiment, under which to group data from multiple files. Multiple experiments can be processed at once, and the pipeline will analyse their data separately.
- `data_filename`: filename of the platereader-exported file.
- `compound_source_filename`: filename of the compound map containing the compounds in the plates of the platereder file.
- `compound_source_plate_id`: the plate name from the `compound_source_filename` which was specifically used to generate data in this platereader file.

Additional columns can be included. These will be included in output tables, so can be used for downstrwam analysis.

Here is an example of a couple of lines from a sample sheet:

| experiment_id | data_filename | compound_source_filename | compound_source_plate_id | strain_name | timepoint |
| ------------- | ------------- | ------------------------ | ------------------------ | ----------- | --------- |
| expt01        | 170423_plate 7 rep 1.txt      | compounds.xlsx | Plate 7 | WT | 14 |
| expt01        | 170423_plate 8 rep 1.txt      | compounds.xlsx | Plate 8 | WT | 14 |

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under four directories:

- `normalized`: Normalized data extracted from platereader-exported files.
- `qc`: QC plots and tables for Z\'-factor and SSMD.
- `plots`: Plate heatmaps, replicates, and histograms.
- `hits`: Statsistical testign tables and plots.

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-platereader-screen/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [hts-tools](https://github.com/scbirlab/hts-tools)