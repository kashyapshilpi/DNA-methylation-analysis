# DNA Methylation Analysis Repository

Welcome to the DNA Methylation Analysis repository! This repository contains scripts and documentation for the analysis of DNA methylation data. DNA methylation plays a crucial role in regulating gene expression and is an essential mechanism in various biological processes. This project provides a comprehensive guide to the analysis of DNA methylation data, including preprocessing, quality control, differential methylation analysis, and visualization.

## Table of Contents

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Data](#data)
4. [Analysis Workflow](#analysis-workflow)
5. [Usage](#usage)
6. [Results](#results)
7. [License](#license)
8. [Contributing](#contributing)

## Introduction

DNA methylation is an epigenetic modification that involves the addition of a methyl group to the DNA molecule. It plays a critical role in gene regulation and can provide insights into various biological processes and diseases. This repository offers a step-by-step guide to analyzing DNA methylation data, from raw data processing to advanced analysis.

## Dependencies

Before running the DNA methylation analysis pipeline, ensure you have the following tools and software installed on your system:

- Bismark (version 0.15.0): A tool for mapping bisulfite-converted DNA sequences and methylation calling.
- GeneDMRs (version 1.1.0): A package for identifying differentially methylated regions.
- R (version >= 3.0): A programming language for statistical computing and graphics.
- Additional R packages (e.g., `methylKit`, `minfi`, `DMRcate`) for advanced analysis and visualization.

Make sure to update the versions with the appropriate ones you are using.

## Data

The raw DNA methylation data used in this analysis is not included in this repository. Please ensure you have access to the data and place it in the appropriate input directory before running the pipeline.

## Analysis Workflow

The analysis pipeline can be divided into several key steps:

1. **Preprocessing**: Process raw DNA methylation data, including quality control and alignment to a reference genome.
2. **Methylation Calling**: Call methylation levels at individual CpG sites.
3. **Differential Methylation Analysis**: Identify differentially methylated CpG sites or regions between experimental conditions.
4. **Functional Analysis**: Associate differentially methylated regions with genes and perform functional enrichment analysis.
5. **Visualization**: Visualize methylation patterns, differentially methylated regions, and results of functional analysis.

## Usage

To perform DNA methylation analysis, follow these steps:

1. Clone this repository to your local machine: `git clone https://github.com/kashyapshilpi/DNA-methylation-analysis.git`
2. Place the raw DNA methylation data in the appropriate input directory.
3. Install the required dependencies listed in the "Dependencies" section.
4. Modify the configuration file if necessary to specify parameters and options.
5. Execute the analysis script: `bash methylkit.sh`, `bash Dna_methylation_script.sh`, `Rscript GeneDMARs.R`, `Rscript MethylKit_Script_15032023.R`
6. The results, including plots and analysis outputs, will be generated in the output directory.

## Results

The results of the DNA methylation analysis will be stored in the output directory. This will include methylation profiles, differential methylation analysis results, visualization plots, and functional analysis outcomes.

## License

This project is licensed under the [MIT License](license.txt).

## Contributing

If you wish to contribute to this project, feel free to open issues, submit pull requests, or suggest improvements. We welcome your contributions!

Thank you for using this DNA Methylation Analysis repository. If you have any questions or encounter any issues, please don't hesitate to contact us.

Happy analyzing!
