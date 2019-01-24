To use this workflow you need to...

- Install MiniConda: Use duckduckgo or whatever to find and install it on your system

- When conda is installed you need to install snakemake. Use the following command:

    conda install -c bioconda -c conda-forge snakemake

- So now we are ready and can analyze your data. In config.yaml you can change locations to pointing to your files or change other things. Every config is explained.

- When everything is done you can start analyzing. Use

    snakemake --use-conda

    in the cmdline where the _Snakefile_ is located to start analyzing.

- Results should be located in the results folder.
