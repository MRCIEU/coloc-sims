# Colocalisation simulations

## Setup

Create `config.json` that specifies where the data and results directories are

```json
{
    "datadir": "/path/to/data/dir",
    "resultsdir": "/path/to/results/dir"
}
```

Install the R packages using

```r
renv::renv_restore()
```


## Running simulations

Run `snakemake`.

