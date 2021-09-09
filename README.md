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



---

```r
set.seed(1234)
source("scripts/simulate.r")


ldobjfile <- "~/data/ld_files/ldmat/AFR_1kg_chr22/ldobj_22_22877795_23154057.rds"
ldobj <- readRDS(ldobjfile)

ldobj <- list(ld=ldobj$ld[1:1000, 1:1000], map=ldobj$map[1:1000,], nref=ldobj$nref)
d <- summary_data(ldobj, 1000, 1000, 0, 0, 1, 0.8, 0.8, radius=20000)
res <- run_coloc(d, susie.args=list(nref=ldobj$nref, check_R=FALSE))

print(res$res$summary)
sensitivity(res$res,"H4 > 0.9",row=1,dataset1=res$d1,dataset2=res$d2)
sensitivity(res$res,"H4 > 0.9",row=2,dataset1=res$d1,dataset2=res$d2)

view_res(res)
```
