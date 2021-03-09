import json

with open("config.json", "r") as f:
  config = json.load(f)

DATADIR = config['datadir']
RESULTSDIR = config['resultsdir']

rule all:
	input:

rule data:
	output:
		'{DATADIR}/AFR_1kg_chr22/ldobj_22_21417539_22877794.rds'
	shell:
		'Rscript scripts/get_ldref.r {DATADIR}'

rule run_simulations:
	input:
		data.output
	output:
		'{RESULTSDIR}/res.rdata'
	shell:
		'Rscript scripts/simulate.r {data.output}'


rule analysis:
	input:
		run_simulations.output,
		'docs/analysis.rmd'
	output:
		'docs/analysis.html'
	shell:
		"""
		'Rscript -e "rmarkdown::render('docs/analysis.rmd')"
		"""
