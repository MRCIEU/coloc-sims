library(tidyverse)
library(coloc)
library(simulateGP)
library(parallel)
library(jsonlite)

sample_causal_variants <- function(map, ncausal, h2, S=0, radius=500000)
{
	# get range given radius
	maprange <- range(map$pos) %>% diff
	if(maprange < (radius*2))
	{
		radius <- floor((maprange-1)/2)
		warning("radius is too small, resetting to ", radius)
	}
	nsnp <- nrow(map)
	rang <- map$pos[nsnp] - map$pos[1]
	window <- c(map$pos[1] + radius, map$pos[nsnp] - radius)
	windowi <- c(which(map$pos > window[1])[1], which(map$pos > window[2])[1]-1)

	if(ncausal==0)
	{
		map$selected <- 1:nrow(map) %in% window
		map$causal <- FALSE
		map$beta <- 0
		return(map)
	}

	# get the first one somewhere within range
	causals <- rep(0, ncausal)
	causals[1] <- sample(windowi[1]:windowi[2], 1)

	# define window
	pos <- map$pos[causals[1]]
	window <- which(map$pos > (pos - radius) & map$pos < (pos + radius))

	# sample others from within window
	if(ncausal > 1)
	{
		causals[2:ncausal] <- sample(window[!window == causals[1]], ncausal - 1, replace=FALSE)
	}
	print(causals)
	p <- simulateGP::generate_gwas_params(map[causals,], h2 = h2, S=S)
	map$selected <- 1:nrow(map) %in% window
	map$causal <- 1:nrow(map) %in% causals
	map$beta <- 0
	map$beta[map$causal] <- p$beta
	return(map)
}

# 1. Choose causal effects for traits 1 and 2
# 2. Add LD
# 3. Generate sample estimates given sample size
summary_data <- function(ldobj, nid1, nid2, ndistinct1, ndistinct2, nshared, h2_1, h2_2, S_1=0, S_2=0, radius=500000)
{
	# How many causal variants in total
	ntotal_1 <- ndistinct1 + nshared
	ntotal_2 <- ndistinct2 + nshared

	distinct1 <- sample_causal_variants(ldobj$map, ndistinct1, h2_1/ntotal_1*ndistinct1, S=S_1, radius=radius)
	distinct2 <- sample_causal_variants(ldobj$map, ndistinct2, h2_2/ntotal_2*ndistinct2, S=S_2, radius=radius)
	shared1 <- sample_causal_variants(ldobj$map, nshared, h2_1/ntotal_1*nshared, S=S_1, radius=radius)
	shared2 <- shared1
	shared2$beta <- sqrt(shared2$beta^2 / (h2_1/h2_2)) * sign(shared2$beta)
	trait1 <- distinct1
	trait1$beta <- trait1$beta + shared1$beta
	trait2 <- distinct2
	trait2$beta <- trait2$beta + shared2$beta

	trait1 <- trait1 %>%
		add_ld_to_params(ldobj=ldobj) %>%
		generate_gwas_ss(nid1)
	trait2 <- trait2 %>%
		add_ld_to_params(ldobj=ldobj) %>%
		generate_gwas_ss(nid2)

	ld <- ldobj$ld
	rownames(ld) <- colnames(ld) <- trait1$snp

	return(list(trait1=trait1, trait2=trait2, ld=ld))
}

run_coloc <- function(summary_data, susie.args)
{
	d1 <- list(
		pvalues = summary_data$trait1$pval,
		N=summary_data$trait1$n,
		MAF=summary_data$trait1$af,
		beta=summary_data$trait1$bhat,
		varbeta=summary_data$trait1$se^2,
		snp=summary_data$trait1$snp,
		position=summary_data$trait1$pos,
		sdY=1,
		type="quant",
		LD=summary_data$ld
	)
	d2 <- list(
		pvalues = summary_data$trait2$pval,
		N=summary_data$trait2$n,
		MAF=summary_data$trait2$af,
		beta=summary_data$trait2$bhat,
		varbeta=summary_data$trait2$se^2,
		snp=summary_data$trait2$snp,
		position=summary_data$trait2$pos,
		sdY=1,
		type="quant",
		LD=summary_data$ld
	)

	s1 <- do.call("runsusie", c(list(d=d1, suffix=1), susie.args))
	s2 <- do.call("runsusie", c(list(d=d2, suffix=1), susie.args))

	# coloc.abf(d1, d2)
	res <- coloc.susie(s1, s2)
	return(list(d1=d1, d2=d2, s1=s1, s2=s2, res=res))
}


evaluate_performance <- function(input, output, regionsize)
{
	# count how many distinct and shared causal variants for the two traits within a given window size
}


simulation <- function(param, ld)
{
	ld <- ld[[param$region]]
	map1 <- sample_causal_variants(ld$map, param$ncausal, param$rsq_trait1)
	ss1 <- simulate_summary_data(ld$ld, map1, param$nid)
	if(param$coloc)
	{
		map2 <- map1
		map2$b <- map2$b * sqrt(param$rsq_trait2) / sqrt(param$rsq_trait1)
	} else {
		map2 <- sample_causal_variants(ld$map, param$ncausal, param$rsq_trait2)
	}
	ss2 <- simulate_summary_data(ld$ld, map2, param$nid)
	param$coloc_result <- run_coloc(ss1, ss2, param$coverage, "sparse")$summary %>% {which.max(.[-1])}
	return(param)
}


view_res <- function(res)
{
	out <- lapply(1:2, function(x)
	{
		cs <- summary(res[[paste0("s", x)]])$cs
		res[[paste0("d", x)]]$susie <- 1:length(res[[paste0("d", x)]]$beta) %in% cs$variable
		tibble(
			pos=res[[paste0("d", x)]]$position,
			pvalues=res[[paste0("d", x)]]$pvalues,
			susie=res[[paste0("d", x)]]$susie,
			trait=x
		)
	}) %>% bind_rows() %>%
	ggplot(., aes(pos, -log10(pvalues))) +
	geom_point(aes(colour=susie)) +
	facet_grid(trait ~ .)
	print(out)
}

###################

args <- commandArgs(T)
ldobjfile <- args[1]

ldobj <- readRDS(ldobjfile)
ldobj <- list(ld=ldobj$ld[1:1000, 1:1000], map=ldobj$map[1:1000,], nref=ldobj$nref)

source("scripts/simulate.r")
d <- summary_data(ldobj, 1000, 1000, 0, 0, 1, 0.8, 0.8, radius=20000)

res <- run_coloc(d, susie.args=list(nref=ldobj$nref, check_R=FALSE))


print(res$res$summary)
sensitivity(res$res,"H4 > 0.9",row=1,dataset1=res$d1,dataset2=res$d2)
sensitivity(res$res,"H4 > 0.9",row=2,dataset1=res$d1,dataset2=res$d2)

view_res(res)


