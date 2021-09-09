library(here)
library(tidyverse)

# Load the functions
source(here("scripts", "simulation-functions.r"))

# Load the pre-generated LD region map file
load(here("data", "ld", "maps.rdata"))

# Get the list of regions
regions <- unique(mapEUR$region)

# Generate the set of simulations to run
simulation_parameters <- expand.grid(
	ldregion=regions,
	ndistinct1 = c(1,2,3),
	ndistinct2 = c(1,2,3),
	nshared = c(1,2,3),
	nid = 10000,
	nrep = 1:100,
	h2_1 = c(0.1, 0.8),
	h2_2 = c(0.1, 0.8),
	S_1 = c(),
	S_2 = c()
)
simulation_parameters$simid <- 1:nrow(simulation_parameters)

# These are command line arguments to help split the set of simulations up across a cluster
args <- commandArgs(T)
chunk <- as.numeric(args[1]) # which chunk to run (starting from 0)
chunksize <- as.numeric(args[2]) # How large is each chunk
nchunk <- ceiling(nrow(simulation_parameters) / chunksize)
start <- chunk * chunksize + 1
end <- min((chunk + 1) * chunksize, nrow(simulation_parameters))

message("total:", nrow(simulation_parameters))
message("start: ", start)
message("end: ", end)

simulation_parameters <- simulation_parameters[start:end,]

# If your simulation ends early then you can try to resume where it left off. 
output <- here("results", paste0(chunk, ".rdata"))
output_int <- paste0(output, ".int")

if(file.exists(output))
{
  message("already complete")
  q()
}

if(file.exists(output_int))
{
  message("Previous run already exists, resuming...")
  load(output_int)
  a <- max(sapply(oinst1, function(x) x$simid))
  message(sum(simulation_parameters$simid > a), " out of ", nrow(simulation_parameters), " remaining")
  j <- which(simulation_parameters$simid > a)[1]
} else {
  message("New run")
  oinst1 <- list()
  omr1 <- list()
  j <- 1
}


res <- list()
for(i in j:nrow(param))
{
  message(i, " of ", nrow(param))
  res[[i]] <- tryCatch(do.call(sim, args=param[i,]), error=function(e) {print(e); return(NULL)})
  # intermediate save
  save(res, file=output_int)
}

save(res, file=output)
