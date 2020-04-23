
args = commandArgs(trailingOnly=TRUE)

setwd(Sys.getenv("freeze"))


getPEERData = function(prefix) {
  #read in R data
  dat_sub = readRDS(paste0(prefix,".rds"))
  #read in PEER factors
  peer = read.csv(paste0(prefix,".peer_factors.csv"), header=FALSE)
  names(peer) = paste0("PF", seq(10))
  dat_sub[["peer"]] = peer
  #save
  saveRDS(dat_sub, paste0(prefix,".rds"))
}


getPEERData(args[1])
