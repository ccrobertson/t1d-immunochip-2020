#Tutorial for PEERS: https://github.com/PMBio/peer/wiki/Tutorial
#Description of PEERS: https://www.nature.com/articles/nprot.2011.457

import peer
import scipy as sp
import os
import argparse
import sys

print("Parsing arguments")
parser = argparse.ArgumentParser()
parser.add_argument("DATASET", help="prefix for DATASET to use")
args = parser.parse_args()

print("Reading data")
expr = sp.loadtxt(args.DATASET + '.counts.csv', delimiter=',')
#covs = sp.loadtxt(args.DATASET + '.samples.csv', delimiter=',')
expr.shape
#covs.shape

print("Setting PEER parameters")
model = peer.PEER()
model.setPhenoMean(expr)
model.getPhenoMean().shape
#model.setCovariates(covs)

print("Running PEER")
sys.stdout.flush()
model.setNk(10)
model.getNk()
model.update()

print("Saving PEER results")
factors = model.getX()
factors.shape
sp.savetxt(args.DATASET + ".peer_factors.csv", factors, delimiter=",")

weights = model.getW()
weights.shape
sp.savetxt(args.DATASET + ".peer_weights.csv", weights, delimiter=",")

residuals = model.getResiduals()
residuals.shape
sp.savetxt(args.DATASET + ".peer_residuals.csv", residuals, delimiter=",")

print("Adding PEER results to R object")
cmd = "Rscript ${scripts}/caqtl/add_peer_factors_to_robject.R " + args.DATASET
print(cmd)
os.system(cmd)
