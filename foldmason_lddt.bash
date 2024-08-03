#!/bin/bash -e

export PATH=$PATH:/z/sw/foldmason/foldmason/bin/

cd pair

foldmason createdb . pairdb
foldmason msa2lddt pairdb pair.afa
