#!/bin/bash

# basic example to create a new sequence

org='fr2' # fr2 = Fugu genome

# delete previous model if exist
rm -rf data/$org

# creating the genome model
perl bin/createModel.pl -v -m $org --rm_tmp

# create a new sequence
#perl bin/create/FakeSequence.pl -m $org
