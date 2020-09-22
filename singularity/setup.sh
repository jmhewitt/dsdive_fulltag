#!/bin/bash

echo Downloading singularity image...

curl -O https://research-singularity-registry.oit.duke.edu/jmh182/rmovement.sif


echo Making libs directory...

if [ ! -d "./libs" ];
then
  mkdir libs
fi

echo Finished setting up singularity image
