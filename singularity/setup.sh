#!/bin/bash

echo Downloading singularity image...

curl -O https://research-singularity-registry.oit.duke.edu/jmh182/rmovement.sif


echo Making libs directory...

if [ ! -d "./libs" ];
then
  mkdir libs
fi

echo Installing packages into singularity image...

singularity exec rmovement.sif install2.r -l libs -s dgof
singularity exec rmovement.sif install2.r -l libs -s kableExtra
singularity exec rmovement.sif install2.r -l libs -s KSgeneral

echo Finished setting up singularity image
