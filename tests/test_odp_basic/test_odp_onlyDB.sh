#!/bin/bash

# This script is used to run the odp_basic test case
snakemake --cores 9 -p --snakefile ../../scripts/odp_onlyDB --rerun-incomplete
