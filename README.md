# read_depth_genotyper
Get genotypes for regions using WSSD and SUNK

## Install
`git clone https://github.com/bnelsj/read_depth_genotyper --recursive`

## Genotype
`snakesub -j 200 -w 60 -kT get_tables`

## Plot
`snakesub -j 200 -w 60 -kT`
