# Scripts for calculating heterozygosity

Expected heterozygosity is calculated for each population for the autosomes and for the X chromosome separately.

First the different configurations are counted (for example, how many positions are homozygous reference).

Then, the counts of these configurations are summarized and used to calculate heterozygosity.

The last step is to calculate the X to autosome heterozygosity ratio.

## Input autosomes

Start from `/proj/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/HC_BPresolution/3maskrecal.realn/allsites/3_geno01_hwefiltering` (?)

## Input X chromosome

The PAR region is removed. **TODO Elaborate**

Start from `/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/Gwenna_Xchr/filtering_after_VQSR`
