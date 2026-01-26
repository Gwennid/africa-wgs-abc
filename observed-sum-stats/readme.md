# Calculating observed summary statistics for ABC

Note that the scripts provided here were written to run on a cluster and hence all paths would need to be updated in order to run.

The entire process is done via the script `observed-sum-stats/select-n-windows_calculate-sumstats_5Ju-5Ka-5NP-5WP-5EP.sh`. This is a template that requires adjustments to run, in particular the paths, name of the population combination and individual names (replace the variable "PLACEHOLDER" by a list of 15 individuals corresponding to 5 wRHG, 5 eRHG and 5 RHGn. This script calls to a number of scripts that are provided in `observed-sum-stats/helper_scripts`.

## Select genomic regions

We prepared a call-set of high-quality regions by applying the 1000 Genomes phase 3 accessibility mask and filtering out indels on the variant callset comprising all individuals. This resulted in a “fragmented” genome, with fragments of high-quality sites and fragments that are filtered out. We then agglutinated high-quality fragments spanning 1 Mb in total agglutinated-length. We then calculated the mapped bp distance between the first and the last fragment of each 1Mb window and then randomly chose 100 windows spanning less than 1.2 mapped Mb and located on chromosomes 1 to 10, among the 824 such windows identified on chromosomes 1 to 22. We then extracted these 100 independent autosomal loci of 1 Mb each for the 54 sets of 25 individuals.

## Calculate summary statistics

We are calculating the same summary statistics as in the simulated data. However, since this is observed data, we need to manipulate it differently. Importantly, we infer ancestral state for a subset of variants to calculate the summary statistics that are dependent on the ancestral state for example the site frequency spectrum. To do that, we used an alignment of three non-human apes (gorilla – reference genome gorGor5, chimpanzee – reference genome panTro6 - and orangutan – reference genome ponAbe3). We considered that the ancestral state of a site could be reliably inferred if the site had the same allele (A, C, T or G) in the three apes. Only these sites were used to calculate the proportion of homozygous ancestral positions per population and the unfolded standardized SFS.
