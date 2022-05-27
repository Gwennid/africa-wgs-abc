## Simulating genetic datasets for ABC

We simulated genetic datasets with fastsimcoal2 (Excoffier and Foll, [2011](https://doi.org/10.1093/bioinformatics/btr124), Excoffier et al., [2013](https://doi.org/10.1371/journal.pgen.1003905)). The genetic model was based on the model in (Jay et al., [2019](https://doi.org/10.1093/molbev/msz038)).

Vectors of parameters were generated with a custom Python script. For each simulation, a fastsimcoal input file (.est) was created by taking the values in a vector and replace the placeholders in the template .est file. The resulting .est file was fed to fastsimcoal together with the template (.tpl) file. The .tpl files, the templates .est, as well as the vectors of parameters for each model are in the subfolder `continuous` and `pulse`.

The outputs from fastsimcoal, in .arp format, were converted to plink TPED format and to a set of Python objects in order to calculate summary statistics. This is done with script `xxx`.

Summary statistics were calculated with plink (Purcell et al., [2007](https://doi.org/10.1086/519795)) v1.90b4.9, R v3.6.1, Python v2.7.15, bash, awk, scripts by Jay et al., [2019](https://doi.org/10.1093/molbev/msz038) (accessible [here](https://gitlab.inria.fr/ml_genetics/public/demoseq/-/tree/master)), and the [asd software](https://github.com/szpiech/asd).

## TODO

- [ ] Include: custom Python script to generate vectors of parameters (do I have it?)
- [ ] Include: script to convert fsc output (or should it go somewhere else?
- [ ] Example of script to perform the full simulation?
- [ ] What else?
