This folder contains the files used for simulating genetic data under a continuous migration model with fastsimcoal. The folders are as follow:

- `def`: the three subfolders `high`, `intermediate` and `low` contain files with vectors of parameters (population size, migration rates etc). Each row corresponds to a simulation. `high` corresponds to migration rates in [0,0.05]; `intermediate` in [0,0.0125]; and `low` in [0,0.001]. The `.def_header` files contain the header row of the `.def` and are identical across migration intensity.
- `est`: these files are given as an example. They contain, for each parameter, information on the prior distribution (for example uniform distribution between 20 and 20000 for population size). These were not used in the simulations.
- `tpl`: these files contains the genetic architecture and demographic model used by fastsimcoal. They are combined with a vector from the `.def` to simulate data.
