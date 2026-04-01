# PerturbGen Reproducibility

This repository contains the code to reproduce the analyses and benchmarking experiments performed in the PerturbGen [preprint](https://www.biorxiv.org/content/early/2026/03/05/2026.03.04.709254).
The PerturbGen source code can be found [here](https://github.com/Lotfollahi-lab/Perturbgen).

# Notebook structure

All analyses are implemented as Jupyter notebooks organized by experiment.  
Within each folder, notebooks should be run in consecutive order, as indicated by their numerical prefixes.

- `LPS_analysis/`  
  Reproduces **Figure 2**

- `HSPC_analysis/tf_screen/`  
  Reproduces the **transcription factor perturbation benchmark (Figure 3)**

- `HSPC_analysis/perturbation_atlas/`  
  Demonstrates how perturbation atlases were constructed (**Figure 3**)

- `HSPC_analysis/open_targets/`  
  Maps perturbation programs to genetic traits (**Figure 4**)

- `HSPC_analysis/etv6/`  
  Reproduces the **monogenic disease analysis (Figure 4)**

- `SkO_analysis/`  
  Reproduces analyses for **Figures 5 and 6**

### Notes

- Notebooks are intended to be run from top to bottom in each directory  
- Intermediate outputs are saved within their respective folders or in a shared `results/` directory  
- Some analyses may require substantial compute resources (e.g. GPU / high memory)  

## Reference
```
@article {Chi Hao Ly2026.03.04.709254,
	author = {Chi Hao Ly, Kevin and Miraki Feriz, Adib and Isobe, Tomoya and Vahidi, Amirhossein and Vaghari, Delshad and Rostron, Anthony et al.},
	title = {Predicting how perturbations reshape cellular trajectories with PerturbGen},
	year = {2026},
	doi = {10.64898/2026.03.04.709254},
	URL = {https://www.biorxiv.org/content/early/2026/03/05/2026.03.04.709254},
	journal = {bioRxiv}
}
```
