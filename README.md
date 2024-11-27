# Method_ranking

This julia package was created for the automation of chomatography methods selection.

The project is 



### Installation
```julia

# Install dependencies
using Pkg
Pkg.add("PyCall")
Pkg.add("Conda")

using PyCall
using Conda
Conda.add("padelpy")        # Package to calculate molecular fingerprints
Conda.add("pubchempy")      # Package to calculate canonical SMILES from InCHiKey
Conda.add("catboost")       # Package to run CatBoost classification models for prediction of retention behaviour

using Pkg
Pkg.build("PyCall")         # Build PyCall with the installed packages

# --Restart julia--

using Pkg
Pkg.add(url="https://github.com/pockos56/Method_ranking.jl")

```

### Documentation

#### Ranking methods based on their suitability for a specific set of chemicals

E.g. For 10 chemicals:

1st way: Providing canonical SMILES
```julia
using Method_ranking
comp_list_1 = ["RYYVLZVUVIJVGH-UHFFFAOYSA-N", "IISBACLAFKSPIT-UHFFFAOYSA-N", "COLNVLDHVKWLRT-QMMMGPOBSA-N", "ZYGHJZDHTFUPRJ-UHFFFAOYSA-N", "CMPQUABWPXYYSH-UHFFFAOYSA-N", "MUMGGOZAMZWBJJ-DYKIIFRCSA-N", "CBCKQZAAMUWICA-UHFFFAOYSA-N", "OYGQVDSRYXATEL-UHFFFAOYSA-N", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", "UFWIBTONFRDIAS-UHFFFAOYSA-N"]
scores_1 = rank_methods(comp_list)

```

2nd way: Providing InChIKeys
```julia
using Method_ranking
comp_list_2 = ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O", "C1=CC=C(C=C1)C[C@@H](C(=O)O)N", "C1=CC=C2C(=C1)C=CC(=O)O2", "C1=CC=C(C=C1)OP(=O)(O)O", "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C", "C1=CC(=CC=C1N)N", "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F", "CC(=O)OC1=CC=CC=C1C(=O)O", "C1=CC=C2C=CC=CC2=C1"]
scores_2 = rank_methods(comp_list)

```
Note