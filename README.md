# SubstitutionModels.jl

**Release status:**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3988662.svg)](https://doi.org/10.5281/zenodo.3988662)
[![Latest Release](https://img.shields.io/github/release/BioJulia/SubstitutionModels.jl.svg)](https://github.com/BioJulia/SubstitutionModels.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/SubstitutionModels.jl/blob/master/LICENSE)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://BioJulia.github.io/SubstitutionModels.jl/stable)
![BioJulia maintainer: jangevaare](https://img.shields.io/badge/BioJulia%20Maintainer-jangevaare-orange.svg)

**Development status:**

[![test-lts](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-lts.yml/badge.svg)](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-lts.yml)
[![test-stable](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-stable.yml/badge.svg)](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-stable.yml)
[![test-nightly](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-nightly.yml/badge.svg)](https://github.com/BioJulia/SubstitutionModels.jl/actions/workflows/test-nightly.yml)
[![codecov.io](http://codecov.io/github/BioJulia/SubstitutionModels.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/SubstitutionModels.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://BioJulia.github.io/SubstitutionModels.jl/latest)

## Description

SubstitutionModels.jl provides facilities to model the substitution process of
biological sequences. Such models are essential for the analysis of sequence
evolution, phylogenetics, and simulation.

We first aim to provide the most common substitution models
used in the literature, but aim to build an extendable framework using julia's
type system and traits, so as custom model types can be created and used.

## Installation

The current release version can be installed
from the Julia REPL:

```julia
using Pkg
Pkg.add("SubstitutionModels")
```

## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [contributing section of the documentation](http://biojulia.dev/Contributing/latest/)
for more information.

If you have a question about
contributing or using this package, you are encouraged to use the
[Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio).
