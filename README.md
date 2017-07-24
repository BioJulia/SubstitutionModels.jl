# SubstitutionModels.jl

**Development status:**

[![Build Status](https://travis-ci.org/BioJulia/SubstitutionModels.jl.svg?branch=master)](https://travis-ci.org/BioJulia/SubstitutionModels.jl) [![Coverage Status](https://coveralls.io/repos/BioJulia/SubstitutionModels.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/BioJulia/SubstitutionModels.jl?branch=master) [![codecov.io](http://codecov.io/github/BioJulia/SubstitutionModels.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/SubstitutionModels.jl?branch=master) [![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/SubstitutionModels.jl/blob/master/LICENSE) [![](https://img.shields.io/badge/docs-master-blue.svg)](https://biojulia.github.io/SubstitutionModels.jl/master) ![BioJulia maintainer: jangevaare](https://img.shields.io/badge/BioJulia%20Maintainer-jangevaare-orange.svg) ![BioJulia maintainer: Ward9250](https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg)

## Description

SubstitutionModels.jl provides facilities to model the substitution process of
biological sequences. Such models are essential for the analysis of sequence
evolution, phylogenetics, and simulation.

We first aim to provide the most common substitution models
used in the literature, but aim to build an extendable framework using julia's
type system and traits, so as custom model types can be created and used.

## Installation

Until SubstitutionModels.jl has had an official release, the current development version can be installed
from the Julia REPL:

```julia
julia> Pkg.clone("https://github.com/BioJulia/SubstitutionModels.jl")
```

## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [contributing section of the documentation](biojulia.github.io/SubstitutionModels.jl/stable/contributing)
for more information.

If you have a question about
contributing or using this package, you are encouraged to use the
[Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio).
