# SubstitutionModels.jl

**Release status:**

[![Latest Release](https://img.shields.io/github/release/BioJulia/SubstitutionModels.jl.svg)](https://github.com/BioJulia/SubstitutionModels.jl/releases/latest)
[![SubstitutionModels](http://pkg.julialang.org/badges/SubstitutionModels_0.6.svg)](http://pkg.julialang.org/?pkg=SubstitutionModels)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/SubstitutionModels.jl/blob/master/LICENSE)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://BioJulia.github.io/SubstitutionModels.jl/stable)
![BioJulia maintainer: jangevaare](https://img.shields.io/badge/BioJulia%20Maintainer-jangevaare-orange.svg)
![BioJulia maintainer: Ward9250](https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg)

**Development status:**

[![Build Status](https://travis-ci.org/BioJulia/SubstitutionModels.jl.svg?branch=master)](https://travis-ci.org/BioJulia/SubstitutionModels.jl)
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

Until SubstitutionModels.jl has had an official release, the current development version can be installed
from the Julia REPL:

```julia
julia> Pkg.clone("https://github.com/BioJulia/SubstitutionModels.jl")
```
