# Changelog

Tagged entries were reconstructed from the changes between successive tags.
The first release is summarized from its source tree. The `v0.5.0+docs` tag is
included separately because it only changes documentation deployment.

## Unreleased (version set to 0.6.1)

- Remove the optional LinearAlgebra version bound so packages using
  SubstitutionModels can resolve their test environments on Julia 1.0.

## [0.6.0](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.5.1...v0.6.0)

- **Breaking:** Remove direct nucleotide indexing on arbitrary arrays. Use
  `a[nucleotide_index(nt)]` or wrap the array with `NucleotideView(a)`.
- Add the exported `NucleotideView` constructor for four-element vectors and
  4×4 matrices with one-based axes. The wrapper shares storage with its parent,
  supports nucleotide indexing and slices, and handles overlapping assignments.
  Mutation requires a mutable parent; `P` and `Q` retain their static matrix types.
- Print compact model representations as constructor expressions, with descriptive
  output through `text/plain` display. Remove hard-coded terminal escape sequences.
- Limit numerical fallback to expected decomposition failures and unusable
  stationary frequencies; propagate unrelated model errors and interrupts.
- Update documentation tooling and test Julia 1.13, 1.10, and the minimum supported
  Julia 1.0. Keep nightly tests nonblocking.
- Allow overall coverage to decrease by up to 2.5 percentage points and make
  patch coverage informational.

## [0.5.1](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.5.0+docs...v0.5.1)

- Fix unchecked nucleotide indexing: ambiguous symbols and gaps throw
  `ArgumentError`, and out-of-bounds accesses throw `BoundsError`. Rejected writes
  leave the array unchanged.
- Add `nucleotide_index` for explicit conversion to integer indices. Deprecate
  direct nucleotide indexing on arbitrary arrays ahead of its removal in 0.6.0.
- Use checked integer indexing in model constructors, including with `safe=false`.
- Allow parameter and frequency arrays to have different container types in
  constructors and `convert`.
- Document migration to `NucleotideView(a)` in 0.6.0.
- Add regression tests and normal-bounds CI runs. Refresh GitHub Actions, add
  Dependabot updates and Julia caching, connect the Codecov token, and run CI
  on Linux with Ubuntu 24.04 runners.

## [0.5.0+docs](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.5.0...v0.5.0+docs)

- Trigger documentation builds when a GitHub release is published.
- No package source changes relative to 0.5.0.

## [0.5.0](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.4.2...v0.5.0)

- **Breaking:** Make `safe` a keyword argument consistently across model
  constructors. Pass it as `safe=false` rather than a positional Boolean.
- Update `convert` to forward the `safe` keyword correctly.
- Require Julia 1.0 or later; add compatibility with BioSymbols 5 and StaticArrays 1.
- Move testing from Travis CI and AppVeyor to GitHub Actions, and repair
  documentation builds and constructor examples.

## [0.4.2](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.4.1...v0.4.2)

- Correct Julia compatibility bounds to include Julia 1.x alongside Julia 0.7.
- Add Zenodo metadata and update CI to test Julia 1.5.
- No package source changes relative to 0.4.1.

## [0.4.1](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.4.0...v0.4.1)

- Declare compatibility with BioSymbols 3 and 4, and StaticArrays 0.12.
- Expand constructor documentation and tests, including construction with
  validation disabled.
- Add CompatHelper, TagBot, and GitHub Actions documentation deployment.
- No package source changes relative to 0.4.0.

## [0.4.0](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.3.0...v0.4.0)

- Add constructors accepting parameter and base-frequency arrays, with parameter
  counts selecting absolute or relative forms through model-family constructors.
- Add `convert` methods for model construction.
- Add optional constructor validation controls and checks for array lengths.
- Reorganize constructor implementations and expand tests across model families.
- Declare `LinearAlgebra` as a runtime dependency.

## [0.3.0](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.2.3...v0.3.0)

- Replace `REQUIRE` with `Project.toml`, including dependency and test metadata.
- Add Julia 1.1 to CI.
- No package source changes relative to 0.2.3.

## [0.2.3](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.2.2...v0.2.3)

- Add `Q(model, true)` and `P(model, t, true)` to scale rates to one expected
  substitution per site per unit time.
- Correct K80 transition-probability formulas to agree with their rate matrices.
- Expand tests comparing transition probabilities with matrix exponentials,
  including scaled calculations.

## [0.2.2](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.2.1...v0.2.2)

- Correct base-frequency terms in the F84, HKY85, and TN93 transition-probability
  formulas, for both absolute and relative forms.
- Correct the displayed second transition-rate parameter in `TN93abs`.
- Add regression tests checking that transition-probability rows sum to one.

## [0.2.1](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.2.0...v0.2.1)

- Fix the sign of off-diagonal JC69 transition probabilities in both model forms.
- Repair generic transition calculations for arrays of times: use Julia's `eigen`
  interface, the inverse eigenvector matrix, and a reversible-model similarity
  transform, with a direct matrix-exponential fallback.
- Add tests for the array-of-times calculation and pin the documentation builder
  used in CI.

## [0.2.0](https://github.com/BioJulia/SubstitutionModels.jl/compare/v0.1.0...v0.2.0)

- Update the package for Julia 0.7 and 1.0, replacing older linear algebra APIs.
- Extend nucleotide indexing from static matrices to arbitrary arrays, including
  single-index access and mutation. This interface is deprecated in 0.5.1 and
  removed in 0.6.0.
- Reorganize the manual, fix matrix documentation rendering, and add contribution
  guidelines and issue/PR templates.

## [0.1.0](https://github.com/BioJulia/SubstitutionModels.jl/tree/v0.1.0)

- Initial release for Julia 0.6 with JC69, K80, F81, F84, HKY85, TN93, and GTR
  nucleotide substitution models in absolute and relative rate forms.
- Provide rate matrices through `Q` and transition-probability matrices through
  `P`, using static 4×4 matrices and model-specific formulas where available.
- Support nucleotide indexing of static matrices, model display, documentation,
  and tests.
