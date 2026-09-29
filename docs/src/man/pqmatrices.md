# Q and P Matrices

Evolutionary analyses of sequences are conducted on a wide variety of time scales.

Thus, it is convenient to express these models in terms of the instantaneous
rates of change between different states. This representation of the model is
typically called the model's Q Matrix.

```@docs
Q
```

If we are given a starting state at one position in a DNA sequence, the model's
Q matrix and a branch length expressing the expected number of changes to have
occurred since the ancestor, then we can derive the probability of the
descendant sequence having each of the four states.

This transformation from the instantaneous rate matrix (Q Matrix), to a
probability matrix for a given time period (P Matrix), is described
[here](@ref pcomp).

```@docs
P
```

## Nucleotide indices

Matrix rows and columns use A/C/G/T order for DNA and A/C/G/U for RNA.
Use `nucleotide_index` to convert a single, unambiguous nucleotide to its position:

```julia
using SubstitutionModels, BioSymbols
f = [0.21, 0.29, 0.23, 0.27]
f[nucleotide_index(DNA_A)] # 0.21
p = P(JC69(), 0.1)
p[nucleotide_index(DNA_A), nucleotide_index(DNA_G)]
```

`nucleotide_index` throws `ArgumentError` for ambiguity symbols (including `N`)
and gaps. An index outside the array throws `BoundsError`, including when a
model constructor is called with `safe=false`.

In 0.5.1, direct indexing such as `f[DNA_A]` and `p[DNA_A, DNA_G]` remains
available but is deprecated. Julia displays the migration warning when
deprecation warnings are enabled, for example with `--depwarn=yes`.
These methods will be removed in 0.6.0. Convert each nucleotide index explicitly
as above. A single nucleotide index on a matrix still means a linear index:
`p[DNA_G]` means `p[3]`, not a row or column.

In 0.6.0, use the exported `NucleotideView(a)` constructor to opt into nucleotide
indexing on a four-element vector or a 4×4 matrix:

```julia
# Requires SubstitutionModels 0.6.0
v = NucleotideView(f)
v[DNA_A] # 0.21
v[RNA_C] = 0.3 # updates f[2]
q = NucleotideView(p)
q[DNA_A, DNA_G]
q[DNA_A, :] # copy the A row
view(q, :, RNA_G) # view the G column
parent(q) === p # true
```

The wrapper shares its parent's storage and requires one-based axes. Mutation
requires a mutable parent; use `copy` to obtain a wrapper over separate storage.
A single index on a wrapped matrix is linear. Ambiguities and gaps still throw
`ArgumentError`. This constructor is available starting in 0.6.0; for 0.5.1,
use `nucleotide_index` as shown above.

`P` and `Q` continue to return ordinary static numerical matrices.

```@docs
nucleotide_index
```
