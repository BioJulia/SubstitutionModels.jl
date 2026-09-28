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

Direct nucleotide indexing on arbitrary arrays was deprecated in 0.5.1 and
removed in 0.6.0. Convert each index as above, or wrap the array explicitly:

```julia
v = NucleotideView(f)
v[DNA_A] # 0.21; v shares f's storage
v[RNA_C] = 0.3 # updates f[2]
q = NucleotideView(p)
q[DNA_A, DNA_G]
q[DNA_A, :] # copy the A row
view(q, :, RNA_G) # view the G column
parent(q) === p # true
```

The wrapper accepts four-element vectors and 4×4 matrices with one-based axes.
It preserves integer indexing, iteration, and the parent's mutability. Wrapping
an immutable static matrix does not make it mutable. `copy(v)` wraps a separate
copy of the parent. Array operations such as slicing and broadcasting return
ordinary numerical arrays. A single nucleotide index on a matrix remains linear:
`q[DNA_G]` means `p[3]`.

`P` and `Q` still return static numerical matrices. Downstream packages that
restrict SubstitutionModels to 0.5 must also update their dependency bounds
when migrating to 0.6.

```@docs
nucleotide_index
NucleotideView
```
