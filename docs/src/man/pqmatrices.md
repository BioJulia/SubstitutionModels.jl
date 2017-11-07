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
