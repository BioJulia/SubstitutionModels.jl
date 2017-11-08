# SubstitutionModels.jl

## What are subsitution models?

Substitution models are phenomenological descriptions of the evolution of DNA,
as a string of four discrete states. A, T (U in RNA), C, and G.

These models are Markovian, do not explicitly depict the mechanisms such as
mutation or natural selection.
Instead they describe the relative rates of different changes, and those rates
are assumed the capture the action of those mechanisms.

For example, when examining aligned sequences, a high number of transition
substitutions are observed vs. the number of transversion substitutions.
The mechanisms that cause this differential include mutation biases, and the
action of purifying selection.
However, the Kimura (K80) substitution model merely attempts to capture the
*effect* of both of those mechanisms, rather than modelling the mechanisms
themselves, by using a parameter that reflects the relative rate of transitions
to transversions.

## Representation of models in code:

```@docs
NucleicAcidSubstitutionModel
```
