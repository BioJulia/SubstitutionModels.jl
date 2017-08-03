# SubstitutionModels.jl Documentation

```@docs
NucleicAcidSubstitutionModel
```

```@docs
Q
```

```@docs
P
```

## Included models

Absolute and relative rate forms of the following popular substitution models are
currently included in SubstitutionModels.jl:
* `JC69`
* `K80`
* `F81`
* `F84`
* `HKY85`
* `TN93`
* `GTR`

## Custom substitution models
The set of substitution models included in this package is easily extended with
user defined types. When defining a new substitution model, the minimum
requirement is that it is a subtype of `NucleicAcidSubstitutionModel`, and that
it has a valid method for the `Q` function.

There are two means of accomplishing this. The first method involves
recognizing that most substitution models are special cases of the Generalized
Time Reversible (GTR) model, which has a Q matrix of the form:

```math
Q = \begin{bmatrix}
-(\delta \pi_{\text{C}} + \eta \pi_{\text{G}} + \beta \pi_{\text{T}}) & \delta \pi_{\text{C}} & \eta \pi_{\text{G}} & \beta \pi_{\text{T}} \\
\delta \pi_{\text{A}} & -(\delta \pi_{\text{A}} + \epsilon \pi_{\text{G}} + \alpha \pi_{\text{T}}) & \epsilon \pi_{\text{G}} & \alpha \pi_{\text{T}} \\
\eta \pi_{\text{A}} & \epsilon \pi_{\text{C}} & -(\eta \pi_{\text{A}} + \epsilon \pi_{\text{C}} + \gamma \pi_{\text{T}}) & \gamma \pi_{\text{T}} \\
\beta \pi_{\text{A}} & \alpha \pi_{\text{C}} & \gamma \pi_{\text{G}} & -(\beta \pi_{\text{A}} + \alpha \pi_{\text{C}} + \gamma \pi_{\text{G}})
\end{bmatrix}.
```

Seeing this, a substitution model can be described by defining methods for the
following internal functions:
* `_α`,
* `_β`,
* `_γ`,
* `_δ`,
* `_ϵ`, and
* `_η`.

If this substitution model allows for unequal base frequencies, methods for
`_πA`, `_πC`, `_πG`, and `_πT` will also need to be defined. With these,
SubstitutionModels.jl will calculate the correct Q and P matrices.

The second method of describing a new substitution model's Q matrix is to do so
directly by defining a new method for the `Q` function.

### P matrix calculation for user defined substitution models
P matrices are found as follows:

```math
P = \text{expm} \left(Q \times t \right).
```

For many substitution models, the P matrix has a known form. Rather than
calculating an approximate P matrix using julia's `expm` function, an exact
method for the `P` function may be optionally defined. This is highly
recommended when it is possible, as it typically results in substantial
computational savings.
