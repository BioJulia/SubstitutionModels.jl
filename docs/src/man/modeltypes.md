# Nucleic Acid Substitution Models

## Included Substitution Models
Absolute and relative rate forms of the following common substitution models are
currently included in SubstitutionModels.jl:
* `JC69`
* `K80`
* `F81`
* `F84`
* `HKY85`
* `TN93`
* `GTR`

## Construction
For the convenience of the user, several construction methods exist for each substitution model. Each models are represented by immutable structs; i.e. once constructed, model parameters can not be changed.

### Method 1
The most basic construction method detects whether the absolute or relative rate form of a substitution model is being referenced based on the number of parameters, and constructs it:
```jldoctest
K80(1e-2)

# output
Kimura 1980 model (relative rate form)
κ = 0.01
```

Substitution model construction involves a number of checks for parameter validity, though these may be optionally bypassed:

```jldoctest
K80(-1e-2, false)

# output
Kimura 1980 model (relative rate form)
κ = -0.01
```

The exact model type may be specific as well, to avoid the software detection of absolute or relative rate specification. The specific forms are the parent models' name followed by `abs` for absolute or `rel` for relative rate form.

```jldoctest
K80rel<:K80, K80abs<:K80

# output
(true, true)
```

```jldoctest
K80abs(1e-2, 1.5e-2)

# output
Kimura 1980 model (absolute rate form)
α = 0.01, β = 0.015
```

For models requiring specification of base frequencies, those can be provided in the same manner as the other model parameters:

```jldoctest
fieldnames(F81abs)

# output
(:β, :πA, :πC, :πG, :πT)
```

```jldoctest
F81abs(1e-2, 0.2475, 0.2425, 0.2575, 0.2525)

# output
Felsenstein 1981 model (absolute rate form)
β = 0.01, π = [0.2475, 0.2425, 0.2575, 0.2525]
```

### Method 2
It may be familar or convenient for some SubstitutionModels.jl users to construct substitution model objects using a parameter vector of some kind.
This is also a supported construction method. When specifying a model with an equal base frequency, a single parameter vector is required for construction.
For models that support different base frequencies, a parameter vector and a base frequency vector are both required.

```jldoctest
F81abs([1e-2], [0.2475, 0.2425, 0.2575, 0.2525])

# output
Felsenstein 1981 model (absolute rate form)
β = 0.01, π = [0.2475, 0.2425, 0.2575, 0.2525]
```

The unsafe construction method is still available, as is auto-detection of absolute or relative rate types of substitution models.

```jldoctest
K80([1e-2, 1.5e-2], true)

# output
Kimura 1980 model (absolute rate form)
α = 0.01, β = 0.015
```

### Method 3
Lastly, `convert` methods are also provided for each substitution model type:

```jldoctest
convert(K80, [1e-2, 1.5e-2], true)

# output
Kimura 1980 model (absolute rate form)
α = 0.01, β = 0.015
```

## User defined substitution models
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

Pull requests that provide additional substitution models are welcome.

### [P matrix calculation for user defined substitution models](@id pcomp)

P matrices are calculated as follows through matrix exponentiation:

```math
P = \text{exp} \left(Q \times t \right).
```

SubstitutionModels.jl provides generic methods that perform this exponentiation
for any given substitution model (either defined in this package, or user defined).

However, for many well known and well defined substitution models, the P matrix
has a known form.

If that is the case, whilst calculating an approximate P matrix using the
generic method, will suffice, it is advised to overload the `P` function with an
exact method for the `P` matrix of the model.

Several provided models such as the Jukes and Cantor 69 model have their own
exact `P` method.

Providing such an exact method typically results in substantial computational
savings, as the computational effort of matrix exponentiation is spared.
