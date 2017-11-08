var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SubstitutionModels.jl-1",
    "page": "Home",
    "title": "SubstitutionModels.jl",
    "category": "section",
    "text": "Release status:(Image: Latest Release) (Image: SubstitutionModels) (Image: License) (Image: ) (Image: BioJulia maintainer: jangevaare) (Image: BioJulia maintainer: Ward9250)Development status:(Image: Build Status) (Image: codecov.io) (Image: )"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "SubstitutionModels.jl provides facilities to model the substitution process of biological sequences. Such models are essential for the analysis of sequence evolution, phylogenetics, and simulation.We first aim to provide the most common substitution models used in the literature, but aim to build an extendable framework using julia's type system and traits, so as custom model types can be created and used."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The current release version can be installed from the Julia REPL:julia> Pkg.add(\"SubstitutionModels\")"
},

{
    "location": "index.html#Contributing-and-Questions-1",
    "page": "Home",
    "title": "Contributing and Questions",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features. Please go to the contributing section of the documentation for more information.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site."
},

{
    "location": "man/models.html#",
    "page": "Substitution models",
    "title": "Substitution models",
    "category": "page",
    "text": ""
},

{
    "location": "man/models.html#SubstitutionModels.jl-1",
    "page": "Substitution models",
    "title": "SubstitutionModels.jl",
    "category": "section",
    "text": ""
},

{
    "location": "man/models.html#What-are-subsitution-models?-1",
    "page": "Substitution models",
    "title": "What are subsitution models?",
    "category": "section",
    "text": "Substitution models are phenomenological descriptions of the evolution of DNA, as a string of four discrete states. A, T (U in RNA), C, and G.These models are Markovian, do not explicitly depict the mechanisms such as mutation or natural selection. Instead they describe the relative rates of different changes, and those rates are assumed the capture the action of those mechanisms.For example, when examining aligned sequences, a high number of transition substitutions are observed vs. the number of transversion substitutions. The mechanisms that cause this differential include mutation biases, and the action of purifying selection. However, the Kimura (K80) substitution model merely attempts to capture the effect of both of those mechanisms, rather than modelling the mechanisms themselves, by using a parameter that reflects the relative rate of transitions to transversions."
},

{
    "location": "man/models.html#SubstitutionModels.NucleicAcidSubstitutionModel",
    "page": "Substitution models",
    "title": "SubstitutionModels.NucleicAcidSubstitutionModel",
    "category": "Type",
    "text": "NucleicAcidSubstitutionModel is an abstract type that contains all models describing a substitution process impacting biological sequences of DNA or RNA with continous time Markov models.\n\n\n\n"
},

{
    "location": "man/models.html#Representation-of-models-in-code:-1",
    "page": "Substitution models",
    "title": "Representation of models in code:",
    "category": "section",
    "text": "NucleicAcidSubstitutionModel"
},

{
    "location": "man/pqmatrices.html#",
    "page": "Q and P matrices",
    "title": "Q and P matrices",
    "category": "page",
    "text": ""
},

{
    "location": "man/pqmatrices.html#SubstitutionModels.Q",
    "page": "Q and P matrices",
    "title": "SubstitutionModels.Q",
    "category": "Function",
    "text": "Generate a Q matrix for a NucleicAcidSubstitutionModel, of the form:\n\nQ = beginbmatrix\nQ_A A  Q_A C  Q_A G  Q_A T \nQ_C A  Q_C C  Q_C G  Q_C T \nQ_G A  Q_G C  Q_G G  Q_G T \nQ_T A  Q_T C  Q_T G  Q_T T endbmatrix\n\n\n\n"
},

{
    "location": "man/pqmatrices.html#SubstitutionModels.P",
    "page": "Q and P matrices",
    "title": "SubstitutionModels.P",
    "category": "Function",
    "text": "Generate a P matrix for a NucleicAcidSubstitutionModel, of the form:\n\nP = beginbmatrix\nP_A A  P_A C  P_A G  P_A T \nP_C A  P_C C  P_C G  P_C T \nP_G A  P_G C  P_G G  P_G T \nP_T A  P_T C  P_T G  P_T T endbmatrix\n\nfor specified time\n\n\n\n"
},

{
    "location": "man/pqmatrices.html#Q-and-P-Matrices-1",
    "page": "Q and P matrices",
    "title": "Q and P Matrices",
    "category": "section",
    "text": "Evolutionary analyses of sequences are conducted on a wide variety of time scales.Thus, it is convenient to express these models in terms of the instantaneous rates of change between different states. This representation of the model is typically called the model's Q Matrix.QIf we are given a starting state at one position in a DNA sequence, the model's Q matrix and a branch length expressing the expected number of changes to have occurred since the ancestor, then we can derive the probability of the descendant sequence having each of the four states.This transformation from the instantaneous rate matrix (Q Matrix), to a probability matrix for a given time period (P Matrix), is described here.P"
},

{
    "location": "man/modeltypes.html#",
    "page": "Provided models & custom models",
    "title": "Provided models & custom models",
    "category": "page",
    "text": ""
},

{
    "location": "man/modeltypes.html#Provided-and-custom-models-1",
    "page": "Provided models & custom models",
    "title": "Provided and custom models",
    "category": "section",
    "text": "Absolute and relative rate forms of the following popular substitution models are currently included in SubstitutionModels.jl:JC69\nK80\nF81\nF84\nHKY85\nTN93\nGTR"
},

{
    "location": "man/modeltypes.html#Custom-substitution-models-1",
    "page": "Provided models & custom models",
    "title": "Custom substitution models",
    "category": "section",
    "text": "The set of substitution models included in this package is easily extended with user defined types. When defining a new substitution model, the minimum requirement is that it is a subtype of NucleicAcidSubstitutionModel, and that it has a valid method for the Q function.There are two means of accomplishing this. The first method involves recognizing that most substitution models are special cases of the Generalized Time Reversible (GTR) model, which has a Q matrix of the form:Q = beginbmatrix\n-(delta pi_textC + eta pi_textG + beta pi_textT)  delta pi_textC  eta pi_textG  beta pi_textT \ndelta pi_textA  -(delta pi_textA + epsilon pi_textG + alpha pi_textT)  epsilon pi_textG  alpha pi_textT \neta pi_textA  epsilon pi_textC  -(eta pi_textA + epsilon pi_textC + gamma pi_textT)  gamma pi_textT \nbeta pi_textA  alpha pi_textC  gamma pi_textG  -(beta pi_textA + alpha pi_textC + gamma pi_textG)\nendbmatrixSeeing this, a substitution model can be described by defining methods for the following internal functions:_α,\n_β,\n_γ,\n_δ,\n_ϵ, and\n_η.If this substitution model allows for unequal base frequencies, methods for _πA, _πC, _πG, and _πT will also need to be defined. With these, SubstitutionModels.jl will calculate the correct Q and P matrices.The second method of describing a new substitution model's Q matrix is to do so directly by defining a new method for the Q function."
},

{
    "location": "man/modeltypes.html#pcomp-1",
    "page": "Provided models & custom models",
    "title": "P matrix calculation for user defined substitution models",
    "category": "section",
    "text": "P matrices are calculated as follows through matrix exponentiation:P = textexpm left(Q times t right)SubstitutionModels.jl provides generic methods that perform this exponentiation for any given substitution model (either defined in this package, or user defined).However, for many well known and well defined substitution models, the P matrix has a known form.If that is the case, whilst calculating an approximate P matrix using the generic method, will suffice, it is advised to overload the P function with an exact method for the P matrix of the model.Several provided models such as the Jukes and Cantor 69 model have their own exact P method.Providing such an exact method typically results in substantial computational savings, as the computational effort of matrix exponentiation is spared."
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site.Detailed guidance for contributing to all BioJulia packages is provided at the BioJulia Contribution Documentation.Here we list specific details about contributing and maintainership pertaining specifically to the SubstitutionModels.jl package."
},

{
    "location": "contributing.html#Named-maintainers-1",
    "page": "Contributing",
    "title": "Named maintainers",
    "category": "section",
    "text": "The named maintainers of this package are Justin Angevarre, & Ben J. Ward. It is their responsibility to make final choices about pull requests and issues, although because of our community structure, you will find other maintainers assisting them."
},

{
    "location": "contributing.html#Branching-model-1",
    "page": "Contributing",
    "title": "Branching model",
    "category": "section",
    "text": "The branching model used to develop and make releases of this package is that which is summarized in the BioJulia Contribution Documentation"
},

]}
