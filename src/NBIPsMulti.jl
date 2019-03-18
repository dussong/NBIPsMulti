"""
# `NBIPsMulti.jl`

Package for specifying interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...) with multispecies.

See `NBodyIPs` for the single species case

See `IPFitting` for the associated fitting and testing framework.
"""
module NBIPsMulti

using Reexport

include("types_and_wrappers.jl")

include("eval_nbody_multi.jl")

include("invariants_multi.jl")

include("descriptors_multi.jl")

include("environ_multi.jl")
@reexport using NBIPsMulti.EnvIPsmulti

end # module
