module AdaptivePredicates

include("macros.jl")
include("arithmetic.jl")
include("predicates.jl")

export orient2

function find_epsilon()
    every_other = true
    half = 0.5
    epsilon = 1.0
    splitter = 1.0
    check = 1.0
    # /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    # /*   one without causing roundoff.  (Also check if the sum is equal to   */
    # /*   the previous sum, for machines that round up instead of using exact */
    # /*   rounding.  Not that this library will work on such machines anyway. */
    cond = true
    while cond
        lastcheck = check
        epsilon *= half
        if (every_other)
            splitter *= 2.0
        end
        every_other = !every_other
        check = 1.0 + epsilon
        cond = ((check != 1.0) && (check != lastcheck))
    end
    splitter += 1.0
    return epsilon, splitter
end

const epsilon, splitter = find_epsilon()
@assert epsilon == eps(1.0) / 2
const resulterrbound = (3.0 + 8.0 * epsilon) * epsilon
const ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon
const ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon
const ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon
const o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon
const o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon
const o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon
const iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon
const iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon
const iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon
const isperrboundA = (16.0 + 224.0 * epsilon) * epsilon
const isperrboundB = (5.0 + 72.0 * epsilon) * epsilon
const isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon

end # module AdaptivePredicates
