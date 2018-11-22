module SphericalSplines

using LinearAlgebra

export InterpolationSpline, SmoothingSpline


"""
    SplineSolution

A callable type representing a spherical spline interpolation solution.

"""
struct SplineSolution
    c::Float64
    a::Array{Float64,1}
    directions::Array{Float64,2}
end



"""
    (S::SplineSolution)(direction::Vector)

Return the value of the spline solution in the given direction.

"""
function (S::SplineSolution)(direction::Vector)
    result = S.c 
    N = size(S.a)[1]
    for k = 1:N
        result -= S.a[k] * Green(direction, S.directions[k,:])
    end
    return result
end



"""
    InterpolationSpline(directions, values)

Return a `SplineSolution` function giving an interpolation spline for the given directions
and values. This spline solution goes through all the given values.

"""
function InterpolationSpline(directions, values)
    c, a = interpolation_solution(directions, values)
    return SplineSolution(c, a, directions)
end



"""
    SmoothingSpline(directions, values, weights, smoothing)

Return a `SplineSolution` function giving a smoothing spline solution. This solution does
not pass exactly through each value, and how closely it approximates them depends on the
weights and the smoothing parameter.

"""
function SmoothingSpline(directions, values, weights, smoothing)
    c, a = smoothing_solution(directions, values, weights, smoothing)
    return SplineSolution(c, a, directions)
end



"""

    Green(a::Vector, b::Vector)

Green's function for the iterated Beltrami operator, where `a` and `b` are two direction
vectors.

"""
function Green(a::Vector, b::Vector)
    d = dot(a, b)
    if d ≈ 1.0
        return 1.0 / (4π)
    elseif d ≈ -1.0
        return 1.0 / (4π) - π/24
    end
    
    l2 = log(2)
    
    if d ≈ 0.0
        return (1 - π^2/12 + l2^2 * (-0.5)) / 4π
    end
    
    l2 = log(2)

    t1 = log(1 - d) * (log(1 + d) - l2)
    t2 = Li2(0.5 * (1 - d))
    t3 = l2^2 - l2 * log(1 + d)

    return (1 - t1 - t2 - t3) / 4π
end



"""
    Li2(x::Real ; tol=0.00001)

Compute dilogarithm function to desired tolerance using simple series representation. 
TODO: is there a better way?

"""
function Li2(x::Real ; tol=0.00001)
    s = 0.0
    sp = 0.0
    k = 0
    while true
        k += 1
        sp = s
        s += x^k / k^2
        if abs(sp - s) < tol
            return s
        end
    end
end



"""
    Gmatrix(directions)

Given a list of direction vectors (each row of `directions` corresponding to a vector),
computes the G matrix used for solving the spline coefficients.

TODO: Could be optimized a little for performance if necessary.

"""
function Gmatrix(directions)
    N = size(directions)[1]
    G = zeros(Float64, (N, N))
    for n = 1:N
        for m = 1:N
            G[n,m] = Green(directions[n, :], directions[m, :])
        end
    end
    return G
end



"""
    interpolation_solution(directions, y)

Given a list of directions as an `(N,3)` array, and an `(N,)` array of values, solve the
coefficients of the spherical interpolation spline.

"""
function interpolation_solution(directions, y)
    A = ones(size(directions)[1])'
    G = Gmatrix(directions)
    return solve_coefs(A, G, y)
end



"""
    smoothing_solution(directions, y, weights, smoothing)

Given a list of directions as an `(N,3)` array, and an `(N,)` array of values, compute the
coefficients of the spherical smoothing spline.

"""
function smoothing_solution(directions, y, weights, smoothing)
    A = ones(size(directions)[1])'
    G = Gmatrix(directions) + smoothing .* Diagonal(weights.^2)
    return solve_coefs(A, G, y)
end



"""
    solve_coefs(A, G, y)

Solve for spline coefficients `c` and `a`, given the matrices `A`, `G` and the values `y`.

"""
function solve_coefs(A, G, y)
    Gi = inv(G)
    c = inv(A * Gi * A') * A * Gi * y
    a = Gi * A' * c - Gi * y
    return c, a
end


end # module
