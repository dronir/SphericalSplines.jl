module SphericalSplines

using LinearAlgebra
using GSL

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
    (S::SplineSolution)(direction)

Return the value of the interpolation spline solution in the given direction.

"""
function (S::SplineSolution)(direction::Vector)
    theta, phi = angles(direction)
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
function InterpolationSpline(directions::Array, values::Vector)
    c, a = interpolation_solution(directions, values)
    return SplineSolution(c, a, directions)
end



"""
    SmoothingSpline(directions, values, weights, smoothing)

Return a `SplineSolution` function giving a smoothing spline solution. This solution does
not pass exactly through each value, and how closely it approximates them depends on the
weights and the smoothing parameter.

"""
function SmoothingSpline(directions::AbstractArray, values::Vector, weights::Vector, smoothing::Real)
    c, a = smoothing_solution(directions, values, weights, smoothing)
    return SplineSolution(c, a, directions)
end



"""
    angles(v::Vector)

Get direction angles from vector.

"""
angles(v::Vector) = (acos(v[3]), atan(v[2], v[1]))



"""
    sphY(l::Integer, m::Integer, theta::Real, phi::Real)

Real-valued spherical harmonic basis function of degree `l`, order `m` at `(theta,
phi)`.

Computed using the GSL function `sf_legendre_sphPlm_e`.

"""
function sphY(l::Integer, m::Integer, theta::Real, phi::Real)
    ct = cos(theta)
    res = sqrt(2) * sf_legendre_sphPlm_e(l, abs(m), ct).val * (-1)^m
    return m<0 ? -res * sin(m*phi) : res * cos(m*phi)
end


"""

    fundamental_system(directions::Array{Float64,2})

Create a 'fundamental system', a matrix where each column has the spherical harmonics of


"""
function fundamental_system(directions::Array)
    N = size(directions)[1]
    if size(directions)[2] != 3
        error("Expected shape (:,3), got $(size(directions))")
    end
    M = floor(Int, sqrt(N))
    result = zeros(M^2, N)
    for n = 1:N
        theta, phi = angles(directions[n,:])
        k = 0
        for l = 0:M-1
            for m = -l:l
                k += 1
                result[k,n] = sphY(l, m, theta, phi)
            end
        end
    end
    return result
end


"""

    Green(a, b)

Green's function for the iterated Beltrami operator, where `a` and `b` are two directions.

"""
function Green(a::Vector, b::Vector)
    d = dot(a, b)
    if d ≈ 1
        return 1.0 / (4π)
    end
    if d ≈ -1
        return 1.0 / (4π) - π/24
    end
    l2 = log(2)

    t1 = log(1 - d) * (log(1 + d) - l2)
    t2 = Li2(0.5 * (1 - d))
    t3 = l2^2 - l2 * log(1 + d)

    return (1 - t1 - t2 - t3) / 4π
end



"""
    Li2(x::Real ; tol=0.0001)

Compute dilogarithm function to desired tolerance using simple series representation. 
TODO: is there a better way?

"""
function Li2(x::Real ; tol=0.0001)
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
    Gmatrix(directions::Array)

Given a list of direction vectors (each row of `directions` corresponding to a vector),
computes the G matrix used for solving the spline coefficients.

"""
function Gmatrix(directions::Array)
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
    interpolation_solution(directions::Array, y::Vector)

Given a list of directions as an `(N,3)` array, and an `(N,)` array of values, solve the
coefficients of the spherical interpolation spline.

"""
function interpolation_solution(directions, y)
    A = ones(size(directions)[1])'
    G = Gmatrix(directions)
    M = size(A)[1]

    Gi = inv(G)
    c = inv(A * Gi * A') * A * Gi * y
    a = Gi * A' * c - Gi * y
    
    return c, a
end


"""
    smoothing_solution(directions::Array, y::Vector)

Given a list of directions as an `(N,3)` array, and an `(N,)` array of values, compute the
coefficients of the spherical smoothing spline.

"""
function smoothing_solution(directions::Array, y::Vector, weights::Vector, smoothing::Real)
    A = ones(size(directions)[1])'
    B = smoothing .* Diagonal(weights.^2)
    G = Gmatrix(directions) + B
    M = size(A)[1]

    Gi = inv(G)
    c = inv(A * Gi * A') * A * Gi * y
    a = Gi * A' * c - Gi * y
    
    return c, a
end


end # module
