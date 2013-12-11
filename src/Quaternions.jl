module Quaternions

export Quaternion, LipschitzInteger, LipschitzInt, ispure, qrotation, rotationmatrix

import Base: convert, promote_rule, show, real, imag, conj, abs, abs2, inv, norm,
             +, -, *, /, ^, log, exp

immutable Quaternion{T<:Real} <: Number
    q0::T
    q1::T
    q2::T
    q3::T
end

typealias LipschitzInteger{T<:Integer} Quaternion{T}
typealias LipschitzInt LipschitzInteger{Int}

Quaternion{T<:Real}(x::T) =
    Quaternion(x, convert(T,0), convert(T,0), convert(T,0))

Quaternion{T<:Real}(z::Complex{T}) =
    Quaternion(real(z), imag(z), convert(T,0), convert(T,0))

# Conversions to quaternions
convert{T}(::Type{Quaternion{T}}, x::Real) =
    Quaternion(convert(T,x), convert(T,0), convert(T,0), convert(T,0))
convert{T}(::Type{Quaternion{T}}, z::Complex) =
    Quaternion(convert(T,real(z)), convert(T,imag(z)), convert(T,0), convert(T,0))

convert{T}(::Type{Quaternion{T}}, z::Quaternion) =
    Quaternion(convert(T,z.q0), convert(T,z.q1),
               convert(T,z.q2), convert(T,z.q3))

# Conversions from quaternions
# 2 x 2 compĺex matrix representation
convert{T}(::Type{Matrix{Complex{T}}}, z::Quaternion) =
    [z.q0 + z.q1*im z.q2 + z.q3*im ;
    -z.q2 + z.q3*im z.q0 - z.q1*im ]

convert{T}(::Type{Matrix{T}}, z::Quaternion) =
    [z.q0  z.q1  z.q2  z.q3 ;
    -z.q1  z.q0 -z.q3  z.q2 ;
    -z.q2  z.q3  z.q0 -z.q1 ;
    -z.q3 -z.q2  z.q1  z.q0 ]

promote_rule{T,S}(::Type{Complex{T}}, ::Type{Quaternion{S}}) =
    Quaternion{promote_type(T,S)}
promote_rule{S}(::Type{Bool}, ::Type{Quaternion{S}}) = Quaternion{S}
promote_rule{T<:Real,S}(::Type{T}, ::Type{Quaternion{S}}) =
    Quaternion{promote_type(T,S)}

function show(io::IO, z::Quaternion)
    show(io, z.q0)
    i = z.q1
    if sign(i) == -1
        i = -i
        print(io, " - ")
    else
        print(io, " + ")
    end
    show(io, i)
    print(io, "i")
    j = z.q2
    if sign(j) == -1
        j = -j
        print(io, " - ")
    else
        print(io, " + ")
    end
    show(io, j)
    print(io, "j")
    k = z.q3
    if sign(k) == -1
        k = -k
        print(io, " - ")
    else
        print(io, " + ")
    end
    show(io, k)
    print(io, "k")
end

real(z::Quaternion) = z.q0
imag{T}(z::Quaternion{T}) = Quaternion(zero(T), z.q1, z.q2, z.q3)

conj(z::Quaternion) = Quaternion(z.q0, -z.q1, -z.q2, -z.q3)
abs(z::Quaternion) = sqrt(z.q0*z.q0 + z.q1*z.q1 + z.q2*z.q2 + z.q3*z.q3)
abs2(z::Quaternion) = z.q0*z.q0 + z.q1*z.q1 + z.q2*z.q2 + z.q3*z.q3
inv(z::Quaternion) = conj(z)/abs2(z)
ispure(z::Quaternion) = z.q0 == 0 && (z.q1 != 0 || z.q2 != 0 || z.q3 != 0)

(-)(z::Quaternion) = Quaternion(-z.q0, -z.q1, -z.q2, -z.q3)

(/)(z::Quaternion, x::Real) = Quaternion(z.q0/x, z.q1/x, z.q2/x, z.q3/x)

(+)(z::Quaternion, w::Quaternion) = Quaternion(z.q0 + w.q0, z.q1 + w.q1,
                                               z.q2 + w.q2, z.q3 + w.q3)
(-)(z::Quaternion, w::Quaternion) = Quaternion(z.q0 - w.q0, z.q1 - w.q1,
                                               z.q2 - w.q2, z.q3 - w.q3)
(*)(z::Quaternion, w::Quaternion) = Quaternion(z.q0*w.q0 - z.q1*w.q1 - z.q2*w.q2 - z.q3*w.q3,
                                               z.q0*w.q1 + z.q1*w.q0 + z.q2*w.q3 - z.q3*w.q2,
                                               z.q0*w.q2 - z.q1*w.q3 + z.q2*w.q0 + z.q3*w.q1,
                                               z.q0*w.q3 + z.q1*w.q2 - z.q2*w.q1 + z.q3*w.q0)
(/)(z::Quaternion, w::Quaternion) = z*inv(w)

function exp(z::Quaternion)
    a = sqrt(z.q1*z.q1 + z.q2*z.q2 + z.q3*z.q3)
    return exp(z.q0)*(cos(a) + imag(z)*sin(a)/a)
end

function log(z::Quaternion)
    a = sqrt(z.q1*z.q1 + z.q2*z.q2 + z.q3*z.q3)
    b = abs(z)
    return log(b) + imag(z)/a*acos(z.q0/b)
end

^{T<:Real}(z::Quaternion{T}, x::Quaternion{T}) = exp(log(z)*x)

norm(q::Quaternion) = sqrt(q.q0*q.q0 + q.q1*q.q1 + q.q2*q.q2 + q.q3*q.q3)

## Rotations

function qrotation{T<:FloatingPoint}(axis::Vector{T}, theta)
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    u = unit(axis)
    thetaT = convert(eltype(u), theta)
    s = sin(thetaT/2)
    Quaternion(cos(thetaT/2), s*u[1], s*u[2], s*u[3])
end
qrotation{T<:Real}(rotvec::Vector{T}, theta) = qrotation(float(rotvec), theta)

# Variant of the above where the norm encodes theta
function qrotation{T<:FloatingPoint}(rotvec::Vector{T})
    if length(rotvec) != 3
        error("Must be a 3-vector")
    end
    theta = norm(rotvec)
    if theta > 0
        s = sin(theta/2)/theta  # divide by theta to make rotvec a unit vector
        return Quaternion(cos(theta/2), s*rotvec[1], s*rotvec[2], s*rotvec[3])
    end
    Quaternion(one(T), zero(T), zero(T), zero(T))
end
qrotation{T<:Real}(rotvec::Vector{T}) = qrotation(float(rotvec))


rotationmatrix(q::Quaternion) = rotationmatrix_normalized(unit(q))

rotationmatrix_normalized(q::Quaternion) =
#     [1-2*q.q1*q.q1-2*q.q2*q.q2  2*(q.q0*q.q1-q.q2*q.q3)    2*(q.q0*q.q2+q.q1*q.q3);
#      2*(q.q0*q.q1+q.q2*q.q3)    1-2*q.q0*q.q0-2*q.q2*q.q2  2*(q.q1*q.q2-q.q0*q.q3);
#      2*(q.q0*q.q2-q.q1*q.q3)    2*(q.q0*q.q3+q.q1*q.q2)    1-2*q.q0*q.q0-2*q.q1*q.q1]

    [1-2*q.q2*q.q2-2*q.q3*q.q3  2*(q.q1*q.q2+q.q0*q.q3)    2*(q.q1*q.q3-q.q0*q.q2);
     2*(q.q1*q.q2-q.q0*q.q3)    1-2*q.q1*q.q1-2*q.q3*q.q3  2*(q.q2*q.q3+q.q0*q.q1);
     2*(q.q1*q.q3+q.q0*q.q2)    2*(q.q2*q.q3-q.q0*q.q1)    1-2*q.q1*q.q1-2*q.q2*q.q2]


function unit{T}(v::Vector{T})
    nv = norm(v)
    if nv > 0
        return v/nv
    end
    zeros(promote_type(T,typeof(nv)), length(v))
end

function unit{T}(q::Quaternion{T})
    nq = norm(q)
    if nq > 0
        return q/nq
    end
    TT = promote_type(T, typeof(nq))
    z = zero(TT)
    Quaternion(z,z,z,z)
end


end #module
