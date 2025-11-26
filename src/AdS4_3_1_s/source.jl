# no source
Base.@kwdef mutable struct NoSource{T} <: Source
	time :: T = 0.0
	step :: Integer = 0
	
end

Sz(t, x, y, ::NoSource) = 1.0
Sz_x(t, x, y, ::NoSource) = 0.0
Sz_y(t, x, y, ::NoSource) = 0.0

Sz_xx(t, x, y, ::NoSource) = 0.0
Sz_yy(t, x, y, ::NoSource) = 0.0
Sz_xxx(t, x, y, ::NoSource) = 0.0
Sz_yyy(t, x, y, ::NoSource) = 0.0
Sz_xxxx(t, x, y, ::NoSource) = 0.0
Sz_yyyy(t, x, y, ::NoSource) = 0.0
Sz_xy(t, x, y, ::NoSource) = 0.0
Sz_xxy(t, x, y, ::NoSource) = 0.0
Sz_xyy(t, x, y, ::NoSource) = 0.0
Sz_xxyy(t, x, y, ::NoSource) = 0.0

Sz_t(t, x, y, ::NoSource) = 0.0
Sz_tx(t, x, y, ::NoSource) = 0.0
Sz_txx(t, x, y, ::NoSource) = 0.0
Sz_ty(t, x, y, ::NoSource) = 0.0
Sz_tyy(t, x, y, ::NoSource) = 0.0
Sz_txy(t, x, y, ::NoSource) = 0.0

Sz_tt(t, x, y, ::NoSource) = 0.0


# Sourced-driven turbulence: sequence of random Fourier blocks with smooth transitions.
#
# follows the paper
# "Driven black holes: from Kolmogorov scaling to turbulent wakes.pdf"
#
# Notes:
#  - All random draws happen once at construction time
#  - b == 0 interval transitions from zero to block 1 smoothly
#  - Mixed time-space derivatives include the b == 0 case
#  - Use Random.seed!(s) before construction for reproducibility, or pass seed kwarg


Base.@kwdef mutable struct RandomFourierSequence{T} <: Source
    # time (user-updated by driver)
    time::T = 0.0

    # geometry / sampling
    MM::Int          # number of blocks
    M::Int           # modes per block
    delta::T         # duration of each block interval (Δ)
    L::T             # physical domain size L (periodic domain [0,L]×[0,L])
    kradius::T       # base wavenumber magnitude |K|

    # block data: Vector of length MM, each containing vectors of length M
    C::Vector{Vector{T}} = Vector{Vector{T}}()    # amplitudes normalized per block
    kx::Vector{Vector{T}} = Vector{Vector{T}}()    # mode kx per block
    ky::Vector{Vector{T}} = Vector{Vector{T}}()   # mode ky per block
    phi::Vector{Vector{T}} = Vector{Vector{T}}()   # phase per block

    step::Int = 0    # optional bookkeeping
end

"""
Constructor:
  RandomFourierSequence(MM, M; kradius=1.0, delta=1.0, L=1.0, seed=nothing)

Generates MM independent blocks, each with M modes.
Set seed (integer) for reproducibility.
"""
function RandomFourierSequence(MM::Int, M::Int; kradius=1.0, delta=1.0, L=1.0, seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    C  = Vector{Vector{Float64}}(undef, MM)
    kx = Vector{Vector{Float64}}(undef, MM)
    ky = Vector{Vector{Float64}}(undef, MM)
    phi = Vector{Vector{Float64}}(undef, MM)

    for b in 1:MM
        Craw = rand(M)
        C[b] = Craw ./ sum(Craw) # normalize so sum_i C_i = 1

        θ = 2π .* rand(M)
        kx[b] = kradius .* cos.(θ)
        ky[b] = kradius .* sin.(θ)

        phi[b] = 2π .* rand(M)
    end

    return RandomFourierSequence(
        time = 0.0,
        MM = MM,
        M = M,
        delta = delta,
        L = L,
        kradius = kradius,
        C = C,
        kx = kx,
        ky = ky,
        phi = phi,
    )
end

# -----------------------
# low-level mode evaluator
# -----------------------
@inline function _mode_arg(two_pi_over_L, kx, ky, x, y, phi)
    return two_pi_over_L * (kx * x + ky * y) + phi
end

@inline function mode_with_phase(two_pi_over_L, kx, ky, x, y, phi, shift)
    return cos(_mode_arg(two_pi_over_L, kx, ky, x, y, phi) + shift)
end

# -----------------------
# block spatial sum F^{(b)}(x,y)
# -----------------------
@inline function spatial_block(x::Float64, y::Float64, b::Int, RS::RandomFourierSequence)
    s = 0.0
    two_pi_over_L = 2π / RS.L
    @inbounds @simd for m in 1:RS.M
        s += RS.C[b][m] * cos(two_pi_over_L * (RS.kx[b][m]*x + RS.ky[b][m]*y) + RS.phi[b][m])
    end
    return s
end

# -----------------------
# generalized derivative of a single block:
# block_dpq returns ∂_x^p ∂_y^q F^{(b)}(x,y)
# chain rule implies multiplication by (2π/L)^(p+q) * kx^p * ky^q
# and a phase shift of (π/2)*(p+q)
# -----------------------
@inline function block_dpq(x::Float64, y::Float64, b::Int, p::Int, q::Int, RS::RandomFourierSequence)
    s = 0.0
    two_pi_over_L = 2π / RS.L
    shift = (p + q) * (π/2)
    prefactor_scale = (two_pi_over_L)^(p + q)
    @inbounds @simd for m in 1:RS.M
        kxm = RS.kx[b][m]
        kym = RS.ky[b][m]
        Cm  = RS.C[b][m]
        # multiply by (kx^p)*(ky^q) and by (2π/L)^(p+q)
        s += Cm * (kxm^p) * (kym^q) * prefactor_scale *
             mode_with_phase(two_pi_over_L, kxm, kym, x, y, RS.phi[b][m], shift)
    end
    return s
end

# convenient wrappers for p,q pairs (up to order 4)
@inline block_x(   x,y,b,RS) = block_dpq(x,y,b,1,0,RS)
@inline block_xx(  x,y,b,RS) = block_dpq(x,y,b,2,0,RS)
@inline block_xxx( x,y,b,RS) = block_dpq(x,y,b,3,0,RS)
@inline block_xxxx(x,y,b,RS) = block_dpq(x,y,b,4,0,RS)

@inline block_y(   x,y,b,RS) = block_dpq(x,y,b,0,1,RS)
@inline block_yy(  x,y,b,RS) = block_dpq(x,y,b,0,2,RS)
@inline block_yyy( x,y,b,RS) = block_dpq(x,y,b,0,3,RS)
@inline block_yyyy(x,y,b,RS) = block_dpq(x,y,b,0,4,RS)

@inline block_xy(  x,y,b,RS) = block_dpq(x,y,b,1,1,RS)
@inline block_xxy( x,y,b,RS) = block_dpq(x,y,b,2,1,RS)
@inline block_xyy( x,y,b,RS) = block_dpq(x,y,b,1,2,RS)
@inline block_xxyy(x,y,b,RS) = block_dpq(x,y,b,2,2,RS)

# -----------------------
# returns (b, i1, i2, θ, w1, w2, dθdt)
# -----------------------
@inline function interp_data(t::Float64, RS::RandomFourierSequence)
    δ = RS.delta
    b = floor(Int, t / δ)           # current interval index (0-based)
    i1 = mod(b, RS.MM) + 1         # current block index (1-based)
    i2 = mod(b + 1, RS.MM) + 1     # next block index (1-based)

    θ = (π/2) * ((t - b * δ) / δ)   # angle ∈ [0, π/2]
    w1 = cos(θ)
    w2 = sin(θ)
    dθdt = (π/2) / δ

    return b, i1, i2, θ, w1, w2, dθdt
end

# -----------------------
# combined F_dpq = ∂_x^p ∂_y^q S(t,x,y)
# where S = cosθ * F^{(i1)} + sinθ * F^{(i2)},
# and for b==0 the "previous" block is zero by convention
# -----------------------
@inline function F_dpq(t::Float64, x::Float64, y::Float64, p::Int, q::Int, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)

    if b == 0
        # previous block is zero; only w2 * F^{(1)}
        return w2 * block_dpq(x, y, 1, p, q, RS)
    end

    return w1 * block_dpq(x, y, i1, p, q, RS) +
           w2 * block_dpq(x, y, i2, p, q, RS)
end

# -----------------------
# Primary Sz (value) and spatial derivatives up to 4th order
# All arguments use Float64 to avoid ambiguous dispatch in inner loops
# -----------------------
@inline function Sz(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    # return S(t,x,y) with same logic as F_dpq with p=q=0
    return F_dpq(t, x, y, 0, 0, RS)
end

@inline Sz_x(   t,x,y,RS) = F_dpq(t,x,y,1,0,RS)
@inline Sz_xx(  t,x,y,RS) = F_dpq(t,x,y,2,0,RS)
@inline Sz_xxx( t,x,y,RS) = F_dpq(t,x,y,3,0,RS)
@inline Sz_xxxx(t,x,y,RS) = F_dpq(t,x,y,4,0,RS)

@inline Sz_y(   t,x,y,RS) = F_dpq(t,x,y,0,1,RS)
@inline Sz_yy(  t,x,y,RS) = F_dpq(t,x,y,0,2,RS)
@inline Sz_yyy( t,x,y,RS) = F_dpq(t,x,y,0,3,RS)
@inline Sz_yyyy(t,x,y,RS) = F_dpq(t,x,y,0,4,RS)

@inline Sz_xy(  t,x,y,RS) = F_dpq(t,x,y,1,1,RS)
@inline Sz_xxy( t,x,y,RS) = F_dpq(t,x,y,2,1,RS)
@inline Sz_xyy( t,x,y,RS) = F_dpq(t,x,y,1,2,RS)
@inline Sz_xxyy(t,x,y,RS) = F_dpq(t,x,y,2,2,RS)

# -----------------------
# Time derivative Sz_t and mixed time-space derivatives
# (the spatial blocks themselves have no explicit t-dependence)
# Sz_t = θ'(t) * (-sinθ F_i1 + cosθ F_i2)
# if b==0: F_i1 == 0 and F_i2 == block 1
# -----------------------
@inline function Sz_t(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)

    if b == 0
        F2 = block_dpq(x, y, 1, 0, 0, RS)
        return  cos(θ) * dθdt * F2
    else
        F1 = block_dpq(x, y, i1, 0, 0, RS)
        F2 = block_dpq(x, y, i2, 0, 0, RS)
        return (-sin(θ) * dθdt) * F1 + (cos(θ) * dθdt) * F2
    end
end

@inline function Sz_tx(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)
    if b == 0
        Fx2 = block_dpq(x, y, 1, 1, 0, RS)
        return cos(θ) * dθdt * Fx2
    else
        Fx1 = block_dpq(x, y, i1, 1, 0, RS)
        Fx2 = block_dpq(x, y, i2, 1, 0, RS)
        return (-sin(θ) * dθdt) * Fx1 + (cos(θ) * dθdt) * Fx2
    end
end

@inline function Sz_ty(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)
    if b == 0
        Fy2 = block_dpq(x, y, 1, 0, 1, RS)
        return cos(θ) * dθdt * Fy2
    else
        Fy1 = block_dpq(x, y, i1, 0, 1, RS)
        Fy2 = block_dpq(x, y, i2, 0, 1, RS)
        return (-sin(θ) * dθdt) * Fy1 + (cos(θ) * dθdt) * Fy2
    end
end

@inline function Sz_txx(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)
    if b == 0
        Fxx2 = block_dpq(x, y, 1, 2, 0, RS)
        return cos(θ) * dθdt * Fxx2
    else
        Fxx1 = block_dpq(x, y, i1, 2, 0, RS)
        Fxx2 = block_dpq(x, y, i2, 2, 0, RS)
        return (-sin(θ) * dθdt) * Fxx1 + (cos(θ) * dθdt) * Fxx2
    end
end

@inline function Sz_tyy(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)
    if b == 0
        Fyy2 = block_dpq(x, y, 1, 0, 2, RS)
        return cos(θ) * dθdt * Fyy2
    else
        Fyy1 = block_dpq(x, y, i1, 0, 2, RS)
        Fyy2 = block_dpq(x, y, i2, 0, 2, RS)
        return (-sin(θ) * dθdt) * Fyy1 + (cos(θ) * dθdt) * Fyy2
    end
end

@inline function Sz_txy(t::Float64, x::Float64, y::Float64, RS::RandomFourierSequence)
    b, i1, i2, θ, w1, w2, dθdt = interp_data(t, RS)
    if b == 0
        Fxy2 = block_dpq(x, y, 1, 1, 1, RS)
        return cos(θ) * dθdt * Fxy2
    else
        Fxy1 = block_dpq(x, y, i1, 1, 1, RS)
        Fxy2 = block_dpq(x, y, i2, 1, 1, RS)
        return (-sin(θ) * dθdt) * Fxy1 + (cos(θ) * dθdt) * Fxy2
    end
end




# Quench

Base.@kwdef mutable struct Quench{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	t0 :: T = 10.0
	sigmat :: T = 1.0 
	step :: Integer = 0
	
end

Sz(t, x, y, GS ::Quench) = 1/sqrt(1 + GS.Amp/MathConstants.e^((t - GS.t0)^2/sigmat^2))
Sz_x(t, x, y, GS ::Quench) = 0
Sz_xx(t, x, y, GS ::Quench) =0
Sz_xxx(t, x, y, GS ::Quench) = 0
Sz_xxxx(t, x, y, GS ::Quench) = 0
Sz_y(t, x, y, GS ::Quench) = 0
Sz_xy(t, x, y, GS ::Quench) = 0
Sz_xxy(t, x, y, GS ::Quench) = 0
Sz_xyy(t, x, y, GS ::Quench) = 0
Sz_xxyy(t, x, y, GS ::Quench) = 0
Sz_yy(t, x, y, GS ::Quench) = 0
Sz_yyy(t, x, y, GS ::Quench) = 0
Sz_yyyy(t, x, y, GS ::Quench) = 0
Sz_t(t, x, y, GS ::Quench) = (GS.Amp*(t - GS.t0))/(MathConstants.e^((t - GS.t0)^2/sigmat^2)*sigmat^2*(1 + GS.Amp/MathConstants.e^((t - GS.t0)^2/sigmat^2))^1.5)
Sz_tt(t, x, y, GS ::Quench) = (GS.Amp*(MathConstants.e^((t - GS.t0)^2/sigmat^2)*(sigmat^2 - 2*t^2 + 4*t*GS.t0 - 2*GS.t0^2) + GS.Amp*(sigmat^2 + t^2 - 2*t*GS.t0 + GS.t0^2)))/(sigmat^4*(MathConstants.e^((t - GS.t0)^2/sigmat^2) + GS.Amp)^2*sqrt(1 + GS.Amp/MathConstants.e^((t - GS.t0)^2/sigmat^2)))
Sz_tx(t, x, y, GS ::Quench) = 0
Sz_ty(t, x, y, GS ::Quench) = 0
Sz_txx(t, x, y, GS ::Quench) = 0
Sz_tyy(t, x, y, GS ::Quench) = 0
Sz_txy(t, x, y, GS ::Quench) = 0


#Test constant source
Base.@kwdef mutable struct SpatialConstantSource{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	t0 :: T = 10.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::SpatialConstantSource) = 1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/2.
Sz_x(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xx(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xxx(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xxxx(t, x, y, GS ::SpatialConstantSource) = 0
Sz_y(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xxy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xyy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_xxyy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_yy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_yyy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_yyyy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_t(t, x, y, GS ::SpatialConstantSource) = -0.5*(GS.Amp*sech((t - GS.t0)/GS.tau)^2)/GS.tau
Sz_tt(t, x, y, GS ::SpatialConstantSource) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*tanh((t - GS.t0)/GS.tau))/GS.tau^2
Sz_tx(t, x, y, GS ::SpatialConstantSource) = 0
Sz_ty(t, x, y, GS ::SpatialConstantSource) = 0
Sz_txx(t, x, y, GS ::SpatialConstantSource) = 0
Sz_tyy(t, x, y, GS ::SpatialConstantSource) = 0
Sz_txy(t, x, y, GS ::SpatialConstantSource) = 0


# localized source in space and time.
Base.@kwdef mutable struct PertSource{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 10.0
	L :: T = 100.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::PertSource) = 1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))
Sz_x(t, x, y, GS ::PertSource) = (GS.Amp*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmax^2)
Sz_xx(t, x, y, GS ::PertSource) =-0.5*(GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^4)
Sz_xxx(t, x, y, GS ::PertSource) = (GS.Amp*(-3*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^6)
Sz_xxxx(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(3*GS.L^4*GS.sigmax^4 - 6*GS.L^2*GS.sigmax^2*(x - GS.x0)^2 + (x - GS.x0)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmax^8)
Sz_y(t, x, y, GS ::PertSource) = (GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmay^2)
Sz_xy(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(x - GS.x0)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^2*GS.sigmay^2)
Sz_xxy(t, x, y, GS ::PertSource) = (GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^4*GS.sigmay^2)
Sz_xyy(t, x, y, GS ::PertSource) = (GS.Amp*(x - GS.x0)*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^2*GS.sigmay^4)
Sz_xxyy(t, x, y, GS ::PertSource) = (GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(GS.L^2*GS.sigmay^2 - (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmax^4*GS.sigmay^4)
Sz_yy(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmay^4)
Sz_yyy(t, x, y, GS ::PertSource) = (GS.Amp*(-3*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmay^6)
Sz_yyyy(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(3*GS.L^4*GS.sigmay^4 - 6*GS.L^2*GS.sigmay^2*(y - GS.y0)^2 + (y - GS.y0)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmay^8)
Sz_t(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*sech((t - GS.t0)/GS.tau)^2)/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.tau)
Sz_tt(t, x, y, GS ::PertSource) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*tanh((t - GS.t0)/GS.tau))/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.tau^2)
Sz_tx(t, x, y, GS ::PertSource) = (GS.Amp*(x - GS.x0)*sech((t - GS.t0)/GS.tau)^2)/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmax^2*GS.tau)
Sz_ty(t, x, y, GS ::PertSource) = (GS.Amp*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2)/(2*MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmay^2*GS.tau)
Sz_txx(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*sech((t - GS.t0)/GS.tau)^2)/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^4*GS.tau)
Sz_tyy(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2)/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmay^4*GS.tau)
Sz_txy(t, x, y, GS ::PertSource) = -0.5*(GS.Amp*(x - GS.x0)*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2)/(MathConstants.e^(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^2*GS.sigmay^2*GS.tau)

#exp source
Base.@kwdef mutable struct ExpSource{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 10.0
	L :: T = 100.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::ExpSource) = 1/sqrt(1 + (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi^2*GS.sigmax^2))

Sz_x(t, x, y, GS ::ExpSource) = (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*((2*pi*x)/GS.L - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi*GS.L*GS.sigmax^4*(1 + (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi^2*GS.sigmax^2))^1.5)

Sz_xx(t, x, y, GS ::ExpSource) = (2*pi^3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(-8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(-GS.sigmax^2 + GS.x0^2)) + GS.Amp*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(2*GS.sigmax^2 + GS.x0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^4*GS.sigmax^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_xxx(t, x, y, GS ::ExpSource) = (2*pi^4*GS.Amp*(2*pi*x - GS.L*GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(-3*GS.sigmax^2 + GS.x0^2)) + GS.Amp^2*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(6*GS.sigmax^2 + GS.x0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*x^2 - 20*pi*x*GS.L*GS.x0 + GS.L^2*(3*GS.sigmax^2 + 5*GS.x0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^6*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_xxxx(t, x, y, GS ::ExpSource) = (2*MathConstants.e^((-2*((-2*pi*x)/GS.L + GS.x0)^2)/GS.sigmax^2 - (2*((-2*pi*y)/GS.L + GS.y0)^2)/GS.sigmay^2)*pi^5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(105*GS.Amp^3*((-2*pi*x)/GS.L + GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + (180*GS.Amp^2*GS.sigmax^2*(-2*pi*x + GS.L*GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.L^2 - (180*GS.Amp^2*(-2*pi*x + GS.L*GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.L^4 + 36*GS.Amp*GS.sigmax^4*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - (216*GS.Amp*GS.sigmax^2*(-2*pi*x + GS.L*GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2)/GS.L^2 + (84*GS.Amp*(-2*pi*x + GS.L*GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2)/GS.L^4 - 24*GS.sigmax^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + (48*GS.sigmax^2*(-2*pi*x + GS.L*GS.x0)^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3)/GS.L^2 - (8*(-2*pi*x + GS.L*GS.x0)^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3)/GS.L^4))/(GS.L^4*GS.sigmax^16*((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)^4.5)

Sz_y(t, x, y, GS ::ExpSource) = (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*((2*pi*y)/GS.L - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi*GS.L*GS.sigmax^2*GS.sigmay^2*(1 + (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi^2*GS.sigmax^2))^1.5)

Sz_xy(t, x, y, GS ::ExpSource) = (-2*pi^3*GS.Amp*(2*pi*x - GS.L*GS.x0)*(2*pi*y - GS.L*GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^4*GS.sigmax^2*GS.sigmay^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_xxy(t, x, y, GS ::ExpSource) = (2*pi^4*GS.Amp*(2*pi*y - GS.L*GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(-GS.sigmax^2 + GS.x0^2)) + GS.Amp^2*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(2*GS.sigmax^2 + GS.x0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*x^2 - 20*pi*x*GS.L*GS.x0 + GS.L^2*(GS.sigmax^2 + 5*GS.x0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^4*GS.sigmay^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_xyy(t, x, y, GS ::ExpSource) = (2*pi^4*GS.Amp*(2*pi*x - GS.L*GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(-GS.sigmay^2 + GS.y0^2)) + GS.Amp^2*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(2*GS.sigmay^2 + GS.y0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*y^2 - 20*pi*y*GS.L*GS.y0 + GS.L^2*(GS.sigmay^2 + 5*GS.y0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^2*GS.sigmay^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_xxyy(t, x, y, GS ::ExpSource) = (2*MathConstants.e^((-2*((-2*pi*x)/GS.L + GS.x0)^2)/GS.sigmax^2 - (2*((-2*pi*y)/GS.L + GS.y0)^2)/GS.sigmay^2)*pi^5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(105*GS.Amp^3*(-2*pi*x + GS.L*GS.x0)^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^3 + 30*GS.Amp^2*GS.L^2*GS.sigmay^2*(-2*pi*x + GS.L*GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 30*GS.Amp^2*GS.L^2*GS.sigmax^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 180*GS.Amp^2*(-2*pi*x + GS.L*GS.x0)^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 12*GS.Amp*GS.L^4*GS.sigmax^2*GS.sigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*GS.Amp*GS.L^2*GS.sigmay^2*(-2*pi*x + GS.L*GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*GS.Amp*GS.L^2*GS.sigmax^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 84*GS.Amp*(-2*pi*x + GS.L*GS.x0)^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 8*GS.L^4*GS.sigmax^2*GS.sigmay^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*GS.L^2*GS.sigmay^2*(-2*pi*x + GS.L*GS.x0)^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*GS.L^2*GS.sigmax^2*(-2*pi*y + GS.L*GS.y0)^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 8*(-2*pi*x + GS.L*GS.x0)^2*(-2*pi*y + GS.L*GS.y0)^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3))/(GS.L^8*GS.sigmax^12*GS.sigmay^4*((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)^4.5)
Sz_yy(t, x, y, GS ::ExpSource) = (-2*pi^3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(-GS.sigmay^2 + GS.y0^2)) - GS.Amp*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(2*GS.sigmay^2 + GS.y0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^4*GS.sigmay^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_yyy(t, x, y, GS ::ExpSource) = (2*pi^4*GS.Amp*(2*pi*y - GS.L*GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(-3*GS.sigmay^2 + GS.y0^2)) + GS.Amp^2*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(6*GS.sigmay^2 + GS.y0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*y^2 - 20*pi*y*GS.L*GS.y0 + GS.L^2*(3*GS.sigmay^2 + 5*GS.y0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmay^6*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_yyyy(t, x, y, GS ::ExpSource) = (2*MathConstants.e^((-2*((-2*pi*x)/GS.L + GS.x0)^2)/GS.sigmax^2 - (2*((-2*pi*y)/GS.L + GS.y0)^2)/GS.sigmay^2)*pi^5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(105*GS.Amp^3*((-2*pi*y)/GS.L + GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + (180*GS.Amp^2*GS.sigmay^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.L^2 - (180*GS.Amp^2*(-2*pi*y + GS.L*GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.L^4 + 36*GS.Amp*GS.sigmay^4*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - (216*GS.Amp*GS.sigmay^2*(-2*pi*y + GS.L*GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2)/GS.L^2 + (84*GS.Amp*(-2*pi*y + GS.L*GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2)/GS.L^4 - 24*GS.sigmay^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + (48*GS.sigmay^2*(-2*pi*y + GS.L*GS.y0)^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3)/GS.L^2 - (8*(-2*pi*y + GS.L*GS.y0)^4*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3)/GS.L^4))/(GS.L^4*GS.sigmax^8*GS.sigmay^8*((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)^4.5)
Sz_t(t, x, y, GS ::ExpSource) = -0.125*(MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*sech((t - GS.t0)/GS.tau)^2)/(pi^2*GS.sigmax^2*GS.tau*(1 + (MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(4*pi^2*GS.sigmax^2))^1.5)
Sz_tt(t, x, y, GS ::ExpSource) = (pi*GS.Amp*sech((t - GS.t0)/GS.tau)^2*(16*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2*tanh((t - GS.t0)/GS.tau) + GS.Amp*(3*sech((t - GS.t0)/GS.tau)^2 + 4*tanh((t - GS.t0)/GS.tau)*(1 + tanh((t - GS.t0)/GS.tau)))))/(2*GS.tau^2*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_tx(t, x, y, GS ::ExpSource) = (pi^2*GS.Amp*(2*pi*x - GS.L*GS.x0)*sech((t - GS.t0)/GS.tau)^2*(8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^2*GS.sigmax^2*GS.tau*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_ty(t, x, y, GS ::ExpSource) = (pi^2*GS.Amp*(2*pi*y - GS.L*GS.y0)*sech((t - GS.t0)/GS.tau)^2*(8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^2*GS.sigmay^2*GS.tau*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2))
Sz_txx(t, x, y, GS ::ExpSource) = -((pi^3*GS.Amp*sech((t - GS.t0)/GS.tau)^2*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(-GS.sigmax^2 + GS.x0^2)) + GS.Amp^2*(4*pi^2*x^2 - 4*pi*x*GS.L*GS.x0 + GS.L^2*(2*GS.sigmax^2 + GS.x0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*x^2 - 20*pi*x*GS.L*GS.x0 + GS.L^2*(GS.sigmax^2 + 5*GS.x0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^4*GS.sigmax^4*GS.tau*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)))
Sz_tyy(t, x, y, GS ::ExpSource) = -((pi^3*GS.Amp*sech((t - GS.t0)/GS.tau)^2*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(-GS.sigmay^2 + GS.y0^2)) + GS.Amp^2*(4*pi^2*y^2 - 4*pi*y*GS.L*GS.y0 + GS.L^2*(2*GS.sigmay^2 + GS.y0^2))*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - 8*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(20*pi^2*y^2 - 20*pi*y*GS.L*GS.y0 + GS.L^2*(GS.sigmay^2 + 5*GS.y0^2))*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^4*GS.sigmay^4*GS.tau*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)))
Sz_txy(t, x, y, GS ::ExpSource) = -((pi^3*GS.Amp*(2*pi*x - GS.L*GS.x0)*(2*pi*y - GS.L*GS.y0)*sech((t - GS.t0)/GS.tau)^2*(64*MathConstants.e^(((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)*pi^4*GS.sigmax^4 - 40*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.Amp*GS.sigmax^2*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))/(GS.L^4*GS.sigmax^2*GS.sigmay^2*GS.tau*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt((MathConstants.e^(-0.5*((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 - ((-2*pi*y)/GS.L + GS.y0)^2/(2*GS.sigmay^2))*(4*MathConstants.e^((((-2*pi*x)/GS.L + GS.x0)^2/GS.sigmax^2 + ((-2*pi*y)/GS.L + GS.y0)^2/GS.sigmay^2)/2)*pi^2*GS.sigmax^2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/GS.sigmax^2)))


#Gaussian source
Base.@kwdef mutable struct GaussianSource{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 20.0
	L :: T = 100.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::GaussianSource) = 1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))
Sz_x(t, x, y, GS ::GaussianSource) = (GS.Amp*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmax^2)
Sz_xx(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^4)
Sz_xxx(t, x, y, GS ::GaussianSource) = (GS.Amp*(-3*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^6)
Sz_xxxx(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(3*GS.L^4*GS.sigmax^4 - 6*GS.L^2*GS.sigmax^2*(x - GS.x0)^2 + (x - GS.x0)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmax^8)
Sz_y(t, x, y, GS ::GaussianSource) = (GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmay^2)
Sz_xy(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(x - GS.x0)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^2*GS.sigmay^2)
Sz_xxy(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^4*GS.sigmay^2)
Sz_xyy(t, x, y, GS ::GaussianSource) = (GS.Amp*(x - GS.x0)*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^2*GS.sigmay^4)
Sz_xxyy(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*(GS.L^2*GS.sigmay^2 - (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmax^4*GS.sigmay^4)
Sz_yy(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmay^4)
Sz_yyy(t, x, y, GS ::GaussianSource) = (GS.Amp*(-3*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^6*GS.sigmay^6)
Sz_yyyy(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(3*GS.L^4*GS.sigmay^4 - 6*GS.L^2*GS.sigmay^2*(y - GS.y0)^2 + (y - GS.y0)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^8*GS.sigmay^8)
Sz_t(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*sech((t - GS.t0)/GS.tau)^2)/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.tau)
Sz_tt(t, x, y, GS ::GaussianSource) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*tanh((t - GS.t0)/GS.tau))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.tau^2)
Sz_tx(t, x, y, GS ::GaussianSource) = (GS.Amp*(x - GS.x0)*sech((t - GS.t0)/GS.tau)^2)/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmax^2*GS.tau)
Sz_ty(t, x, y, GS ::GaussianSource) = (GS.Amp*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2)/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmay^2*GS.tau)
Sz_txx(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmax^2) + (x - GS.x0)^2)*sech((t - GS.t0)/GS.tau)^2)/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^4*GS.tau)
Sz_tyy(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(-(GS.L^2*GS.sigmay^2) + (y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2)/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmay^4*GS.tau)
Sz_txy(t, x, y, GS ::GaussianSource) = -0.5*(GS.Amp*(x - GS.x0)*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2)/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^4*GS.sigmax^2*GS.sigmay^2*GS.tau)


# Lorentz Cauchy distribution. Here x0=y0=0 to simplify, and sigma is already multiplied by gamma in a unique parameter Lsigma
Base.@kwdef mutable struct LC{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	Lsigmax :: T = 0.5
	Lsigmay :: T = 0.5
	t0 :: T = 10.0
	L :: T = 100.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::LC) = 1/sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))
Sz_x(t, x, y, GS ::LC) = (x*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((x^2/GS.Lsigmax + GS.Lsigmax)^2*(1 + y^2/GS.Lsigmay^2)*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))^1.5)
Sz_xx(t, x, y, GS ::LC) =(GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-3*x^4*(y^2 + GS.Lsigmay^2) - 2*x^2*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^4*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau)))))/((x^2 + GS.Lsigmax^2)^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2)
Sz_xxx(t, x, y, GS ::LC) = (3*x*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(4*x^6*(y^2 + GS.Lsigmay^2)^2 + x^4*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(4*y^2 - GS.Lsigmay^2*(-4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) - 2*x^2*GS.Lsigmax^4*(y^2 + GS.Lsigmay^2)*(2*y^2 + GS.Lsigmay^2*(2 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) - GS.Lsigmax^6*(4*y^4 + y^2*GS.Lsigmay^2*(8 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^4*(4 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)^3*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3)
Sz_xxxx(t, x, y, GS ::LC) = (3*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-90*x^4*GS.Amp^2*GS.Lsigmax^4*(x^2 + GS.Lsigmax^2)*GS.Lsigmay^4*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))^2 + 30*x^2*GS.Amp^2*GS.Lsigmax^6*(x^2 + GS.Lsigmax^2)*GS.Lsigmay^4*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))^2 - 55*x^4*GS.Amp^3*GS.Lsigmax^6*GS.Lsigmay^6*(1 + tanh((t - GS.t0)/GS.tau))^3 + 30*x^2*GS.Amp^3*GS.Lsigmax^8*GS.Lsigmay^6*(1 + tanh((t - GS.t0)/GS.tau))^3 - 72*x^2*GS.Amp*GS.Lsigmax^2*(x^2 + GS.Lsigmax^2)*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 3*GS.Amp*GS.Lsigmax^2*(x^2 + GS.Lsigmax^2)^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 64*x^4*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 48*x^2*(x^2 + GS.Lsigmax^2)*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 4*(x^2 + GS.Lsigmax^2)^2*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 144*x^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(GS.Lsigmax*(x^2 + GS.Lsigmax^2)*GS.Lsigmay*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^3*GS.Lsigmay^3*(1 + tanh((t - GS.t0)/GS.tau)))^2))/((x^2 + GS.Lsigmax^2)^8*(y^2 + GS.Lsigmay^2)^4*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))^4.5)
Sz_y(t, x, y, GS ::LC) = (y*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(y^2/GS.Lsigmay + GS.Lsigmay)^2*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))^1.5)
Sz_xy(t, x, y, GS ::LC) = -((x*y*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(2*x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(2*y^2 - GS.Lsigmay^2*(-2 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau)))))/((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2)*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2))
Sz_xxy(t, x, y, GS ::LC) = -((y*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-6*x^6*(y^2 + GS.Lsigmay^2)^2 - x^4*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(10*y^2 + GS.Lsigmay^2*(10 - 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) - 2*x^2*GS.Lsigmax^4*(y^2 + GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2*(1 - 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + GS.Lsigmax^6*(2*y^4 + y^2*GS.Lsigmay^2*(4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^4*(2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) - GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)^2*(y^2 + GS.Lsigmay^2)*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3))
Sz_xyy(t, x, y, GS ::LC) = -((x*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*x^4*(3*y^2 - GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2)^2 - x^2*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(12*y^4 + y^2*GS.Lsigmay^2*(8 - 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - GS.Lsigmay^4*(4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + GS.Lsigmax^4*(-6*y^6 + 2*y^2*GS.Lsigmay^4*(-1 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + y^4*GS.Lsigmay^2*(-10 + 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^6*(2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) - GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2)^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3))
Sz_xxyy(t, x, y, GS ::LC) = -((GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(6*x^8*(3*y^2 - GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2)^3 - 9*x^2*GS.Amp*GS.Lsigmax^6*GS.Lsigmay^2*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))*(7*y^4 - 3*y^2*GS.Lsigmay^2*(-2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - GS.Lsigmay^4*(1 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + x^6*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)^2*(48*y^4 + y^2*GS.Lsigmay^2*(32 - 69*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^4*(-16 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + GS.Lsigmax^8*(-6*y^8 + 9*y^2*GS.Amp*GS.Lsigmay^6*(1 + tanh((t - GS.t0)/GS.tau))*(1 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - GS.Lsigmay^8*(-2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))*(1 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + y^6*GS.Lsigmay^2*(-16 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 3*y^4*GS.Lsigmay^4*(-4 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + 3*GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2)) + 3*x^4*GS.Lsigmax^4*(y^2 + GS.Lsigmay^2)*(12*y^6 - 5*y^4*GS.Lsigmay^2*(-4 + 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 2*y^2*GS.Lsigmay^4*(2 - 21*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + 3*GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2) + GS.Lsigmay^6*(-4 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + 3*GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)^2*(y^2 + GS.Lsigmay^2)^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^4))
Sz_yy(t, x, y, GS ::LC) = (GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(x^2*(-3*y^4 - 2*y^2*GS.Lsigmay^2 + GS.Lsigmay^4) + GS.Lsigmax^2*(-3*y^4 - 2*y^2*GS.Lsigmay^2 + GS.Lsigmay^4*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau)))))/((y^2 + GS.Lsigmay^2)^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2)
Sz_yyy(t, x, y, GS ::LC) = (-3*y*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-4*x^4*(y^2 - GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2)^2 - x^2*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(8*(y^4 - GS.Lsigmay^4) - GS.Amp*GS.Lsigmay^2*(y^2 + 5*GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmax^4*(-4*y^6 + y^4*GS.Lsigmay^2*(-4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 2*y^2*GS.Lsigmay^4*(2 + 3*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^6*(4 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((y^2 + GS.Lsigmay^2)^3*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3)
Sz_yyyy(t, x, y, GS ::LC) = (3*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-90*y^4*GS.Amp^2*GS.Lsigmax^4*(x^2 + GS.Lsigmax^2)*GS.Lsigmay^4*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))^2 + 30*y^2*GS.Amp^2*GS.Lsigmax^4*(x^2 + GS.Lsigmax^2)*GS.Lsigmay^6*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))^2 - 55*y^4*GS.Amp^3*GS.Lsigmax^6*GS.Lsigmay^6*(1 + tanh((t - GS.t0)/GS.tau))^3 + 30*y^2*GS.Amp^3*GS.Lsigmax^6*GS.Lsigmay^8*(1 + tanh((t - GS.t0)/GS.tau))^3 - 72*y^2*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(y^2 + GS.Lsigmay^2)*(1 + tanh((t - GS.t0)/GS.tau))*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 3*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(y^2 + GS.Lsigmay^2)^2*(1 + tanh((t - GS.t0)/GS.tau))*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 64*y^4*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 48*y^2*(y^2 + GS.Lsigmay^2)*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 4*(y^2 + GS.Lsigmay^2)^2*((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 144*y^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(GS.Lsigmax*(x^2 + GS.Lsigmax^2)*GS.Lsigmay*(y^2 + GS.Lsigmay^2) + GS.Amp*GS.Lsigmax^3*GS.Lsigmay^3*(1 + tanh((t - GS.t0)/GS.tau)))^2))/((x^2 + GS.Lsigmax^2)^4*(y^2 + GS.Lsigmay^2)^8*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))^4.5)
Sz_t(t, x, y, GS ::LC) = -0.5*(GS.Amp*sech((t - GS.t0)/GS.tau)^2)/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)*GS.tau*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))^1.5)
Sz_tt(t, x, y, GS ::LC) = (GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^4*(2*(x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2)*sinh((2*(t - GS.t0))/GS.tau) + GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + 2*cosh((2*(t - GS.t0))/GS.tau) + 2*sinh((2*(t - GS.t0))/GS.tau))))/(4*GS.tau^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2)
Sz_tx(t, x, y, GS ::LC) = (x*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^2*(2*x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(2*y^2 - GS.Lsigmay^2*(-2 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau)))))/(2*(x^2 + GS.Lsigmax^2)*GS.tau*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2)
Sz_ty(t, x, y, GS ::LC) = (y*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^2*(2*x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(2*y^2 - GS.Lsigmay^2*(-2 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau)))))/(2*(y^2 + GS.Lsigmay^2)*GS.tau*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^2)
Sz_txx(t, x, y, GS ::LC) = -0.5*(GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^2*(6*x^6*(y^2 + GS.Lsigmay^2)^2 + x^4*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(10*y^2 + GS.Lsigmay^2*(10 - 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + 2*x^2*GS.Lsigmax^4*(y^2 + GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2*(1 - 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) - GS.Lsigmax^6*(2*y^4 + y^2*GS.Lsigmay^2*(4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^4*(2 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) - GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)^2*GS.tau*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3)
Sz_tyy(t, x, y, GS ::LC) = -0.5*(GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^2*(2*x^4*(3*y^2 - GS.Lsigmay^2)*(y^2 + GS.Lsigmay^2)^2 + x^2*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(12*y^4 + y^2*GS.Lsigmay^2*(8 - 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - GS.Lsigmay^4*(4 + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + GS.Lsigmax^4*(6*y^6 + y^4*GS.Lsigmay^2*(10 - 9*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 2*y^2*GS.Lsigmay^4*(-1 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^6*(-2 - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((y^2 + GS.Lsigmay^2)^2*GS.tau*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3)
Sz_txy(t, x, y, GS ::LC) = -0.5*(x*y*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*sech((t - GS.t0)/GS.tau)^2*(4*x^4*(y^2 + GS.Lsigmay^2)^2 + 2*x^2*GS.Lsigmax^2*(y^2 + GS.Lsigmay^2)*(4*y^2 + GS.Lsigmay^2*(4 - 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))) + GS.Lsigmax^4*(4*y^4 - 2*y^2*GS.Lsigmay^2*(-4 + 5*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + GS.Lsigmay^4*(4 - 10*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))))/((x^2 + GS.Lsigmax^2)*(y^2 + GS.Lsigmay^2)*GS.tau*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/((1 + x^2/GS.Lsigmax^2)*(1 + y^2/GS.Lsigmay^2)))*(x^2*(y^2 + GS.Lsigmay^2) + GS.Lsigmax^2*(y^2 + GS.Lsigmay^2*(1 + GS.Amp + GS.Amp*tanh((t - GS.t0)/GS.tau))))^3)





#Gaussian source with moving center location.
Base.@kwdef mutable struct MovingGaussianSource{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	Lsigmax :: T = 0.5
	Lsigmay :: T = 0.5
	t0 :: T = 10.0
	L :: T = 100.0
	tau :: T = 1.0
	v :: T = 10
	step :: Integer = 0
	
end

Sz(t, x, y, GS ::MovingGaussianSource) = 1/sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))
Sz_x(t, x, y, GS ::MovingGaussianSource) = (2*pi^2*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^2*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^1.5)
Sz_xx(t, x, y, GS ::MovingGaussianSource) =(2*pi^2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(GS.Lsigmax^2 - 4*pi^2*(x + t*GS.v - GS.x0)^2) + GS.Amp*(GS.Lsigmax^2 + 2*pi^2*(x + t*GS.v - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmax^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_xxx(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(30*pi^2*GS.Amp^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2 + 9*GS.Amp*GS.Lsigmax^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 36*pi^2*GS.Amp*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 6*GS.Lsigmax^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 8*pi^2*(x + t*GS.v - GS.x0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2))/(MathConstants.e^(6*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^6*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^3.5)
Sz_xxxx(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(420*pi^4*GS.Amp^3*(x + t*GS.v - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + 180*pi^2*GS.Amp^2*GS.Lsigmax^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 720*pi^4*GS.Amp^2*(x + t*GS.v - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 9*GS.Amp*GS.Lsigmax^4*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 216*pi^2*GS.Amp*GS.Lsigmax^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 336*pi^4*GS.Amp*(x + t*GS.v - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 6*GS.Lsigmax^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 48*pi^2*GS.Lsigmax^2*(x + t*GS.v - GS.x0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 32*pi^4*(x + t*GS.v - GS.x0)^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3))/(MathConstants.e^(8*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^8*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^4.5)
Sz_y(t, x, y, GS ::MovingGaussianSource) = (2*pi^2*GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmay^2*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^1.5)
Sz_xy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(x + t*GS.v - GS.x0)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(-2*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmax^2*GS.Lsigmay^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_xxy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(30*pi^2*GS.Amp^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2 + 3*GS.Amp*GS.Lsigmax^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 36*pi^2*GS.Amp*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 2*GS.Lsigmax^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 8*pi^2*(x + t*GS.v - GS.x0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2))/(MathConstants.e^(6*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^4*GS.Lsigmay^2*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^3.5)
Sz_xyy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(-2*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(GS.Lsigmay^2 - 4*pi^2*(y - GS.y0)^2) + GS.Amp^2*(GS.Lsigmay^2 + 2*pi^2*(y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Amp*(GS.Lsigmay^2 + 20*pi^2*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmax^2*GS.Lsigmay^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_xxyy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(420*pi^4*GS.Amp^3*(x + t*GS.v - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^3 + 30*pi^2*GS.Amp^2*GS.Lsigmay^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 30*pi^2*GS.Amp^2*GS.Lsigmax^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 720*pi^4*GS.Amp^2*(x + t*GS.v - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 3*GS.Amp*GS.Lsigmax^2*GS.Lsigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*pi^2*GS.Amp*GS.Lsigmay^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*pi^2*GS.Amp*GS.Lsigmax^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 336*pi^4*GS.Amp*(x + t*GS.v - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 2*GS.Lsigmax^2*GS.Lsigmay^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*pi^2*GS.Lsigmay^2*(x + t*GS.v - GS.x0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*pi^2*GS.Lsigmax^2*(y - GS.y0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 32*pi^4*(x + t*GS.v - GS.x0)^2*(y - GS.y0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3))/(MathConstants.e^(8*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^4*GS.Lsigmay^4*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^4.5)
Sz_yy(t, x, y, GS ::MovingGaussianSource) = (2*pi^2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(GS.Lsigmay^2 - 4*pi^2*(y - GS.y0)^2) + GS.Amp*(GS.Lsigmay^2 + 2*pi^2*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmay^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_yyy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(-2*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(3*GS.Lsigmay^2 - 4*pi^2*(y - GS.y0)^2) + GS.Amp^2*(3*GS.Lsigmay^2 + 2*pi^2*(y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Amp*(3*GS.Lsigmay^2 + 20*pi^2*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmay^6*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_yyyy(t, x, y, GS ::MovingGaussianSource) = (4*pi^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(420*pi^4*GS.Amp^3*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + 180*pi^2*GS.Amp^2*GS.Lsigmay^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 720*pi^4*GS.Amp^2*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 9*GS.Amp*GS.Lsigmay^4*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 216*pi^2*GS.Amp*GS.Lsigmay^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 336*pi^4*GS.Amp*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 6*GS.Lsigmay^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 48*pi^2*GS.Lsigmay^2*(y - GS.y0)^2*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 32*pi^4*(y - GS.y0)^4*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3))/(MathConstants.e^(8*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmay^8*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^4.5)
Sz_t(t, x, y, GS ::MovingGaussianSource) = (GS.Amp*(-(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2) + 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))))/(2*GS.Lsigmax^2*GS.tau*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_tt(t, x, y, GS ::MovingGaussianSource) = (GS.Amp*(3*GS.Amp*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 4*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(4*pi^2*GS.Lsigmax^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*sech((t - GS.t0)/GS.tau)^2 + GS.Lsigmax^4*sech((t - GS.t0)/GS.tau)^2*tanh((t - GS.t0)/GS.tau) + 2*pi^2*GS.Lsigmax^2*GS.tau^2*GS.v^2*(1 + tanh((t - GS.t0)/GS.tau)) - 8*pi^4*GS.tau^2*GS.v^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau)))))/(4*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^4*GS.tau^2*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^2.5)
Sz_tx(t, x, y, GS ::MovingGaussianSource) = -((pi^2*GS.Amp*(3*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) - 2*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*(x + t*GS.v - GS.x0)*sech((t - GS.t0)/GS.tau)^2 + GS.Lsigmax^2*GS.tau*GS.v*(1 + tanh((t - GS.t0)/GS.tau)) - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau)))))/(MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^4*GS.tau*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^2.5))
Sz_ty(t, x, y, GS ::MovingGaussianSource) = (pi^2*GS.Amp*(y - GS.y0)*(4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0) + GS.Lsigmax^2*(-1 + tanh((t - GS.t0)/GS.tau)))*(1 + tanh((t - GS.t0)/GS.tau))*(-2*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmax^2*GS.Lsigmay^2*GS.tau*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_txx(t, x, y, GS ::MovingGaussianSource) = (GS.Amp*(-30*pi^4*GS.Amp^2*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) - 3*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^2*GS.Amp*GS.Lsigmax^2*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) + 12*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^4*GS.Amp*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) + 24*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^4*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*(x + t*GS.v - GS.x0)*sech((t - GS.t0)/GS.tau)^2 + GS.Lsigmax^2*GS.tau*GS.v*(1 + tanh((t - GS.t0)/GS.tau)) - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))) + 2*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^2*(pi^2*GS.Lsigmax^4*sech((t - GS.t0)/GS.tau)^2 - 4*pi^4*GS.Lsigmax^2*(x + t*GS.v - GS.x0)^2*sech((t - GS.t0)/GS.tau)^2 - 12*pi^4*GS.Lsigmax^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)) + 16*pi^6*GS.tau*GS.v*(x + t*GS.v - GS.x0)^3*(1 + tanh((t - GS.t0)/GS.tau)))))/(MathConstants.e^(6*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^6*GS.tau*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^3.5)
Sz_tyy(t, x, y, GS ::MovingGaussianSource) = (pi^2*GS.Amp*(-(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2) + 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))*(-2*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*(GS.Lsigmay^2 - 4*pi^2*(y - GS.y0)^2) + GS.Amp^2*(GS.Lsigmay^2 + 2*pi^2*(y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) - MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Amp*(GS.Lsigmay^2 + 20*pi^2*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.Lsigmax^2*GS.Lsigmay^4*GS.tau*(MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))))
Sz_txy(t, x, y, GS ::MovingGaussianSource) = (-2*GS.Amp*(y - GS.y0)*(15*pi^4*GS.Amp^2*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))^2*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) - 6*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^4*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) - 6*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^2*GS.Amp*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(pi^2*GS.Lsigmax^2*sech((t - GS.t0)/GS.tau)^2 - 4*pi^4*GS.tau*GS.v*(x + t*GS.v - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))) - 6*MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^4*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))*(GS.Lsigmax^2*(x + t*GS.v - GS.x0)*sech((t - GS.t0)/GS.tau)^2 + GS.Lsigmax^2*GS.tau*GS.v*(1 + tanh((t - GS.t0)/GS.tau)) - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))) + 4*MathConstants.e^(4*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*pi^4*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^2*(GS.Lsigmax^2*(x + t*GS.v - GS.x0)*sech((t - GS.t0)/GS.tau)^2 + GS.Lsigmax^2*GS.tau*GS.v*(1 + tanh((t - GS.t0)/GS.tau)) - 4*pi^2*GS.tau*GS.v*(x + t*GS.v - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau)))))/(MathConstants.e^(6*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2))*GS.Lsigmax^4*GS.Lsigmay^2*GS.tau*(1 + (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/MathConstants.e^(2*pi^2*((x + t*GS.v - GS.x0)^2/GS.Lsigmax^2 + (y - GS.y0)^2/GS.Lsigmay^2)))^3.5)



Base.@kwdef mutable struct GaussianTilt{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 20.0
	L :: T = 100.0
	tau :: T = 1.0
	phi :: T = 0.0
	step :: Integer = 0
end

#tilted gaussian
Sz(t, x, y, GS ::GaussianTilt) = 1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2)))
Sz_x(t, x, y, GS ::GaussianTilt) = (GS.Amp*((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(4*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^2)
Sz_xx(t, x, y, GS ::GaussianTilt) = (GS.Amp*(2*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2) - ((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4)
Sz_xxx(t, x, y, GS ::GaussianTilt) = (GS.Amp*((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*(-6*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2) + ((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(16*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^6)
Sz_xxxx(t, x, y, GS ::GaussianTilt) = (GS.Amp*(-12*GS.L^4*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2)^2 + 12*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2)*((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2 - ((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(32*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^8)
Sz_y(t, x, y, GS ::GaussianTilt) = (GS.Amp*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(4*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^2)
Sz_xy(t, x, y, GS ::GaussianTilt) = (GS.Amp*(-(((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)) + 2*GS.L^2*(GS.sigmax^(-2) - GS.sigmay^(-2))*sin(2*GS.phi))*(1 + tanh((t - GS.t0)/GS.tau)))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4)
Sz_xxy(t, x, y, GS ::GaussianTilt) = (GS.Amp*(-2*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2)*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2) + ((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2) - 4*GS.L^2*(GS.sigmax^(-2) - GS.sigmay^(-2))*((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*sin(2*GS.phi))*(1 + tanh((t - GS.t0)/GS.tau)))/(16*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^6)
Sz_xyy(t, x, y, GS ::GaussianTilt) = (GS.Amp*(GS.L^2*GS.sigmax^2*GS.sigmay^2*(cos(GS.phi)^2*GS.sigmax^2 + GS.sigmay^2*sin(GS.phi)^2)*(GS.sigmax^2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) - cos(GS.phi)*GS.sigmay^2*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))) + 2*cos(GS.phi)*GS.L^2*GS.sigmax^2*GS.sigmay^2*(GS.sigmax^2 - GS.sigmay^2)*sin(GS.phi)*(cos(GS.phi)*GS.sigmax^2*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) + GS.sigmay^2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))) - (GS.sigmax^2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) - cos(GS.phi)*GS.sigmay^2*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))*(cos(GS.phi)*GS.sigmax^2*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) + GS.sigmay^2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^6*GS.sigmax^6*GS.sigmay^6)
Sz_xxyy(t, x, y, GS ::GaussianTilt) = -0.5*(GS.Amp*(GS.L^4*GS.sigmax^4*GS.sigmay^4*(cos(GS.phi)^2*GS.sigmay^2 + GS.sigmax^2*sin(GS.phi)^2)*(cos(GS.phi)^2*GS.sigmax^2 + GS.sigmay^2*sin(GS.phi)^2) - GS.L^2*GS.sigmax^2*GS.sigmay^2*(cos(GS.phi)^2*GS.sigmax^2 + GS.sigmay^2*sin(GS.phi)^2)*(GS.sigmax^2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) - cos(GS.phi)*GS.sigmay^2*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))^2 - 4*cos(GS.phi)*GS.L^2*GS.sigmax^2*GS.sigmay^2*(GS.sigmax^2 - GS.sigmay^2)*sin(GS.phi)*(GS.sigmax^2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) - cos(GS.phi)*GS.sigmay^2*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))*(cos(GS.phi)*GS.sigmax^2*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) + GS.sigmay^2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))) - GS.L^2*GS.sigmax^2*GS.sigmay^2*(cos(GS.phi)^2*GS.sigmay^2 + GS.sigmax^2*sin(GS.phi)^2)*(cos(GS.phi)*GS.sigmax^2*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) + GS.sigmay^2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))^2 + (GS.sigmax^2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) - cos(GS.phi)*GS.sigmay^2*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))^2*(cos(GS.phi)*GS.sigmax^2*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)) + GS.sigmay^2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))^2 + (GS.L^4*GS.sigmax^4*GS.sigmay^4*(GS.sigmax^2 - GS.sigmay^2)^2*sin(2*GS.phi)^2)/2)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^8*GS.sigmax^8*GS.sigmay^8)
Sz_yy(t, x, y, GS ::GaussianTilt) = (GS.Amp*(2*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmay^2 + (2*sin(GS.phi)^2)/GS.sigmax^2) - ((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4)
Sz_yyy(t, x, y, GS ::GaussianTilt) = (GS.Amp*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*(-6*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmay^2 + (2*sin(GS.phi)^2)/GS.sigmax^2) + ((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2)*(1 + tanh((t - GS.t0)/GS.tau)))/(16*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^6)
Sz_yyyy(t, x, y, GS ::GaussianTilt) = (GS.Amp*(-12*GS.L^4*((2*cos(GS.phi)^2)/GS.sigmay^2 + (2*sin(GS.phi)^2)/GS.sigmax^2)^2 + 12*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmay^2 + (2*sin(GS.phi)^2)/GS.sigmax^2)*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2 - ((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^4)*(1 + tanh((t - GS.t0)/GS.tau)))/(32*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^8)
Sz_t(t, x, y, GS ::GaussianTilt) = -0.5*(GS.Amp*sech((t - GS.t0)/GS.tau)^2)/(exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.tau)
Sz_tt(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*tanh((t - GS.t0)/GS.tau))/(exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.tau^2)
Sz_tx(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2))/(4*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^2*GS.tau)
Sz_ty(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2))/(4*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^2*GS.tau)
Sz_txx(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(2*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmax^2 + (2*sin(GS.phi)^2)/GS.sigmay^2) - ((2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 - (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4*GS.tau)
Sz_tyy(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(2*GS.L^2*((2*cos(GS.phi)^2)/GS.sigmay^2 + (2*sin(GS.phi)^2)/GS.sigmax^2) - ((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)^2))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4*GS.tau)
Sz_txy(t, x, y, GS ::GaussianTilt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(-(((-2*sin(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*cos(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)*((2*cos(GS.phi)*(cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi)))/GS.sigmay^2 + (2*sin(GS.phi)*(cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi)))/GS.sigmax^2)) + 2*GS.L^2*(GS.sigmax^(-2) - GS.sigmay^(-2))*sin(2*GS.phi)))/(8*exp(((cos(GS.phi)*(y - GS.y0) + (-x + GS.x0)*sin(GS.phi))^2/GS.sigmay^2 + (cos(GS.phi)*(x - GS.x0) + (y - GS.y0)*sin(GS.phi))^2/GS.sigmax^2)/(2*GS.L^2))*GS.L^4*GS.tau)


Base.@kwdef mutable struct Gaussiangtt{T} <: Source
	time :: T = 0.0
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 20.0
	L :: T = 100.0
	tau :: T = 1.0
	step :: Integer = 0
end

Sz(t, x, y, GS ::Gaussiangtt) = 1/sqrt(1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_x(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmax^2*(1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))^1.5)
Sz_xx(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(-4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*(GS.L^2*GS.sigmax^2 - (x - GS.x0)^2) + GS.Amp*(2*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(2*GS.L^4*GS.sigmax^4*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_xxx(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(3*GS.L^2*GS.sigmax^2 - (x - GS.x0)^2) + GS.Amp^2*(6*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(3*GS.L^2*GS.sigmax^2 + 5*(x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^6*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_xxxx(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(-360*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp^2*GS.L^2*GS.sigmax^2*(x - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2 + 360*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp^2*(x - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2 + 180*GS.Amp^3*GS.L^2*GS.sigmax^2*(x - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^3 - 75*GS.Amp^3*(x - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + 24*GS.L^4*GS.sigmax^4*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 48*GS.L^2*GS.sigmax^2*(x - GS.x0)^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*(x - GS.x0)^4*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 36*GS.Amp*GS.L^4*GS.sigmax^4*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 216*GS.Amp*GS.L^2*GS.sigmax^2*(x - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 84*GS.Amp*(x - GS.x0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2))/(8*Sqrt(2)*exp((2*((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2))/GS.L^2)*GS.L^8*GS.sigmax^8*(2 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))^4.5)
Sz_y(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau)))/(exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.L^2*GS.sigmay^2*(1 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))^1.5)
Sz_xy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(x - GS.x0)*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(2*GS.L^4*GS.sigmax^2*GS.sigmay^2*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_xxy(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(GS.L^2*GS.sigmax^2 - (x - GS.x0)^2) + GS.Amp^2*(2*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(GS.L^2*GS.sigmax^2 + 5*(x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^4*GS.sigmay^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_xyy(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(x - GS.x0)*(1 + tanh((t - GS.t0)/GS.tau))*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(GS.L^2*GS.sigmay^2 - (y - GS.y0)^2) + GS.Amp^2*(2*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(GS.L^2*GS.sigmay^2 + 5*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmax^2*GS.sigmay^4*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_xxyy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(105*GS.Amp^3*(x - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^3 + 8*GS.L^4*GS.sigmax^2*GS.sigmay^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 8*GS.L^2*GS.sigmay^2*(x - GS.x0)^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 8*GS.L^2*GS.sigmax^2*(y - GS.y0)^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*(x - GS.x0)^2*(y - GS.y0)^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 30*GS.Amp^2*GS.L^2*GS.sigmay^2*(x - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 30*GS.Amp^2*GS.L^2*GS.sigmax^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) - 180*GS.Amp^2*(x - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))) + 12*GS.Amp*GS.L^4*GS.sigmax^2*GS.sigmay^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*GS.Amp*GS.L^2*GS.sigmay^2*(x - GS.x0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 36*GS.Amp*GS.L^2*GS.sigmax^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 84*GS.Amp*(x - GS.x0)^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2))/(8*Sqrt(2)*exp((2*((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2))/GS.L^2)*GS.L^8*GS.sigmax^4*GS.sigmay^4*(2 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))^4.5)
Sz_yy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(-4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*(GS.L^2*GS.sigmay^2 - (y - GS.y0)^2) + GS.Amp*(2*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(2*GS.L^4*GS.sigmay^4*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_yyy(t, x, y, GS ::Gaussiangtt) = -0.25*(GS.Amp*(y - GS.y0)*(1 + tanh((t - GS.t0)/GS.tau))*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(3*GS.L^2*GS.sigmay^2 - (y - GS.y0)^2) + GS.Amp^2*(6*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(3*GS.L^2*GS.sigmay^2 + 5*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^6*GS.sigmay^6*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_yyyy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))*(-360*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp^2*GS.L^2*GS.sigmay^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^2 + 360*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp^2*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^2 + 180*GS.Amp^3*GS.L^2*GS.sigmay^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))^3 - 75*GS.Amp^3*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))^3 + 24*GS.L^4*GS.sigmay^4*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 - 48*GS.L^2*GS.sigmay^2*(y - GS.y0)^2*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 8*(y - GS.y0)^4*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3 + 36*GS.Amp*GS.L^4*GS.sigmay^4*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 - 216*GS.Amp*GS.L^2*GS.sigmay^2*(y - GS.y0)^2*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2 + 84*GS.Amp*(y - GS.y0)^4*(1 + tanh((t - GS.t0)/GS.tau))*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2))/(8*Sqrt(2)*exp((2*((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2))/GS.L^2)*GS.L^8*GS.sigmay^8*(2 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))^4.5)
Sz_t(t, x, y, GS ::Gaussiangtt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2)/(sqrt(2)*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.tau*(2 - (GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)))^1.5)
Sz_tt(t, x, y, GS ::Gaussiangtt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(-8*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*tanh((t - GS.t0)/GS.tau) + GS.Amp*(3*sech((t - GS.t0)/GS.tau)^2 + 4*tanh((t - GS.t0)/GS.tau)*(1 + tanh((t - GS.t0)/GS.tau)))))/(2*GS.tau^2*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_tx(t, x, y, GS ::Gaussiangtt) = -0.5*(GS.Amp*(x - GS.x0)*sech((t - GS.t0)/GS.tau)^2*(4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^2*GS.sigmax^2*GS.tau*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_ty(t, x, y, GS ::Gaussiangtt) = -0.5*(GS.Amp*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2*(4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau))))/(GS.L^2*GS.sigmay^2*GS.tau*(-2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) + GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^2*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_txx(t, x, y, GS ::Gaussiangtt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(GS.L^2*GS.sigmax^2 - (x - GS.x0)^2) + GS.Amp^2*(2*GS.L^2*GS.sigmax^2 + (x - GS.x0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(GS.L^2*GS.sigmax^2 + 5*(x - GS.x0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(4*GS.L^4*GS.sigmax^4*GS.tau*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_tyy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*sech((t - GS.t0)/GS.tau)^2*(-16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2)*(GS.L^2*GS.sigmay^2 - (y - GS.y0)^2) + GS.Amp^2*(2*GS.L^2*GS.sigmay^2 + (y - GS.y0)^2)*sech((t - GS.t0)/GS.tau)^2*(cosh((2*(t - GS.t0))/GS.tau) + sinh((2*(t - GS.t0))/GS.tau)) + 4*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(GS.L^2*GS.sigmay^2 + 5*(y - GS.y0)^2)*(1 + tanh((t - GS.t0)/GS.tau))))/(4*GS.L^4*GS.sigmay^4*GS.tau*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
Sz_txy(t, x, y, GS ::Gaussiangtt) = (GS.Amp*(x - GS.x0)*(y - GS.y0)*sech((t - GS.t0)/GS.tau)^2*(16*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/GS.L^2) + 20*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)) + GS.Amp^2*(1 + tanh((t - GS.t0)/GS.tau))^2))/(4*GS.L^4*GS.sigmax^2*GS.sigmay^2*GS.tau*(2*exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2)) - GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))^3*sqrt(4 - (2*GS.Amp*(1 + tanh((t - GS.t0)/GS.tau)))/exp(((x - GS.x0)^2/GS.sigmax^2 + (y - GS.y0)^2/GS.sigmay^2)/(2*GS.L^2))))
@inline getSourcetime(ff::Source) = ff.time


