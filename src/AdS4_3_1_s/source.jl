# no source
Base.@kwdef mutable struct NoSource{T} <: Source
	time :: T = 0.0
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


#Gaussian source
Base.@kwdef mutable struct GaussianSource{T} <: Source
	time :: T = 0.0
end

Sz(t, x, y, ::GaussianSource) = exp(-1/2*x^2 - y^2/2)*(2 - tanh(12 - t))
Sz_x(t, x, y, ::GaussianSource) = -(exp(-1/2*x^2 - y^2/2)*x*(2 - tanh(12 - t)))
Sz_y(t, x, y, ::GaussianSource) = -(exp(-1/2*x^2 - y^2/2)*y*(2 - tanh(12 - t)))

Sz_xx(t, x, y, ::GaussianSource) = (-exp(-1/2*x^2 - y^2/2) + exp(-1/2*x^2 - y^2/2)*x^2)*(2 - tanh(12 - t))
Sz_yy(t, x, y, ::GaussianSource) = (-exp(-1/2*x^2 - y^2/2) + exp(-1/2*x^2 - y^2/2)*y^2)*(2 - tanh(12 - t))
Sz_xxx(t, x, y, ::GaussianSource) = (3*exp(-1/2*x^2 - y^2/2)*x - exp(-1/2*x^2 - y^2/2)*x^3)*(2 - tanh(12 - t))
Sz_yyy(t, x, y, ::GaussianSource) = (3*exp(-1/2*x^2 - y^2/2)*y - exp(-1/2*x^2 - y^2/2)*y^3)*(2 - tanh(12 - t))
Sz_xxxx(t, x, y, ::GaussianSource) = (3*exp(-1/2*x^2 - y^2/2) - 6*exp(-1/2*x^2 - y^2/2)*x^2 + exp(-1/2*x^2 - y^2/2)*x^4)*(2 - tanh(12 - t))
Sz_yyyy(t, x, y, ::GaussianSource) = (3*exp(-1/2*x^2 - y^2/2) - 6*exp(-1/2*x^2 - y^2/2)*y^2 + exp(-1/2*x^2 - y^2/2)*y^4)*(2 - tanh(12 - t))
Sz_xy(t, x, y, ::GaussianSource) = exp(-1/2*x^2 - y^2/2)*x*y*(2 - tanh(12 - t))
Sz_xxy(t, x, y, ::GaussianSource) = (exp(-1/2*x^2 - y^2/2)*y - exp(-1/2*x^2 - y^2/2)*x^2*y)*(2 - tanh(12 - t))
Sz_xyy(t, x, y, ::GaussianSource) = x*(exp(-1/2*x^2 - y^2/2) - exp(-1/2*x^2 - y^2/2)*y^2)*(2 - tanh(12 - t))
Sz_xxyy(t, x, y, ::GaussianSource) = (exp(-1/2*x^2 - y^2/2) - exp(-1/2*x^2 - y^2/2)*y^2 + x^2*(-exp(-1/2*x^2 - y^2/2) + exp(-1/2*x^2 - y^2/2)*y^2))*(2 - tanh(12 - t))

Sz_t(t, x, y, ::GaussianSource) = exp(-1/2*x^2 - y^2/2)*sech(12 - t)^2
Sz_tx(t, x, y, ::GaussianSource) = -(exp(-1/2*x^2 - y^2/2)*x*sech(12 - t)^2)
Sz_txx(t, x, y, ::GaussianSource) = (-exp(-1/2*x^2 - y^2/2) + exp(-1/2*x^2 - y^2/2)*x^2)*sech(12 - t)^2
Sz_ty(t, x, y, ::GaussianSource) = -(exp(-1/2*x^2 - y^2/2)*y*sech(12 - t)^2)
Sz_tyy(t, x, y, ::GaussianSource) = (-exp(-1/2*x^2 - y^2/2) + exp(-1/2*x^2 - y^2/2)*y^2)*sech(12 - t)^2
Sz_txy(t, x, y, ::GaussianSource) = exp(-1/2*x^2 - y^2/2)*x*y*sech(12 - t)^2

Sz_tt(t, x, y, ::GaussianSource) = 2*exp(-1/2*x^2 - y^2/2)*sech(12 - t)^2*tanh(12 - t)





@inline getSourcetime(ff::Source) = ff.time


