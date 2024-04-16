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
	sigmax :: T = 0.5
	sigmay :: T = 0.5
end

Sz(t, x, y, ::GaussianSource) = 1 + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_x(t, x, y, ::GaussianSource) = -((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x*(1 - tanh(12 - t)))/(sigmax^3*sigmay))
Sz_y(t, x, y, ::GaussianSource) = -((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y*(1 - tanh(12 - t)))/(sigmax*sigmay^3))

Sz_xx(t, x, y, ::GaussianSource) = ((-(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmax^2) + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^2)/sigmax^4)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_yy(t, x, y, ::GaussianSource) = ((-(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmay^2) + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^4)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_xxx(t, x, y, ::GaussianSource) = (((3*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x)/sigmax^4 - (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^3)/sigmax^6)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_yyy(t, x, y, ::GaussianSource) = (((3*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y)/sigmay^4 - (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^3)/sigmay^6)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_xxxx(t, x, y, ::GaussianSource) = (((3*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2)))/sigmax^4 - (6*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^2)/sigmax^6 + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^4)/sigmax^8)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_yyyy(t, x, y, ::GaussianSource) = (((3*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2)))/sigmay^4 - (6*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^6 + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^4)/sigmay^8)*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_xy(t, x, y, ::GaussianSource) = (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x*y*(1 - tanh(12 - t)))/(sigmax^3*sigmay^3)
Sz_xxy(t, x, y, ::GaussianSource) = (((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y)/(sigmax^2*sigmay^2) - (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^2*y)/(sigmax^4*sigmay^2))*(1 - tanh(12 - t)))/(sigmax*sigmay)
Sz_xyy(t, x, y, ::GaussianSource) = (x*(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmay^2 - (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^4)*(1 - tanh(12 - t)))/(sigmax^3*sigmay)
Sz_xxyy(t, x, y, ::GaussianSource) = (((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmay^2 - (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^4)/sigmax^2 + (x^2*(-(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmay^2) + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^4))/sigmax^4)*(1 - tanh(12 - t)))/(sigmax*sigmay)

Sz_t(t, x, y, ::GaussianSource) = (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*sech(12 - t)^2)/(sigmax*sigmay)
Sz_tx(t, x, y, ::GaussianSource) = -((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x*sech(12 - t)^2)/(sigmax^3*sigmay))
Sz_txx(t, x, y, ::GaussianSource) = ((-(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmax^2) + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x^2)/sigmax^4)*sech(12 - t)^2)/(sigmax*sigmay)
Sz_ty(t, x, y, ::GaussianSource) = -((exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y*sech(12 - t)^2)/(sigmax*sigmay^3))
Sz_tyy(t, x, y, ::GaussianSource) = ((-(exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))/sigmay^2) + (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*y^2)/sigmay^4)*sech(12 - t)^2)/(sigmax*sigmay)
Sz_txy(t, x, y, ::GaussianSource) = (exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*x*y*sech(12 - t)^2)/(sigmax^3*sigmay^3)
Sz_tt(t, x, y, ::GaussianSource) = (2*exp(-1/2*x^2/sigmax^2 - y^2/(2*sigmay^2))*sech(12 - t)^2*tanh(12 - t))/(sigmax*sigmay)




@inline getSourcetime(ff::Source) = ff.time


