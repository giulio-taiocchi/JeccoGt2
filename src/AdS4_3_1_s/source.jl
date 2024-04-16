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

Sz(t, x, y, GS ::GaussianSource) = 1 + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_x(t, x, y, GS ::GaussianSource) = -((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x*(1 - tanh(12 - t)))/(GS.sigmax^3*GS.sigmay))
Sz_y(t, x, y, GS ::GaussianSource) = -((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay^3))

Sz_xx(t, x, y, GS ::GaussianSource) = ((-(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmax^2) + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^2)/GS.sigmax^4)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_yy(t, x, y, GS ::GaussianSource) = ((-(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^4)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_xxx(t, x, y, GS ::GaussianSource) = (((3*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x)/GS.sigmax^4 - (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^3)/GS.sigmax^6)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_yyy(t, x, y, GS ::GaussianSource) = (((3*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y)/GS.sigmay^4 - (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^3)/GS.sigmay^6)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_xxxx(t, x, y, GS ::GaussianSource) = (((3*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2)))/GS.sigmax^4 - (6*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^2)/GS.sigmax^6 + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^4)/GS.sigmax^8)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_yyyy(t, x, y, GS ::GaussianSource) = (((3*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2)))/GS.sigmay^4 - (6*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^6 + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^4)/GS.sigmay^8)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_xy(t, x, y, GS ::GaussianSource) = (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x*y*(1 - tanh(12 - t)))/(GS.sigmax^3*GS.sigmay^3)
Sz_xxy(t, x, y, GS ::GaussianSource) = (((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y)/(GS.sigmax^2*GS.sigmay^2) - (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^2*y)/(GS.sigmax^4*GS.sigmay^2))*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)
Sz_xyy(t, x, y, GS ::GaussianSource) = (x*(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmay^2 - (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^4)*(1 - tanh(12 - t)))/(GS.sigmax^3*GS.sigmay)
Sz_xxyy(t, x, y, GS ::GaussianSource) = (((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmay^2 - (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^4)/GS.sigmax^2 + (x^2*(-(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^4))/GS.sigmax^4)*(1 - tanh(12 - t)))/(GS.sigmax*GS.sigmay)

Sz_t(t, x, y, GS ::GaussianSource) = (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*sech(12 - t)^2)/(GS.sigmax*GS.sigmay)
Sz_tx(t, x, y, GS ::GaussianSource) = -((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x*sech(12 - t)^2)/(GS.sigmax^3*GS.sigmay))
Sz_txx(t, x, y, GS ::GaussianSource) = ((-(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmax^2) + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x^2)/GS.sigmax^4)*sech(12 - t)^2)/(GS.sigmax*GS.sigmay)
Sz_ty(t, x, y, GS ::GaussianSource) = -((exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y*sech(12 - t)^2)/(GS.sigmax*GS.sigmay^3))
Sz_tyy(t, x, y, GS ::GaussianSource) = ((-(exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*y^2)/GS.sigmay^4)*sech(12 - t)^2)/(GS.sigmax*GS.sigmay)
Sz_txy(t, x, y, GS ::GaussianSource) = (exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*x*y*sech(12 - t)^2)/(GS.sigmax^3*GS.sigmay^3)
Sz_tt(t, x, y, GS ::GaussianSource) = (2*exp(-1/2*x^2/GS.sigmax^2 - y^2/(2*GS.sigmay^2))*sech(12 - t)^2*tanh(12 - t))/(GS.sigmax*GS.sigmay)




@inline getSourcetime(ff::Source) = ff.time


