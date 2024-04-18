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
	Amp :: T = 0.01
	x0 :: T = 0.0
	y0 :: T = 0.0
	sigmax :: T = 0.5
	sigmay :: T = 0.5
	t0 :: T = 10.0
end

Sz(t, x, y, GS ::GaussianSource) = 1 + (GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_x(t, x, y, GS ::GaussianSource) = -1/2*(GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)*(1 + tanh(t - GS.t0)))/(Pi*GS.sigmax^3*GS.sigmay)
Sz_y(t, x, y, GS ::GaussianSource) =-1/2*(GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)*(1 + tanh(t - GS.t0)))/(Pi*GS.sigmax*GS.sigmay^3)

Sz_xx(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmax^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^2)/GS.sigmax^4)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_yy(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^4)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_xxx(t, x, y, GS ::GaussianSource) = (GS.Amp*((3*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0))/GS.sigmax^4 - (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^3)/GS.sigmax^6)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_yyy(t, x, y, GS ::GaussianSource) = (GS.Amp*((3*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0))/GS.sigmay^4 - (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^3)/GS.sigmay^6)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_xxxx(t, x, y, GS ::GaussianSource) =(GS.Amp*((3*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2)))/GS.sigmax^4 - (6*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^2)/GS.sigmax^6 + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^4)/GS.sigmax^8)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_yyyy(t, x, y, GS ::GaussianSource) = (GS.Amp*((3*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2)))/GS.sigmay^4 - (6*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^6 + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^4)/GS.sigmay^8)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)

Sz_xy(t, x, y, GS ::GaussianSource) = (GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)*(y - GS.y0)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax^3*GS.sigmay^3)
Sz_xxy(t, x, y, GS ::GaussianSource) =(GS.Amp*((exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0))/(GS.sigmax^2*GS.sigmay^2) - (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^2*(y - GS.y0))/(GS.sigmax^4*GS.sigmay^2))*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)
Sz_xyy(t, x, y, GS ::GaussianSource) = -1/2*(GS.Amp*(x - GS.x0)*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^4)*(1 + tanh(t - GS.t0)))/(Pi*GS.sigmax^3*GS.sigmay)
Sz_xxyy(t, x, y, GS ::GaussianSource) = (GS.Amp*((exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmay^2 - (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^4)/GS.sigmax^2 + ((x - GS.x0)^2*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^4))/GS.sigmax^4)*(1 + tanh(t - GS.t0)))/(2*Pi*GS.sigmax*GS.sigmay)

Sz_t(t, x, y, GS ::GaussianSource) = (GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*sech(t - GS.t0)^2)/(2*Pi*GS.sigmax*GS.sigmay)
Sz_tx(t, x, y, GS ::GaussianSource) = -1/2*(GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)*sech(t - GS.t0)^2)/(Pi*GS.sigmax^3*GS.sigmay)
Sz_txx(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmax^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)^2)/GS.sigmax^4)*sech(t - GS.t0)^2)/(2*Pi*GS.sigmax*GS.sigmay)
Sz_ty(t, x, y, GS ::GaussianSource) = -1/2*(GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)*sech(t - GS.t0)^2)/(Pi*GS.sigmax*GS.sigmay^3)
Sz_tyy(t, x, y, GS ::GaussianSource) = (GS.Amp*(-(exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))/GS.sigmay^2) + (exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(y - GS.y0)^2)/GS.sigmay^4)*sech(t - GS.t0)^2)/(2*Pi*GS.sigmax*GS.sigmay)
Sz_txy(t, x, y, GS ::GaussianSource) = (GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*(x - GS.x0)*(y - GS.y0)*sech(t - GS.t0)^2)/(2*Pi*GS.sigmax^3*GS.sigmay^3)
Sz_tt(t, x, y, GS ::GaussianSource) = -((GS.Amp*exp(-1/2*(x - GS.x0)^2/GS.sigmax^2 - (y - GS.y0)^2/(2*GS.sigmay^2))*sech(t - GS.t0)^2*tanh(t - GS.t0))/(Pi*GS.sigmax*GS.sigmay))




@inline getSourcetime(ff::Source) = ff.time


