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

Sz(t, x, y, GS ::GaussianSource) = 1 + (GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_x(t, x, y, GS ::GaussianSource) = -1/2*(GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)*(1 + tanh(t - GaussianSource.t0)))/(Pi*GaussianSource.sigmax^3*GaussianSource.sigmay)
Sz_y(t, x, y, GS ::GaussianSource) =-1/2*(GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)*(1 + tanh(t - GaussianSource.t0)))/(Pi*GaussianSource.sigmax*GaussianSource.sigmay^3)

Sz_xx(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmax^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^2)/GaussianSource.sigmax^4)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_yy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmay^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^4)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_xxx(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*((3*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0))/GaussianSource.sigmax^4 - (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^3)/GaussianSource.sigmax^6)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_yyy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*((3*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0))/GaussianSource.sigmay^4 - (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^3)/GaussianSource.sigmay^6)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_xxxx(t, x, y, GS ::GaussianSource) =(GaussianSource.Amp*((3*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2)))/GaussianSource.sigmax^4 - (6*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^2)/GaussianSource.sigmax^6 + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^4)/GaussianSource.sigmax^8)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_yyyy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*((3*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2)))/GaussianSource.sigmay^4 - (6*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^6 + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^4)/GaussianSource.sigmay^8)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)

Sz_xy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)*(y - GaussianSource.y0)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax^3*GaussianSource.sigmay^3)
Sz_xxy(t, x, y, GS ::GaussianSource) =(GaussianSource.Amp*((exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0))/(GaussianSource.sigmax^2*GaussianSource.sigmay^2) - (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^2*(y - GaussianSource.y0))/(GaussianSource.sigmax^4*GaussianSource.sigmay^2))*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_xyy(t, x, y, GS ::GaussianSource) = -1/2*(GaussianSource.Amp*(x - GaussianSource.x0)*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmay^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^4)*(1 + tanh(t - GaussianSource.t0)))/(Pi*GaussianSource.sigmax^3*GaussianSource.sigmay)
Sz_xxyy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*((exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmay^2 - (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^4)/GaussianSource.sigmax^2 + ((x - GaussianSource.x0)^2*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmay^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^4))/GaussianSource.sigmax^4)*(1 + tanh(t - GaussianSource.t0)))/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)

Sz_t(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*sech(t - GaussianSource.t0)^2)/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_tx(t, x, y, GS ::GaussianSource) = -1/2*(GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)*sech(t - GaussianSource.t0)^2)/(Pi*GaussianSource.sigmax^3*GaussianSource.sigmay)
Sz_txx(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmax^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)^2)/GaussianSource.sigmax^4)*sech(t - GaussianSource.t0)^2)/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_ty(t, x, y, GS ::GaussianSource) = -1/2*(GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)*sech(t - GaussianSource.t0)^2)/(Pi*GaussianSource.sigmax*GaussianSource.sigmay^3)
Sz_tyy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*(-(exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))/GaussianSource.sigmay^2) + (exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(y - GaussianSource.y0)^2)/GaussianSource.sigmay^4)*sech(t - GaussianSource.t0)^2)/(2*Pi*GaussianSource.sigmax*GaussianSource.sigmay)
Sz_txy(t, x, y, GS ::GaussianSource) = (GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*(x - GaussianSource.x0)*(y - GaussianSource.y0)*sech(t - GaussianSource.t0)^2)/(2*Pi*GaussianSource.sigmax^3*GaussianSource.sigmay^3)
Sz_tt(t, x, y, GS ::GaussianSource) = -((GaussianSource.Amp*exp(-1/2*(x - GaussianSource.x0)^2/GaussianSource.sigmax^2 - (y - GaussianSource.y0)^2/(2*GaussianSource.sigmay^2))*sech(t - GaussianSource.t0)^2*tanh(t - GaussianSource.t0))/(Pi*GaussianSource.sigmax*GaussianSource.sigmay))




@inline getSourcetime(ff::Source) = ff.time


