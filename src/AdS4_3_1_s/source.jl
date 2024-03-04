
Base.@kwdef mutable struct NoSource{T} <: Source
	time :: T = 0.0
end

Sz(t, x, y, ::NoSource) = 1.0
Sz_x(t, x, y, ::NoSource) = 0.0
Sz_y(t, x, y, ::NoSource) = 0.0
Sz_xx(t, x, y, ::NoSource) = 0.0
Sz_yy(t, x, y, ::NoSource) = 0.0
Sz_t(t, x, y, ::NoSource) = 0.0
Sz_tx(t, x, y, ::NoSource) = 0.0
Sz_ty(t, x, y, ::NoSource) = 0.0


@inline getSourcetime(ff::Source) = ff.time


