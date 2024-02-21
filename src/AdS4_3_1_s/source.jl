
Base.@kwdef struct NoSource{T} <: Source{T} 
	time :: T = 0.0
end

Sz(x, y, ::NoSource) = 1
Sz_t(x, y, ::NoSource) = 0


@inline getSourcetime(ff::Source) = ff.time


