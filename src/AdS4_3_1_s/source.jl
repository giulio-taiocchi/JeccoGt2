
Base.@kwdef struct NoSource{T} <: Source 
	time :: T = 0.0
end

Sz(x, y, ::NoSource) = 1
Sz_t(x, y, ::NoSource) = 0





