
function compute_residual_AH!(res::Array, sigma::Array,
                              gauge::Gauge, cache::HorizonCache,
                              sys::System{Outer}, evoleq::EvolutionEquations)
    _, Nx, Ny = size(sys)

    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy
    L = evoleq.L

    # we assume here that the bulk functions and their u-derivatives have
    # already been interpolated to the uAH surface (u = 1/sigma)

    B_uAH      = cache.bulkhorizon.B_uAH
    G_uAH       = cache.bulkhorizon.G_uAH
    S_uAH       = cache.bulkhorizon.S_uAH
    Fx_uAH      = cache.bulkhorizon.Fx_uAH
    Fy_uAH      = cache.bulkhorizon.Fy_uAH
    Sd_uAH      = cache.bulkhorizon.Sd_uAH

    Du_B_uAH   = cache.bulkhorizon.Du_B_uAH
    Du_G_uAH    = cache.bulkhorizon.Du_G_uAH
    Du_S_uAH    = cache.bulkhorizon.Du_S_uAH
    Du_Fx_uAH   = cache.bulkhorizon.Du_Fx_uAH
    Du_Fy_uAH   = cache.bulkhorizon.Du_Fy_uAH
    Du_Sd_uAH   = cache.bulkhorizon.Du_Sd_uAH

    Duu_B_uAH  = cache.bulkhorizon.Duu_B_uAH
    Duu_G_uAH   = cache.bulkhorizon.Duu_G_uAH
    Duu_S_uAH   = cache.bulkhorizon.Duu_S_uAH
    Duu_Fx_uAH  = cache.bulkhorizon.Duu_Fx_uAH
    Duu_Fy_uAH  = cache.bulkhorizon.Duu_Fy_uAH

    # compute AH equation residual
    @fastmath @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            uAH = 1 / sigma[1,i,j]
            u2  = uAH * uAH
            u3  = uAH * uAH * uAH
            u4  = uAH * uAH * uAH * uAH

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)/L
            xi_y  = Dy(gauge.xi, 1,i,j)/L
            xi_xx = Dxx(gauge.xi, 1,i,j)/L/L
            xi_yy = Dyy(gauge.xi, 1,i,j)/L/L
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)/L/L

            sigma0    = sigma[1,i,j]
            sigma0_x  = Dx(sigma, 1,i,j)/L
            sigma0_y  = Dy(sigma, 1,i,j)/L
            sigma0_xx = Dxx(sigma, 1,i,j)/L/L
            sigma0_yy = Dyy(sigma, 1,i,j)/L/L
            sigma0_xy = Dx(Dy, sigma, 1,i,j)/L/L

            B    = B_uAH[1,i,j]
            G     = G_uAH[1,i,j]
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]

            # r derivatives

            Bp        = -u2 * Du_B_uAH[1,i,j]/L
            Gp         = -u2 * Du_G_uAH[1,i,j]/L
            Sp         = -u2 * Du_S_uAH[1,i,j]/L
            Fxp        = -u2 * Du_Fx_uAH[1,i,j]/L
            Fyp        = -u2 * Du_Fy_uAH[1,i,j]/L
            Sdp        = -u2 * Du_Sd_uAH[1,i,j]/L

            Bpp       = 2*u3 * Du_B_uAH[1,i,j]/L  + u4 * Duu_B_uAH[1,i,j]/L/L
            Gpp        = 2*u3 * Du_G_uAH[1,i,j] /L  + u4 * Duu_G_uAH[1,i,j]/L/L
            Spp        = 2*u3 * Du_S_uAH[1,i,j] /L  + u4 * Duu_S_uAH[1,i,j]/L/L
            Fxpp       = 2*u3 * Du_Fx_uAH[1,i,j]/L  + u4 * Duu_Fx_uAH[1,i,j]/L/L
            Fypp       = 2*u3 * Du_Fy_uAH[1,i,j]/L  + u4 * Duu_Fy_uAH[1,i,j]/L/L

            # x and y derivatives

            B_x    = (Dx(B_uAH,   1,i,j) - Bp * sigma0_x)/L
            G_x     = (Dx(G_uAH,    1,i,j) -  Gp * sigma0_x)/L
            S_x     = (Dx(S_uAH,    1,i,j) -  Sp * sigma0_x)/L
            Fx_x    = (Dx(Fx_uAH,   1,i,j) - Fxp * sigma0_x)/L
            Fy_x    = (Dx(Fy_uAH,   1,i,j) - Fyp * sigma0_x)/L
            Sd_x    = (Dx(Sd_uAH,   1,i,j) - Sdp * sigma0_x)/L

            B_y    = (Dy(B_uAH,   1,i,j) - Bp * sigma0_y)/L
            G_y     = (Dy(G_uAH,    1,i,j) -  Gp * sigma0_y)/L
            S_y     = (Dy(S_uAH,    1,i,j) -  Sp * sigma0_y)/L
            Fx_y    = (Dy(Fx_uAH,   1,i,j) - Fxp * sigma0_y)/L
            Fy_y    = (Dy(Fy_uAH,   1,i,j) - Fyp * sigma0_y)/L
            Sd_y    = (Dy(Sd_uAH,   1,i,j) - Sdp * sigma0_y)/L

            Bp_x   = (-u2 * Dx(Du_B_uAH, 1,i,j) - Bpp * sigma0_x + 2 * u3 * sigma0_x * Du_B_uAH[1,i,j])/L
            Gp_x    = (-u2 * Dx(Du_G_uAH,  1,i,j) -  Gpp * sigma0_x + 2 * u3 * sigma0_x *  Du_G_uAH[1,i,j])/L
            Sp_x    = (-u2 * Dx(Du_S_uAH,  1,i,j) -  Spp * sigma0_x + 2 * u3 * sigma0_x *  Du_S_uAH[1,i,j])/L
            Fxp_x   = (-u2 * Dx(Du_Fx_uAH, 1,i,j) - Fxpp * sigma0_x + 2 * u3 * sigma0_x * Du_Fx_uAH[1,i,j])/L
            Fyp_x   = (-u2 * Dx(Du_Fy_uAH, 1,i,j) - Fypp * sigma0_x + 2 * u3 * sigma0_x * Du_Fy_uAH[1,i,j])/L

            Bp_y   = (-u2 * Dy(Du_B_uAH, 1,i,j) - Bpp * sigma0_y + 2 * u3 * sigma0_y * Du_B_uAH[1,i,j])/L
            Gp_y    = (-u2 * Dy(Du_G_uAH,  1,i,j) -  Gpp * sigma0_y + 2 * u3 * sigma0_y *  Du_G_uAH[1,i,j])/L
            Sp_y    = (-u2 * Dy(Du_S_uAH,  1,i,j) -  Spp * sigma0_y + 2 * u3 * sigma0_y *  Du_S_uAH[1,i,j])/L
            Fxp_y   = (-u2 * Dy(Du_Fx_uAH, 1,i,j) - Fxpp * sigma0_y + 2 * u3 * sigma0_y * Du_Fx_uAH[1,i,j])/L
            Fyp_y   = (-u2 * Dy(Du_Fy_uAH, 1,i,j) - Fypp * sigma0_y + 2 * u3 * sigma0_y * Du_Fy_uAH[1,i,j])/L


            vars = (
                sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
                xi    , xi_x    , xi_y    , xi_xx    , xi_yy    , xi_xy,
                B   , G   ,  S    , Fx    , Fy    , Sd ,
                Bp  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
                Bpp , Gpp ,  Spp  , Fxpp  , Fypp  ,
                B_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
                B_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
            )

            res[1,i,j] = AH_eq_res(vars, sys.gridtype)
        end
    end

    res
end

function compute_coeffs_AH!(sigma::Array, gauge::Gauge, cache::HorizonCache,
                            sys::System{Outer}, evoleq::EvolutionEquations)
    _, Nx, Ny = size(sys)

    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy
    L = evoleq.L

    B_uAH      = cache.bulkhorizon.B_uAH
    G_uAH       = cache.bulkhorizon.G_uAH
    S_uAH       = cache.bulkhorizon.S_uAH
    Fx_uAH      = cache.bulkhorizon.Fx_uAH
    Fy_uAH      = cache.bulkhorizon.Fy_uAH
    Sd_uAH      = cache.bulkhorizon.Sd_uAH

    Du_B_uAH   = cache.bulkhorizon.Du_B_uAH
    Du_G_uAH    = cache.bulkhorizon.Du_G_uAH
    Du_S_uAH    = cache.bulkhorizon.Du_S_uAH
    Du_Fx_uAH   = cache.bulkhorizon.Du_Fx_uAH
    Du_Fy_uAH   = cache.bulkhorizon.Du_Fy_uAH
    Du_Sd_uAH   = cache.bulkhorizon.Du_Sd_uAH

    Duu_B_uAH  = cache.bulkhorizon.Duu_B_uAH
    Duu_G_uAH   = cache.bulkhorizon.Duu_G_uAH
    Duu_S_uAH   = cache.bulkhorizon.Duu_S_uAH
    Duu_Fx_uAH  = cache.bulkhorizon.Duu_Fx_uAH
    Duu_Fy_uAH  = cache.bulkhorizon.Duu_Fy_uAH

    axx         = cache.axx
    ayy         = cache.ayy
    axy         = cache.axy
    bx          = cache.bx
    by          = cache.by
    cc          = cache.cc
    #b_vec       = cache.b_vec

    ind2D  = LinearIndices(B_uAH[1,:,:])

    # coefficients of the derivative operators
    @fastmath @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            idx   = ind2D[i,j]

            uAH   = 1 / sigma[1,i,j]
            u2    = uAH * uAH
            u3    = uAH * uAH * uAH
            u4    = uAH * uAH * uAH * uAH

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)/L
            xi_y  = Dy(gauge.xi, 1,i,j)/L
            xi_xx = Dxx(gauge.xi, 1,i,j)/L/L
            xi_yy = Dyy(gauge.xi, 1,i,j)/L/L
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)/L/L

            sigma0    = sigma[1,i,j]
            sigma0_x  = Dx(sigma, 1,i,j)/L
            sigma0_y  = Dy(sigma, 1,i,j)/L
            sigma0_xx = Dxx(sigma, 1,i,j)/L/L
            sigma0_yy = Dyy(sigma, 1,i,j)/L/L
            sigma0_xy = Dx(Dy, sigma, 1,i,j)/L/L

            B    = B_uAH[1,i,j]
            G     = G_uAH[1,i,j]
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]

            # r derivatives

            Bp        = -u2 * Du_B_uAH[1,i,j]/L
            Gp         = -u2 * Du_G_uAH[1,i,j]/L
            Sp         = -u2 * Du_S_uAH[1,i,j]/L
            Fxp        = -u2 * Du_Fx_uAH[1,i,j]/L
            Fyp        = -u2 * Du_Fy_uAH[1,i,j]/L
            Sdp        = -u2 * Du_Sd_uAH[1,i,j]/L

            Bpp       = 2*u3 * Du_B_uAH[1,i,j] /L + u4 * Duu_B_uAH[1,i,j]/L/L
            Gpp        = 2*u3 * Du_G_uAH[1,i,j]/L   + u4 * Duu_G_uAH[1,i,j]/L/L
            Spp        = 2*u3 * Du_S_uAH[1,i,j]/L   + u4 * Duu_S_uAH[1,i,j]/L/L
            Fxpp       = 2*u3 * Du_Fx_uAH[1,i,j]/L  + u4 * Duu_Fx_uAH[1,i,j]/L/L
            Fypp       = 2*u3 * Du_Fy_uAH[1,i,j]  /L+ u4 * Duu_Fy_uAH[1,i,j]/L/L

            # x and y derivatives

            B_x    = (Dx(B_uAH,   1,i,j) - Bp * sigma0_x)/L
            G_x     = (Dx(G_uAH,    1,i,j) -  Gp * sigma0_x)/L
            S_x     = (Dx(S_uAH,    1,i,j) -  Sp * sigma0_x)/L
            Fx_x    = (Dx(Fx_uAH,   1,i,j) - Fxp * sigma0_x)/L
            Fy_x    = (Dx(Fy_uAH,   1,i,j) - Fyp * sigma0_x)/L
            Sd_x    = (Dx(Sd_uAH,   1,i,j) - Sdp * sigma0_x)/L

            B_y    = (Dy(B_uAH,   1,i,j) - Bp * sigma0_y)/L
            G_y     = (Dy(G_uAH,    1,i,j) -  Gp * sigma0_y)/L
            S_y     = (Dy(S_uAH,    1,i,j) -  Sp * sigma0_y)/L
            Fx_y    = (Dy(Fx_uAH,   1,i,j) - Fxp * sigma0_y)/L
            Fy_y    = (Dy(Fy_uAH,   1,i,j) - Fyp * sigma0_y)/L
            Sd_y    = (Dy(Sd_uAH,   1,i,j) - Sdp * sigma0_y)/L

            Bp_x   = (-u2 * Dx(Du_B_uAH/L, 1,i,j) - Bpp * sigma0_x + 2 * u3 * sigma0_x * Du_B_uAH[1,i,j]/L)/L
            Gp_x    = (-u2 * Dx(Du_G_uAH/L,  1,i,j) -  Gpp * sigma0_x + 2 * u3 * sigma0_x *  Du_G_uAH[1,i,j]/L)/L
            Sp_x    = (-u2 * Dx(Du_S_uAH/L,  1,i,j) -  Spp * sigma0_x + 2 * u3 * sigma0_x *  Du_S_uAH[1,i,j]/L)/L
            Fxp_x   = (-u2 * Dx(Du_Fx_uAH/L, 1,i,j) - Fxpp * sigma0_x + 2 * u3 * sigma0_x * Du_Fx_uAH[1,i,j]/L)/L
            Fyp_x   = (-u2 * Dx(Du_Fy_uAH/L, 1,i,j) - Fypp * sigma0_x + 2 * u3 * sigma0_x * Du_Fy_uAH[1,i,j]/L)/L

            Bp_y   = (-u2 * Dy(Du_B_uAH/L, 1,i,j) - Bpp * sigma0_y + 2 * u3 * sigma0_y * Du_B_uAH[1,i,j])/L
            Gp_y    = (-u2 * Dy(Du_G_uAH/L,  1,i,j) -  Gpp * sigma0_y + 2 * u3 * sigma0_y *  Du_G_uAH[1,i,j]/L)/L
            Sp_y    = (-u2 * Dy(Du_S_uAH/L,  1,i,j) -  Spp * sigma0_y + 2 * u3 * sigma0_y *  Du_S_uAH[1,i,j]/L)/L
            Fxp_y   = (-u2 * Dy(Du_Fx_uAH/L, 1,i,j) - Fxpp * sigma0_y + 2 * u3 * sigma0_y * Du_Fx_uAH[1,i,j]/L)/L
            Fyp_y   = (-u2 * Dy(Du_Fy_uAH/L, 1,i,j) - Fypp * sigma0_y + 2 * u3 * sigma0_y * Du_Fy_uAH[1,i,j]/L)/L

            vars = (
                sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
                xi    , xi_x    , xi_y    , xi_xx    , xi_yy    , xi_xy,
                B   ,  G   ,  S    , Fx    , Fy    , Sd ,
                Bp  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
                Bpp , Gpp ,  Spp  , Fxpp  , Fypp  ,
                B_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	        B_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
                Bp_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
                Bp_y, Gp_y,  Sp_y , Fxp_y , Fyp_y ,
            )

            a11, a22, a12, b1, b2, c = AH_eq_coeff(vars, sys.gridtype)

            axx[idx]   = a11
            ayy[idx]   = a22
            axy[idx]   = a12
            bx[idx]    = b1
            by[idx]    = b2
            cc[idx]    = c
        end
    end

    nothing
end


#= Finding the Apparent Horizon

this is a 2D (non-linear) PDE of the type

  αxx f_xx + αyy f_yy + αxy f_xy + βxx (f_x)^2 + βyy (f_y)^2 + βxy (f_x) (f_y)
   + γx f_x + γy f_y + S = 0

we solve this equation with a Newton-Kantorovich method where, starting with a
guess, we solve the associated linear problem (using the same strategy as in
compute_xi_t!) to improve our guess until a sufficiently precise solution is
reached.
=#
function find_AH!(sigma::Array, bulkconstrain::BulkConstrained,
                  bulkevol::BulkEvolved, deriv::BulkDeriv, gauge::Gauge,
                  cache::HorizonCache, sys::System{Outer}, ahf::AHF, evoleq::EvolutionEquations)
    _, Nx, Ny = size(sys)
    bulk = Bulk(bulkevol, bulkconstrain)

    Du  = sys.Du
    Duu = sys.Duu
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy
    L = evoleq.L

    interp = sys.uinterp

    B_uAH      = cache.bulkhorizon.B_uAH
    G_uAH       = cache.bulkhorizon.G_uAH
    S_uAH       = cache.bulkhorizon.S_uAH
    Fx_uAH      = cache.bulkhorizon.Fx_uAH
    Fy_uAH      = cache.bulkhorizon.Fy_uAH
    Sd_uAH      = cache.bulkhorizon.Sd_uAH

    Du_B_uAH   = cache.bulkhorizon.Du_B_uAH
    Du_G_uAH    = cache.bulkhorizon.Du_G_uAH
    Du_S_uAH    = cache.bulkhorizon.Du_S_uAH
    Du_Fx_uAH   = cache.bulkhorizon.Du_Fx_uAH
    Du_Fy_uAH   = cache.bulkhorizon.Du_Fy_uAH
    Du_Sd_uAH   = cache.bulkhorizon.Du_Sd_uAH

    Duu_B_uAH  = cache.bulkhorizon.Duu_B_uAH
    Duu_G_uAH   = cache.bulkhorizon.Duu_G_uAH
    Duu_S_uAH   = cache.bulkhorizon.Duu_S_uAH
    Duu_Fx_uAH  = cache.bulkhorizon.Duu_Fx_uAH
    Duu_Fy_uAH  = cache.bulkhorizon.Duu_Fy_uAH

    axx         = cache.axx
    ayy         = cache.ayy
    axy         = cache.axy
    bx          = cache.bx
    by          = cache.by
    cc          = cache.cc
    b_vec       = cache.b_vec

    Dx_2D       = cache.Dx_2D
    Dxx_2D      = cache.Dxx_2D
    Dy_2D       = cache.Dy_2D
    Dyy_2D      = cache.Dyy_2D
    Dxy_2D      = cache.Dxy_2D
    _Dx_2D      = cache._Dx_2D
    _Dxx_2D     = cache._Dxx_2D
    _Dy_2D      = cache._Dy_2D
    _Dyy_2D     = cache._Dyy_2D
    _Dxy_2D     = cache._Dxy_2D

    # take all u-derivatives. we assume that they've not been computed beforehand
    @sync begin
        @spawn mul!(deriv.Du_B,  Du,  bulk.B)
        @spawn mul!(deriv.Du_G,   Du,  bulk.G)
        @spawn mul!(deriv.Du_S,   Du,  bulk.S)
        @spawn mul!(deriv.Du_Fx,  Du,  bulk.Fx)
        @spawn mul!(deriv.Du_Fy,  Du,  bulk.Fy)
        @spawn mul!(deriv.Du_Sd,  Du,  bulk.Sd)

        @spawn mul!(deriv.Duu_B, Duu, bulk.B)
        @spawn mul!(deriv.Duu_G,  Duu, bulk.G)
        @spawn mul!(deriv.Duu_S,  Duu, bulk.S)
        @spawn mul!(deriv.Duu_Fx, Duu, bulk.Fx)
        @spawn mul!(deriv.Duu_Fy, Duu, bulk.Fy)
        # @spawn mul!(deriv.Duu_Sd, Duu, bulk.Sd)
    end

    res = similar(sigma)
    f0  = similar(b_vec)

    println("INFO (AH): Looking for the apparent horizon...")
    println("    it \t max_res")
    # start relaxation method
    it = 0
    while true

        # interpolate bulk functions (and u-derivatives) to the 1/u = r = sigma surface
        @inbounds Threads.@threads for j in 1:Ny
            @inbounds for i in 1:Nx
                uAH = 1 / sigma[1,i,j]
                u2  = uAH * uAH
                u3  = uAH * uAH * uAH
                u4  = uAH * uAH * uAH * uAH
                #println("Current AH: $uAH")
                B_uAH[1,i,j]       = interp(view(bulk.B,  :,i,j))(uAH)
                G_uAH[1,i,j]        = interp(view(bulk.G,   :,i,j))(uAH)
                S_uAH[1,i,j]        = interp(view(bulk.S,   :,i,j))(uAH)
                Fx_uAH[1,i,j]       = interp(view(bulk.Fx,  :,i,j))(uAH)
                Fy_uAH[1,i,j]       = interp(view(bulk.Fy,  :,i,j))(uAH)
                Sd_uAH[1,i,j]       = interp(view(bulk.Sd,  :,i,j))(uAH)

                Du_B_uAH[1,i,j]    = interp(view(deriv.Du_B,  :,i,j))(uAH)
                Du_G_uAH[1,i,j]     = interp(view(deriv.Du_G,   :,i,j))(uAH)
                Du_S_uAH[1,i,j]     = interp(view(deriv.Du_S,   :,i,j))(uAH)
                Du_Fx_uAH[1,i,j]    = interp(view(deriv.Du_Fx,  :,i,j))(uAH)
                Du_Fy_uAH[1,i,j]    = interp(view(deriv.Du_Fy,  :,i,j))(uAH)
                Du_Sd_uAH[1,i,j]    = interp(view(deriv.Du_Sd,  :,i,j))(uAH)

                Duu_B_uAH[1,i,j]   = interp(view(deriv.Duu_B,  :,i,j))(uAH)
                Duu_G_uAH[1,i,j]    = interp(view(deriv.Duu_G,   :,i,j))(uAH)
                Duu_S_uAH[1,i,j]    = interp(view(deriv.Duu_S,   :,i,j))(uAH)
                Duu_Fx_uAH[1,i,j]   = interp(view(deriv.Duu_Fx,  :,i,j))(uAH)
                Duu_Fy_uAH[1,i,j]   = interp(view(deriv.Duu_Fy,  :,i,j))(uAH)
            end
        end

        compute_residual_AH!(res, sigma, gauge, cache, sys, evoleq)
        max_res = maximum(abs.(res))
        println("    $it \t $max_res")

        if max_res < ahf.epsilon
            break
        end

        # compute axx, ayy, axy, bx, by and cc coefficients of the linearized
        # equation (stored in the cache struct)
        compute_coeffs_AH!(sigma, gauge, cache, sys, evoleq)
        copyto!(cache.b_vec, res)

        if it >= ahf.itmax
            println("INFO (AH): maximum iteration reached.")
            println("INFO (AH): giving up")
            break
        end

        axxId = Diagonal(axx)
        ayyId = Diagonal(ayy)
        axyId = Diagonal(axy)
        bxId  = Diagonal(bx)
        byId  = Diagonal(by)
        ccId  = Diagonal(cc)

        # build the differential operators by multiplying with the coefficients
        # computed in the loop above (note that Dxx_2D, Dyy_2D, etc, are never
        # overwritten)
        mul!(_Dxx_2D, axxId, Dxx_2D)
        mul!(_Dyy_2D, ayyId, Dyy_2D)
        mul!(_Dxy_2D, axyId, Dxy_2D)
        mul!(_Dx_2D,  bxId,  Dx_2D)
        mul!(_Dy_2D,  byId,  Dy_2D)

        # build actual operator to be inverted
        A_mat = _Dxx_2D/L/L + _Dyy_2D/L/L + _Dxy_2D /L/L+ _Dx_2D/L + _Dy_2D /L+ ccId

        # solve system
        A_fact = lu(A_mat)
        ldiv!(f0, A_fact, b_vec)

        # update solution
        @inbounds for idx in eachindex(sigma)
            sigma[idx] -= f0[idx]
        end
        it += 1

        # check if sigma inside domain
        min_uAH = 1/maximum(sigma)
        max_uAH = 1/minimum(sigma)
        println("    min : $min_uAH")
        println("    max : $max_uAH")
        println(" ")
        if ( min_uAH < sys.ucoord[1] || max_uAH > sys.ucoord[end] )
            println("INFO (AH): guess outside domain, min_uAH = $min_uAH, max_uAH = $max_uAH")
            println("INFO (AH): giving up")
            break
        end

    end # end relaxation loop

    nothing
end
