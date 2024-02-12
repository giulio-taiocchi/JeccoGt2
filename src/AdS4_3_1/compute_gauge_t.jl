
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained,
                       bulkevol::BulkEvolved, deriv::BulkDeriv, gauge::Gauge,
                       cache::HorizonCache, sys::System{Outer}, ::ConstantGauge)
    xi_t = getxi(gauge_t)
    fill!(xi_t, 0)
    nothing
end


function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained,
                       bulkevol::BulkEvolved, deriv::BulkDeriv, gauge::Gauge,
                       cache::HorizonCache, sys::System{Outer}, gaugecondition::Advect_xi)
    _, Nx, Ny = size(sys)

    Dx  = sys.Dx
    Dy  = sys.Dy

    xiGF = getxi(gauge)
    xi_t = getxi(gauge_t)

    vx   = gaugecondition.xi_vx
    vy   = gaugecondition.xi_vy

    @fastmath @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_x  = Dx(xiGF, 1,i,j)
            xi_y  = Dy(xiGF, 1,i,j)

            xi_t[1,i,j] = -vx * xi_x - vy * xi_y
        end
    end

    nothing
end


#= Solving for xi_t

this is a 2D PDE of the type

  (axx Dxx + ayy Dyy + Axy Dxy + bx Dx + by Dy + c Id) xi_t = -S

we first build each operator (Dxx, Dyy, etc) through a Kronecker product and
then overwrite each of them with the corresponding coefficient. we then sum them
all up and solve the linear system.

we use SparseMatrices since these are finite differencing operators. also, there
is no reason to use the same operators as the ones we use elsewhere. here, as
we're solving a very large linear system, it may be better to use second order
accurate operators since the resulting matrix is smaller and therefore much
faster to invert.

since this is a gauge condition, i think it shouldn't affect the overall
convergence order.

=#
function compute_xi_t!(gauge_t::Gauge, bulkconstrain::BulkConstrained,
                       bulkevol::BulkEvolved, deriv::BulkDeriv, gauge::Gauge,
                       cache::HorizonCache, sys::System{Outer},
                       gaugecondition::ConstantAH)
    _, Nx, Ny = size(sys)
    bulk = Bulk(bulkevol, bulkconstrain)

    Du  = sys.Du
    Dx  = sys.Dx
    Dxx = sys.Dxx
    Dy  = sys.Dy
    Dyy = sys.Dyy

    interp = sys.uinterp

    B_uAH      = cache.bulkhorizon.B_uAH
    G_uAH       = cache.bulkhorizon.G_uAH
    S_uAH       = cache.bulkhorizon.S_uAH
    Fx_uAH      = cache.bulkhorizon.Fx_uAH
    Fy_uAH      = cache.bulkhorizon.Fy_uAH
    Sd_uAH      = cache.bulkhorizon.Sd_uAH
    Bd_uAH     = cache.bulkhorizon.Bd_uAH
    Gd_uAH      = cache.bulkhorizon.Gd_uAH
    A_uAH       = cache.bulkhorizon.A_uAH
    Du_B_uAH   = cache.bulkhorizon.Du_B_uAH
    Du_G_uAH    = cache.bulkhorizon.Du_G_uAH
    Du_S_uAH    = cache.bulkhorizon.Du_S_uAH
    Du_Fx_uAH   = cache.bulkhorizon.Du_Fx_uAH
    Du_Fy_uAH   = cache.bulkhorizon.Du_Fy_uAH
    Du_Sd_uAH   = cache.bulkhorizon.Du_Sd_uAH
    Du_Bd_uAH  = cache.bulkhorizon.Du_Bd_uAH
    Du_Gd_uAH   = cache.bulkhorizon.Du_Gd_uAH
    Du_A_uAH    = cache.bulkhorizon.Du_A_uAH
    Duu_B_uAH  = cache.bulkhorizon.Duu_B_uAH
    Duu_G_uAH   = cache.bulkhorizon.Duu_G_uAH
    Duu_S_uAH   = cache.bulkhorizon.Duu_S_uAH
    Duu_Fx_uAH  = cache.bulkhorizon.Duu_Fx_uAH
    Duu_Fy_uAH  = cache.bulkhorizon.Duu_Fy_uAH
    Duu_A_uAH   = cache.bulkhorizon.Duu_A_uAH

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

    xi_t = getxi(gauge_t)

    uAH   = gaugecondition.u_AH
    kappa = gaugecondition.kappa

    u2 = uAH * uAH
    u3 = uAH * uAH * uAH
    u4 = uAH * uAH * uAH * uAH

    # take u-derivatives of Sd, Bd,and Gd. these are not needed in the
    # nested system, so they haven't been computed before. since we need them
    # here, compute them now
    @sync begin
        @spawn mul!(deriv.Du_Sd,  Du,  bulkconstrain.Sd)
        @spawn mul!(deriv.Du_Bd, Du,  bulkconstrain.Bd)
        @spawn mul!(deriv.Du_Gd,  Du,  bulkconstrain.Gd)
    end

    # interpolate bulk functions (and u-derivatives) to the u=uAH surface
    @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            B_uAH[1,i,j]       = interp(view(bulk.B,  :,i,j))(uAH)
            G_uAH[1,i,j]        = interp(view(bulk.G,   :,i,j))(uAH)
            S_uAH[1,i,j]        = interp(view(bulk.S,   :,i,j))(uAH)
            Fx_uAH[1,i,j]       = interp(view(bulk.Fx,  :,i,j))(uAH)
            Fy_uAH[1,i,j]       = interp(view(bulk.Fy,  :,i,j))(uAH)
            Sd_uAH[1,i,j]       = interp(view(bulk.Sd,  :,i,j))(uAH)
            Bd_uAH[1,i,j]      = interp(view(bulk.Bd, :,i,j))(uAH)
            Gd_uAH[1,i,j]       = interp(view(bulk.Gd,  :,i,j))(uAH)
            A_uAH[1,i,j]        = interp(view(bulk.A,   :,i,j))(uAH)

            Du_B_uAH[1,i,j]    = interp(view(deriv.Du_B,  :,i,j))(uAH)
            Du_G_uAH[1,i,j]     = interp(view(deriv.Du_G,   :,i,j))(uAH)
            Du_S_uAH[1,i,j]     = interp(view(deriv.Du_S,   :,i,j))(uAH)
            Du_Fx_uAH[1,i,j]    = interp(view(deriv.Du_Fx,  :,i,j))(uAH)
            Du_Fy_uAH[1,i,j]    = interp(view(deriv.Du_Fy,  :,i,j))(uAH)
            Du_Sd_uAH[1,i,j]    = interp(view(deriv.Du_Sd,  :,i,j))(uAH)
            Du_Bd_uAH[1,i,j]   = interp(view(deriv.Du_Bd, :,i,j))(uAH)
            Du_Gd_uAH[1,i,j]    = interp(view(deriv.Du_Gd,  :,i,j))(uAH)
            Du_A_uAH[1,i,j]     = interp(view(deriv.Du_A,   :,i,j))(uAH)

            Duu_B_uAH[1,i,j]   = interp(view(deriv.Duu_B,  :,i,j))(uAH)
            Duu_G_uAH[1,i,j]    = interp(view(deriv.Duu_G,   :,i,j))(uAH)
            Duu_S_uAH[1,i,j]    = interp(view(deriv.Duu_S,   :,i,j))(uAH)
            Duu_Fx_uAH[1,i,j]   = interp(view(deriv.Duu_Fx,  :,i,j))(uAH)
            Duu_Fy_uAH[1,i,j]   = interp(view(deriv.Duu_Fy,  :,i,j))(uAH)
            Duu_A_uAH[1,i,j]    = interp(view(deriv.Duu_A,   :,i,j))(uAH)
        end
    end

    ind2D  = LinearIndices(B_uAH[1,:,:])

    # coefficients of the derivative operators
    @fastmath @inbounds Threads.@threads for j in 1:Ny
        @inbounds for i in 1:Nx
            idx   = ind2D[i,j]

            xi    = gauge.xi[1,i,j]
            xi_x  = Dx(gauge.xi, 1,i,j)
            xi_y  = Dy(gauge.xi, 1,i,j)
            xi_xx = Dxx(gauge.xi, 1,i,j)
            xi_yy = Dyy(gauge.xi, 1,i,j)
            xi_xy = Dx(Dy, gauge.xi, 1,i,j)

            B    = B_uAH[1,i,j]
            G     = G_uAH[1,i,j]
            S     = S_uAH[1,i,j]
            Fx    = Fx_uAH[1,i,j]
            Fy    = Fy_uAH[1,i,j]
            Sd    = Sd_uAH[1,i,j]
            Bd   = Bd_uAH[1,i,j]
            Gd    = Gd_uAH[1,i,j]
            A     = A_uAH[1,i,j]

            # r derivatives

            Bp        = -u2 * Du_B_uAH[1,i,j]
            Gp         = -u2 * Du_G_uAH[1,i,j]
            Sp         = -u2 * Du_S_uAH[1,i,j]
            Fxp        = -u2 * Du_Fx_uAH[1,i,j]
            Fyp        = -u2 * Du_Fy_uAH[1,i,j]
            Sdp        = -u2 * Du_Sd_uAH[1,i,j]
            Bdp       = -u2 * Du_Bd_uAH[1,i,j]
            Gdp        = -u2 * Du_Gd_uAH[1,i,j]
            Ap         = -u2 * Du_A_uAH[1,i,j]

            Bpp       = 2*u3 * Du_B_uAH[1,i,j]  + u4 * Duu_B_uAH[1,i,j]
            Gpp        = 2*u3 * Du_G_uAH[1,i,j]   + u4 * Duu_G_uAH[1,i,j]
            Spp        = 2*u3 * Du_S_uAH[1,i,j]   + u4 * Duu_S_uAH[1,i,j]
            Fxpp       = 2*u3 * Du_Fx_uAH[1,i,j]  + u4 * Duu_Fx_uAH[1,i,j]
            Fypp       = 2*u3 * Du_Fy_uAH[1,i,j]  + u4 * Duu_Fy_uAH[1,i,j]
            App        = 2*u3 * Du_A_uAH[1,i,j]   + u4 * Duu_A_uAH[1,i,j]

            # x and y derivatives

            B_x    = Dx(B_uAH,   1,i,j)
            G_x     = Dx(G_uAH,    1,i,j)
            S_x     = Dx(S_uAH,    1,i,j)
            Fx_x    = Dx(Fx_uAH,   1,i,j)
            Fy_x    = Dx(Fy_uAH,   1,i,j)
            Sd_x    = Dx(Sd_uAH,   1,i,j)
            Bd_x   = Dx(Bd_uAH,  1,i,j)
            Gd_x    = Dx(Gd_uAH,   1,i,j)
            A_x     = Dx(A_uAH,    1,i,j)

            B_y    = Dy(B_uAH,   1,i,j)
            G_y     = Dy(G_uAH,    1,i,j)
            S_y     = Dy(S_uAH,    1,i,j)
            Fx_y    = Dy(Fx_uAH,   1,i,j)
            Fy_y    = Dy(Fy_uAH,   1,i,j)
            Sd_y    = Dy(Sd_uAH,   1,i,j)
            Bd_y   = Dy(Bd_uAH,  1,i,j)
            Gd_y    = Dy(Gd_uAH,   1,i,j)
            A_y     = Dy(A_uAH,    1,i,j)

            Bp_x   = -u2 * Dx(Du_B_uAH, 1,i,j)
            Gp_x    = -u2 * Dx(Du_G_uAH,  1,i,j)
            Sp_x    = -u2 * Dx(Du_S_uAH,  1,i,j)
            Fxp_x   = -u2 * Dx(Du_Fx_uAH, 1,i,j)
            Fyp_x   = -u2 * Dx(Du_Fy_uAH, 1,i,j)
            Ap_x    = -u2 * Dx(Du_A_uAH,  1,i,j)

            Bp_y   = -u2 * Dy(Du_B_uAH, 1,i,j)
            Gp_y    = -u2 * Dy(Du_G_uAH,  1,i,j)
            Sp_y    = -u2 * Dy(Du_S_uAH,  1,i,j)
            Fxp_y   = -u2 * Dy(Du_Fx_uAH, 1,i,j)
            Fyp_y   = -u2 * Dy(Du_Fy_uAH, 1,i,j)
            Ap_y    = -u2 * Dy(Du_A_uAH,  1,i,j)

            Fy_xx   = Dxx(Fy_uAH,  1,i,j)
            Fx_xx   = Dxx(Fx_uAH,  1,i,j)
            A_xx    = Dxx(A_uAH,   1,i,j)
            B_xx    = Dxx(B_uAH,   1,i,j)
            G_xx    = Dxx(G_uAH,   1,i,j)
            S_xx    = Dxx(S_uAH,   1,i,j)

            Fx_yy   = Dyy(Fx_uAH,  1,i,j)
            Fy_yy   = Dyy(Fy_uAH,  1,i,j)
            A_yy    = Dyy(A_uAH,   1,i,j)
            B_yy    = Dyy(B_uAH,   1,i,j)
            G_yy    = Dyy(G_uAH,   1,i,j)
            S_yy    = Dyy(S_uAH,   1,i,j)

            Fx_xy   = Dx(Dy, Fx_uAH, 1,i,j)
            Fy_xy   = Dx(Dy, Fy_uAH, 1,i,j)
            A_xy    = Dx(Dy, A_uAH,  1,i,j)
            B_xy    = Dx(Dy, B_uAH,  1,i,j)
            G_xy    = Dx(Dy, G_uAH,  1,i,j)
            S_xy    = Dx(Dy, S_uAH,  1,i,j)

            vars =  (
                kappa, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
                B   , G   ,  S    , Fx    , Fy    , Sd ,  Bd  , Gd, A   ,
                Bp  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,  Bdp , Gdp,       Ap  ,
                Bpp , Gpp ,        Spp  , Fxpp  , Fypp  ,                                App ,
                B_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x, Bd_x,  Gd_x,      A_x ,
	        B_y ,G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y, Bd_y,  Gd_y,      A_y ,
                Bp_x, Gp_x,        Sp_x , Fxp_x , Fyp_x ,                                Ap_x,
                Bp_y, Gp_y,        Sp_y , Fxp_y , Fyp_y ,                                Ap_y,
                                                          Fy_xx ,                                A_xx,
                                                  Fx_yy ,                                        A_yy,
                                                  Fx_xy , Fy_xy ,                                A_xy,Fx_xx,Fy_yy,B_xx,B_yy,B_xy,G_xx,G_yy,G_xy,S_xx,S_yy,S_xy
            )

            a11, a22, a12, b1, b2, c, SS = xi_t_eq_coeff(vars, sys.gridtype)

            axx[idx]   = a11
            ayy[idx]   = a22
            axy[idx]   = a12
            bx[idx]    = b1
            by[idx]    = b2
            cc[idx]    = c

            b_vec[idx] = -SS
        end
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
    A_mat = _Dxx_2D + _Dyy_2D + _Dxy_2D + _Dx_2D + _Dy_2D + ccId

    # solve system
    A_fact = lu(A_mat)
    ldiv!(A_fact, b_vec)

    copyto!(xi_t, b_vec)

    nothing
end
