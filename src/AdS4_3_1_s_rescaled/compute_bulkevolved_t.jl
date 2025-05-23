
function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System, ::EvolTest0)

    B_t, G_t = unpack(bulkevol_t)
    # B  , G   = unpack(bulkevol)

    fill!(B_t,  0)
    fill!(G_t,   0)

    nothing
end


function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Inner}, evoleq::AffineNull)
    uu  = sys.ucoord
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy
    L = evoleq.L

    Nu, Nx, Ny = size(sys)

    B_t, G_t= unpack(bulkevol_t)
    source = evoleq.source
    test = source.time

    # u = 0
    @fastmath @inbounds @threads for j in 1:Ny
        @inbounds @simd for i in 1:Nx
            xi     = gauge.xi[1,i,j]

            B     = bulkevol.B[1,i,j]
            G      = bulkevol.G[1,i,j]

            B_u   = Du(bulkevol.B, 1,i,j)/L
            G_u    = Du(bulkevol.G,  1,i,j)/L

            Bd_u  = Du(bulkconstrain.Bd, 1,i,j)/L
            Gd_u   = Du(bulkconstrain.Gd,  1,i,j)/L
	    
	    # the u->0 limit should stay the same since the source appears at higher orders

            B_t[1,i,j]  = L*(Bd_u + 2 * B_u  + 3 * B * xi)
            G_t[1,i,j]  = L*(Gd_u + 2 * G_u  + 3 * G * xi)
        end
    end
  
    # remaining inner grid points
    @fastmath @inbounds for j in 1:Ny
    	y = sys.ycoord[j]
        @inbounds for i in 1:Nx
            x = sys.xcoord[i]
            xi    = gauge.xi[1,i,j]
            xi_t  = gauge_t.xi[1,i,j]/L
            @inbounds @simd for a in 2:Nu
            
                S0 = Sz(test, x, y, source)
                S0_x = Sz_x(test, x, y, source)/L
                S0_y = Sz_y(test, x, y, source)/L
                S0_t = Sz_t(test, x, y, source)/L
                S0_xx = Sz_xx(test, x, y, source)/L/L
                S0_yy = Sz_yy(test, x, y, source)/L/L
                S0_tt = Sz_tt(test, x, y, source)/L/L
                u      = uu[a]
                u3     = u * u * u 

                B     = bulkevol.B[a,i,j]
                G      = bulkevol.G[a,i,j]

                Bd    = bulkconstrain.Bd[a,i,j]
                Gd     = bulkconstrain.Gd[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                B_u   = Du(bulkevol.B, a,i,j)/L
                G_u    = Du(bulkevol.G,  a,i,j)/L

		B_t[a,i,j] = L*(2*Bd + ((3*B + u*(B_u))*(S0^2*u^2*(S0_t)^2 - 2*S0^3*u^2*(S0_tt) + u^2*((S0_x)^2 + (S0_y)^2) - S0*u^2*((S0_xx) + (S0_yy)) + S0^4*(A*u^3 + (1 + u*xi)^2 - 2*u^2*(xi_t))))/S0^4)/(2*u)#((3* B + u * B_u) *
                               #(-2 * u * u * xi_t + A * u3 +
                                #(xi * u + 1) * (xi * u + 1))+ 2 * Bd)/(2 * u) #
                                

		G_t[a,i,j] = L*(2*Gd + ((3*G + u*(G_u))*(S0^2*u^2*(S0_t)^2 - 2*S0^3*u^2*(S0_tt) + u^2*((S0_x)^2 + (S0_y)^2) - S0*u^2*((S0_xx) + (S0_yy)) + S0^4*(A*u^3 + (1 + u*xi)^2 - 2*u^2*(xi_t))))/S0^4)/(2*u)    #((3 * G + u * G_u) *
                               #(-2 * u * u * xi_t + A * u3 +
                                #(xi * u + 1) * (xi * u + 1)) + 2 * Gd)/(2 * u)     #                    
            end
        end
    end


    nothing
end


function compute_bulkevolved_t!(bulkevol_t::BulkEvolved,
                                bulkconstrain::BulkConstrained, gauge_t::Gauge,
                                bulkevol::BulkEvolved, boundary::Boundary,
                                gauge::Gauge, sys::System{Outer}, evoleq::AffineNull)
    uu  = sys.ucoord
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy
    L = evoleq.L
    Nu, Nx, Ny = size(sys)

    B_t, G_t = unpack(bulkevol_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t  = gauge_t.xi[1,i,j]/L
            @inbounds @simd for a in 1:Nu
                u      = uu[a]
                u2     = u * u

                Bd    = bulkconstrain.Bd[a,i,j]
                Gd     = bulkconstrain.Gd[a,i,j]
                A      = bulkconstrain.A[a,i,j]

                B_u   = Du(bulkevol.B, a,i,j)/L
                G_u    = Du(bulkevol.G,  a,i,j)/L

		B_t[a,i,j] = L*(Bd + u2 * (A/2 - xi_t) * B_u)
		G_t[a,i,j] = L*(Gd + u2 * (A/2 - xi_t) * G_u)
            end
        end
    end

    nothing
end

function sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t::Gauge,
                           systems::SystemPartition, evoleq::AffineNull)
    Nsys = length(systems)

    @inbounds for i in 1:Nsys-1
        sync_bulkevolved!(bulkevols_t[i], bulkevols_t[i+1], bulkconstrains[i+1],
                          gauge_t, systems[i], systems[i+1], evoleq)
    end

    nothing
end

function sync_bulkevolved!(bulkevol1_t::BulkEvolved, bulkevol2_t::BulkEvolved,
                           bulkconstrain2::BulkConstrained, gauge_t::Gauge,
                           sys1::System{Outer}, sys2::System{Outer}, evoleq::AffineNull)
    u0 = sys2.ucoord[1]

    _, Nx, Ny = size(sys2)
    L = evoleq.L
    B1_t, G1_t = unpack(bulkevol1_t)
    B2_t, G2_t = unpack(bulkevol2_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]/L
            A      = bulkconstrain2.A[1,i,j]
  
  
  	    #CHECK HERE!
  	    
  	    
            # characteristic speed
            c = u0 * u0 * (A/2 - xi_t)

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                B1_t[end,i,j]  = B2_t[1,i,j]
                G1_t[end,i,j]   = G2_t[1,i,j]
            elseif c < 0
                B2_t[1,i,j]  = B1_t[end,i,j]
                G2_t[1,i,j]   = G1_t[end,i,j]
            end
        end
    end

    nothing
end

function sync_bulkevolved!(bulkevol1_t::BulkEvolved, bulkevol2_t::BulkEvolved,
                           bulkconstrain2::BulkConstrained, gauge_t::Gauge,
                           sys1::System{Inner}, sys2::System{Outer}, evoleq::AffineNull)
    u0  = sys2.ucoord[1]
    u02 = u0 * u0
    u03 = u0 * u02
    u04 = u02 * u02
    L = evoleq.L
    _, Nx, Ny = size(sys2)

    B1_t, G1_t = unpack(bulkevol1_t)
    B2_t, G2_t = unpack(bulkevol2_t)

    @fastmath @inbounds for j in 1:Ny
        @inbounds for i in 1:Nx
            xi_t   = gauge_t.xi[1,i,j]/L
            A      = bulkconstrain2.A[1,i,j]

            # characteristic speed
            c = u0 * u0 * (A/2 - xi_t)

            # if c > 0, mode is entering grid1 from grid2;
            # if c < 0, mode is entering grid2 from grid1.
            # we assume here that grids merely touch at the interface
            if c > 0
                B1_t[end,i,j]  = B2_t[1,i,j] / u03
                G1_t[end,i,j]   = G2_t[1,i,j] / u03
            elseif c < 0
                B2_t[1,i,j]  = u03 * B1_t[end,i,j]
                G2_t[1,i,j]   = u03 * G1_t[end,i,j]
            end
        end
    end


    nothing
end
