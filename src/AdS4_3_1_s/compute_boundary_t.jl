
function compute_boundary_t!(boundary_t::Boundary, bulkevol::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System, ::EvolTest0)

    a3_t, fx1_t, fy1_t = unpack(boundary_t)
    # a3  , fx1  , fy1   = unpack(boundary)

    fill!(a3_t,  0)
    fill!(fx1_t, 0)
    fill!(fy1_t, 0)

    nothing
end

function compute_boundary_t!(boundary_t::Boundary, bulk::BulkEvolved,
                             boundary::Boundary, gauge::Gauge, sys::System{Inner},
                             evoleq::AffineNull)
    Du  = sys.Du
    Dx  = sys.Dx
    Dy  = sys.Dy


    _, Nx, Ny = size(sys)

    a3_t, fx1_t, fy1_t = unpack(boundary_t)
    
    source = evoleq.source
    test = source.time
    step = source.step
    fx1_x_vals = zeros(Nx, Ny)
    a3_x_vals = zeros(Nx, Ny)
    b3_x_vals = zeros(Nx, Ny)
    g3_y_vals = zeros(Nx, Ny)
    fy1_y_vals = zeros(Nx, Ny)
    source_vals = zeros(Nx,Ny)
    no_source_terms = zeros(Nx,Ny)

    @fastmath @inbounds for j in 1:Ny
        y = sys.ycoord[j]
        @inbounds for i in 1:Nx
            x = sys.xcoord[i]
            #importing the source and its derivatives
            S0 = Sz(test, x, y, source)
            source_vals[i,j] = S0 
            S0_x = Sz_x(test, x, y, source)
            S0_y = Sz_y(test, x, y, source)
            S0_t = Sz_t(test, x, y, source)
            S0_xx = Sz_xx(test, x, y, source)
            S0_xxx = Sz_xxx(test, x, y, source)
            S0_xxxx = Sz_xxxx(test, x, y, source)
            S0_yy = Sz_yy(test, x, y, source)
            S0_yyy = Sz_yyy(test, x, y, source)
            S0_yyyy = Sz_yyyy(test, x, y, source)
            S0_xy = Sz_xy(test, x, y, source)
            S0_xxyy = Sz_xxyy(test, x, y, source)
            S0_xyy = Sz_xyy(test, x, y, source)
            S0_xxy = Sz_xxy(test, x, y, source)
            S0_tx = Sz_tx(test, x, y, source)
            S0_txx = Sz_txx(test, x, y, source)
            S0_ty = Sz_ty(test, x, y, source)
            S0_tyy = Sz_tyy(test, x, y, source)
            S0_txy = Sz_txy(test, x, y, source)
            
            xi      = gauge.xi[1,i,j]
            xi3     = xi*xi*xi

            xi_x    = Dx(gauge.xi, 1,i,j)
            xi_y    = Dy(gauge.xi, 1,i,j)

	    b13	    = bulk.B[1,i,j]
            b13_x   = Dx(bulk.B, 1,i,j)
            b3_x_vals[i, j] = b13_x
            b13_y   = Dy(bulk.B, 1,i,j)

	    g3	    = bulk.G[1,i,j]
            g3_x    = Dx(bulk.G, 1,i,j)
            g3_y    = Dy(bulk.G, 1,i,j)
            g3_y_vals[i, j] = g3_y
            
	    a3      = boundary.a3[1,i,j]
            a3_x    = Dx(boundary.a3, 1,i,j)
            a3_x_vals[i, j] = a3_x          
            a3_y    = Dy(boundary.a3, 1,i,j)

	    
	    fx1     = boundary.fx1[1,i,j]
            fx1_x   = Dx(boundary.fx1, 1,i,j)
            fy1     = boundary.fy1[1,i,j]
            fy1_y   = Dy(boundary.fy1, 1,i,j)
            fx1_x_vals[i, j] = fx1_x
            fy1_y_vals[i, j] = fy1_y
            no_source_terms[i,j] = (-3*S0^6*(((fx1_x) + (fy1_y)))/(S0^2))
            

            a3_t[1,i,j]  = (20*(S0_x)^4 + 20*(S0_y)^4 - 16*S0*(S0_x)^2*(2*(S0_xx) + (S0_yy)) + S0^2*(5*(S0_xx)^2 + 4*(S0_xy)^2 + 8*(S0_x)*((S0_xxx) + (S0_xyy)) + 6*(S0_xx)*(S0_yy) + 5*(S0_yy)^2) - 8*(S0_y)^2*(-5*(S0_x)^2 + 2*S0*((S0_xx) + 2*(S0_yy))) + 8*S0*(S0_y)*(-4*(S0_x)*(S0_xy) + S0*((S0_xxy) + (S0_yyy))) - S0^3*((S0_xxxx) + 2*(S0_xxyy) + (S0_yyyy)))/(2*S0^8) - 3*((fx1_x) + (fy1_y))/(2*S0^2)-3*a3*(S0_t)/S0#-3//2 * (fx1_x + fy1_y)#
            #a3_t[1,i,j]  = (-6*S0^7*a3*(S0_t) + 20*(S0_x)^4 + 20*(S0_y)^4 - 16*S0*(S0_x)^2*(2*(S0_xx) + (S0_yy)) + S0^2*(5*(S0_xx)^2 + 4*(S0_xy)^2 + 8*(S0_x)*((S0_xxx) + (S0_xyy)) + 6*(S0_xx)*(S0_yy) + 5*(S0_yy)^2) - 8*(S0_y)^2*(-5*(S0_x)^2 + 2*S0*((S0_xx) + 2*(S0_yy))) + 8*S0*(S0_y)*(-4*(S0_x)*(S0_xy) + S0*((S0_xxy) + (S0_yyy))) - S0^3*((S0_xxxx) + 2*(S0_xxyy) + (S0_yyyy)) - 3*S0^6*((fx1_x) + (fy1_y)))/(2*S0^8)

            #fx1_t[1,i,j] = (-6*fx1*(S0_t) - 6*b13*(S0_x) + 6*g3*(S0_y) - S0*(a3_x) + 3*S0*(g3_y) - 3*S0*(b13_x))/(3*S0)#g3_y - b13_x - 1//3 * a3_x #
            #fy1_t[1,i,j] = (-6*fy1*(S0_t) + 6*g3*(S0_x) + 6*b13*(S0_y) - S0*(a3_y) + 3*S0*(g3_x) + 3*S0*(b13_y))/(3*S0)#g3_x + b13_y - 1//3 * a3_y #
            fx1_t[1,i,j] = (-6*fx1*(S0_t) - 6*b13*(S0_x) + 6*g3*(S0_y))/(3*S0) - (a3_x)/3 + (g3_y) - (b13_x)#g3_y - b13_x - 1//3 * a3_x #
            fy1_t[1,i,j] = (-6*fy1*(S0_t) + 6*g3*(S0_x) + 6*b13*(S0_y))/(3*S0) - (a3_y)/3 + (g3_x) + (b13_y)#g3_x + b13_y - 1//3 * a3_y #
        end
    end
	#println("Max of dxfx = ", maximum(fx1_x_vals))
	#println("Min of dxfx = ", minimum(fx1_x_vals))
	#println("-------------------------------------------------------")
	#println("Max of no source terms = ", maximum(no_source_terms))
	#println("Min of no source terms = ", minimum(no_source_terms))
	println("-------------------------------------------------------")
	println("t = ", test)
	println("iteration : ", step)
    output_dir = "/home/giulio/University/PhD/Simulations/testoutput/"
    if step % 100 == 0
        fx_t_filename = joinpath(output_dir, "fx_t$(step).h5")
        h5write(fx_t_filename, "fx1_t", fx1_t)
        fx_filename = joinpath(output_dir, "fx$(step).h5")
        h5write(fx_filename, "fx1", boundary.fx1)
        fy_filename = joinpath(output_dir, "fy$(step).h5")
        h5write(fy_filename, "fy1", boundary.fy1)
        a_filename = joinpath(output_dir, "a$(step).h5")
        h5write(a_filename, "a3", boundary.a3)
        g_filename = joinpath(output_dir, "g$(step).h5")
        h5write(g_filename, "g3", bulk.G)
        b_filename = joinpath(output_dir, "b$(step).h5")
        h5write(b_filename, "b3", bulk.B)
        
        ax_filename = joinpath(output_dir, "ax$(step).h5")
        h5write(ax_filename, "a3x", a3_x_vals)
        bx_filename = joinpath(output_dir, "bx$(step).h5")
        h5write(bx_filename, "b3x", b3_x_vals)
        gy_filename = joinpath(output_dir, "gy$(step).h5")
        h5write(gy_filename, "g3y", g3_y_vals)
        source_filename = joinpath(output_dir, "sigma$(step).h5")
        h5write(source_filename, "sigma", source_vals)
        
    end
	
	
    nothing
end
