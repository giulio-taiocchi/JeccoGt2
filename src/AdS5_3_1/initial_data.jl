Base.@kwdef struct BlackBrane_xi1{T} <: InitialData
    a40           :: T   = 1.0
    # guess for the AH position
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    xi_0          :: T   = 0.0
    xi_Ax         :: T   = 0.0
    xi_nx         :: Int = 1
    xi_Ay         :: T   = 0.0
    xi_ny         :: Int = 1
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    ahf           :: AHF = AHF()
end


abstract type ID_ConstantAH  <: InitialData end

Base.@kwdef struct BlackBrane{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    ahf           :: AHF = AHF()
end



Base.@kwdef struct BlackBraneNumericalphi{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    ahf           :: AHF = AHF()
end

Base.@kwdef struct BoostedBBnumerical{T} <: ID_ConstantAH
    #energy_dens :: T   = 5.0
    AH_pos      :: T   = 1.0
    phi0        :: T   = 0.0
    ahf         :: AHF = AHF()
    IDdir	:: AbstractString = "/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/"
end

Base.@kwdef struct BlackBranePert{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    oophiM2       :: T   = 0.0
    B1_amp        :: T   = 0.0
    B1_nx         :: Int = 1
    B1_ny         :: Int = 2
    B2_amp        :: T   = 0.0
    B2_nx         :: Int = 1
    B2_ny         :: Int = 2
    G_amp         :: T   = 0.0
    G_nx          :: Int = 1
    G_ny          :: Int = 2
    phi2          :: T   = 0.0
    phi5          :: T   = 0.0
    a4_ampx       :: T   = 0.0
    a4_ampy       :: T   = 0.0
    a4_kx         :: Int = 1
    a4_ky         :: Int = 1
    fx2_ampx      :: T   = 0.0
    fx2_ampy      :: T   = 0.0
    fx2_kx        :: Int = 1
    fx2_ky        :: Int = 1
    fy2_ampx      :: T   = 0.0
    fy2_ampy      :: T   = 0.0
    fy2_kx        :: Int = 1
    fy2_ky        :: Int = 1
    xi0           :: T   = 0.0
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    ahf           :: AHF = AHF()
end

Base.@kwdef struct PhiGaussian_u{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    phi0          :: T   = 0.0
    phi2          :: T   = 0.0
    oophiM2       :: T   = 0.0
    amp           :: T   = 0.0
    u0            :: T   = 0.0
    sigma         :: T   = 0.1
    ahf           :: AHF = AHF()
end

Base.@kwdef struct QNM_1D{T} <: InitialData
    energy_dens :: T   = 1.0
    phi0        :: T   = 0.0
    phi2        :: T   = 0.0
    oophiM2     :: T   = 0.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
end

Base.@kwdef struct BlackBraneGaussPert{T} <: InitialData
    energy_dens   :: T     = 1.0
    AH_pos        :: T     = 1.0
    phi0          :: T     = 0.0
    oophiM2       :: T     = 0.0
    phi2          :: T     = 0.0
    xi0           :: T     = 0.0
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    sigma         :: T
    ahf           :: AHF    = AHF()
    recover       :: Symbol = :no
    recover_dir   :: String

end

Base.@kwdef struct BlackBraneNoise{T} <: InitialData
    energy_dens   :: T      = 1.0
    AH_pos        :: T      = 1.0
    phi0          :: T      = 0.0
    oophiM2       :: T      = 0.0
    phi2          :: T      = 0.0
    xi0           :: T      = 0.0
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    pert          :: T
    nx_max        :: Int    = 1
    ny_max        :: Int    = 1
    ahf           :: AHF    = AHF()
end

function (id::InitialData)(bulkconstrains, bulkevols, bulkderivs, boundary::Boundary,
                           gauge::Gauge, horizoncache::HorizonCache, systems::SystemPartition,
                           evoleq::EvolutionEquations)
    _, Nx, Ny = size(systems[end])
    AH_pos    = id.AH_pos
    xi        = getxi(gauge)

    # function to solve the nested system
    nested = Nested(systems, bulkconstrains)

    init_data!(boundary, systems[1],   id)

    init_data!(gauge,    systems[end], id)
    init_data!(bulkevols, gauge, systems, id)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    # find the Apparent Horizon
    sigma = similar(gauge.xi)
    fill!(sigma, 1/AH_pos)  # initial guess
    find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
             horizoncache, systems[end], id.ahf)

    nothing
end



function (id::ID_ConstantAH)(bulkconstrains, bulkevols, bulkderivs, boundary::Boundary,
                           gauge::Gauge, horizoncache::HorizonCache, systems::SystemPartition,
                           evoleq::EvolutionEquations)
    _, Nx, Ny = size(systems[end])
    AH_pos    = id.AH_pos
    xi        = getxi(gauge)

    # function to solve the nested system
    nested = Nested(systems, bulkconstrains)

    init_data!(boundary, systems[1],   id)

    init_data!(gauge,    systems[end], id)
    init_data!(bulkevols, gauge, systems, id)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    # find the Apparent Horizon
    sigma = similar(gauge.xi)
    fill!(sigma, 1/AH_pos)  # initial guess
    find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
             horizoncache, systems[end], id.ahf)

    # assuming that the AH has been found, we now update xi and the bulk variables
    AH_fixed_pos=evoleq.gaugecondition.u_AH
    for j in 1:Ny
        for i in 1:Nx
            xi[1,i,j] += -1 / AH_fixed_pos + sigma[1,i,j]
        end
    end

    init_data!(bulkevols, gauge, systems, id)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    # AH should now be at u = AH_pos

    nothing
end

function init_data!(bulkevols, gauge::Gauge, systems::SystemPartition,
                    id::InitialData)
    # the Ref() makes its argument a scalar with respect to broadcast
    #init_data!.(bulkevols, Ref(gauge), systems, Ref(id))
    counting = 0
    for k in systems
    	init_data!.(Ref(bulkevols.x[counting+1]), Ref(gauge), Ref(k), Ref(id),Ref(counting))
    	counting = counting + 1
    end
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Inner},
                    id::InitialData,counting)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    B1  = getB1(bulk)
    B2  = getB2(bulk)
    G   = getG(bulk)
    phi = getphi(bulk)
    xi  = getxi(gauge)

    phi0 = id.phi0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u       = uu[a]
                x       = xx[i]
                y       = yy[j]
                xi_ij   = xi[1,i,j]
                aux     = 1 + xi_ij * u
                aux3    = aux * aux * aux
                aux4    = aux * aux3
                u_old   = u / aux
                B1_old  = analytic_B1(a, i, j, u_old, x, y, id, counting)
                B2_old  = analytic_B2(a, i, j, u_old, x, y, id, counting)
                G_old   = analytic_G(a, i, j, u_old, x, y, id, counting)

                B1[a,i,j]  = B1_old / aux4
                B2[a,i,j]  = B2_old / aux4
                G[a,i,j]   = G_old  / aux4
            end
        end
    end

    # if phi0 = 0 set phi to zero
    if abs(phi0) < 1e-9
        fill!(phi, 0)
    else
        for j in 1:Ny
            for i in 1:Nx
                for a in 1:Nu
                    u       = uu[a]
                    x       = xx[i]
                    y       = yy[j]
                    xi_ij   = xi[1,i,j]
                    aux     = 1 + xi_ij * u
                    aux3    = aux * aux * aux
                    aux4    = aux * aux3
                    u_old   = u / aux
                    phi_old = analytic_phi(a, i, j, u_old, x, y, id, counting)

                    phi[a,i,j] = xi_ij * xi_ij / (phi0 * phi0 * aux) +
                        phi_old / aux3
                end
            end
        end
    end

    bulk
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Outer},
                    id::InitialData,counting)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    B1  = getB1(bulk)
    B2  = getB2(bulk)
    G   = getG(bulk)
    phi = getphi(bulk)
    xi  = getxi(gauge)

    phi0 = id.phi0

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u         = uu[a]
                x         = xx[i]
                y         = yy[j]
                xi_ij     = xi[1,i,j]
                aux       = 1 + xi_ij * u
                aux3      = aux * aux * aux
                aux4      = aux * aux3
                u_old     = u / aux
                B1_old    = analytic_B1(a, i, j, u_old, x, y, id, counting)
                B2_old    = analytic_B2(a, i, j, u_old, x, y, id, counting)
                G_old     = analytic_G(a, i, j, u_old, x, y, id, counting)
                B1_inner  = B1_old / aux4
                B2_inner  = B2_old / aux4
                G_inner   = G_old  / aux4

                B1[a,i,j]  = u^4 * B1_inner
                B2[a,i,j]  = u^4 * B2_inner
                G[a,i,j]   = u^4 * G_inner
            end
        end
    end

    # if phi0 = 0 set phi to zero
    if abs(phi0) < 1e-9
        fill!(phi, 0)
    else
        for j in 1:Ny
            for i in 1:Nx
                for a in 1:Nu
                    u         = uu[a]
                    x         = xx[i]
                    y         = yy[j]
                    xi_ij     = xi[1,i,j]
                    aux       = 1 + xi_ij * u
                    aux3      = aux * aux * aux
                    aux4      = aux * aux3
                    u_old     = u / aux
                    phi_old   = analytic_phi(a, i, j, u_old, x, y, id, counting)
                    phi_inner = xi_ij * xi_ij / (phi0 * phi0 * aux) +
                        phi_old / aux3

                    phi[a,i,j] = u * phi0 - u^2 * phi0 * xi_ij + u^3 * phi0^3 * phi_inner
                end
            end
        end
    end

    bulk
end


# BlackBrane_xi1

analytic_B1(i, j, k, u, x, y,  id::BlackBrane_xi1, whichsystem)  = 0
analytic_B2(i, j, k, u, x, y, id::BlackBrane_xi1, whichsystem)  = 0
analytic_G(i, j, k, u, x, y,  id::BlackBrane_xi1, whichsystem)   = 0
analytic_phi(i, j, k, u, x, y, id::BlackBrane_xi1, whichsystem) = 0

function init_data!(ff::Boundary, sys::System, id::BlackBrane_xi1)
    a40 = id.a40

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBrane_xi1)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    xi  = getxi(ff)

    xmax  = id.xmax
    xmin  = id.xmin
    ymax  = id.ymax
    ymin  = id.ymin

    for j in 1:Ny
        for i in 1:Nx
            x         = xx[i]
            y         = yy[j]

            xi[1,i,j] = id.xi_0 +
                id.xi_Ax * sin( 2 * π * id.xi_nx * (xmax-x)/(xmax-xmin) ) +
                id.xi_Ay * sin( 2 * π * id.xi_ny * (ymax-y)/(ymax-ymin) )
        end
    end

    ff
end


# BlackBrane initial data

analytic_B1(i, j, k, u, x, y,  id::BlackBrane, whichsystem)  = 0
analytic_B2(i, j, k, u, x, y,  id::BlackBrane, whichsystem)  = 0
analytic_G(i, j, k, u, x, y,  id::BlackBrane, whichsystem)   = 0
analytic_phi(i, j, k, u, x, y,  id::BlackBrane, whichsystem) = 0

function init_data!(ff::Boundary, sys::System, id::BlackBrane)
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBrane)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# BlackBranePert initial data

analytic_phi(i, j, k, u, x, y,  id::BlackBranePert, whichsystem) = id.phi2 / id.phi0^3 + id.phi5 / id.phi0^3*u

function analytic_B1(i, j, k, u, x, y,  id::BlackBranePert, whichsystem)
    # add the perturbation on B1
    pert_amp = id.B1_amp
    xmax     = id.xmax
    xmin     = id.xmin
    ymax     = id.ymax
    ymin     = id.ymin
    # number of maxima in each direction
    nx       = id.B1_nx
    ny       = id.B1_ny

    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
        sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end

function analytic_B2(i, j, k, u, x, y, id::BlackBranePert, whichsystem)
    # add the perturbation on B2
    pert_amp = id.B2_amp
    xmax     = id.xmax
    xmin     = id.xmin
    ymax     = id.ymax
    ymin     = id.ymin
    # number of maxima in each direction
    nx       = id.B2_nx
    ny       = id.B2_ny

    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
        sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end

function analytic_G(i, j, k, u, x, y,  id::BlackBranePert, whichsystem)
    # add the perturbation on G
    pert_amp = id.G_amp
    xmax     = id.xmax
    xmin     = id.xmin
    ymax     = id.ymax
    ymin     = id.ymin
    # number of maxima in each direction
    nx       = id.G_nx
    ny       = id.G_ny

    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
        sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end

function init_data!(ff::Boundary, sys::System{Inner}, id::BlackBranePert)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.oophiM2

    # a4 perturbation amplitude
    ampx     = id.a4_ampx
    ampy     = id.a4_ampy
    fx2_ampx = id.fx2_ampx
    fx2_ampy = id.fx2_ampy
    fy2_ampx = id.fy2_ampx
    fy2_ampy = id.fy2_ampy
    # number of maxima
    kx     = id.a4_kx
    ky     = id.a4_ky
    fx2_kx = id.fx2_kx
    fx2_ky = id.fx2_ky
    fy2_kx = id.fy2_kx
    fy2_ky = id.fy2_ky

    xmax = id.xmax
    xmin = id.xmin
    xmid = (xmax + xmin) / 2
    ymax = id.ymax
    ymin = id.ymin
    ymid = (ymax + ymin) / 2

    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord

    phi04 = phi0 * phi0 * phi0 * phi0
    a40   = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    for j in 1:Ny
        for i in 1:Nx
            x = xx[i]
            y = yy[j]
            a4[1,i,j]  += -a40 * ( ampx * cos(2 * π * kx * (x-xmid)/(xmax-xmin)) +
                                   ampy * cos(2 * π * ky * (y-ymid)/(ymax-ymin)) )

            fx2[1,i,j] += fx2_ampx * cos(2 * π * fx2_kx * (x-xmid)/(xmax-xmin)) +
                           fx2_ampy * cos(2 * π * fx2_ky * (y-ymid)/(ymax-ymin))

            fy2[1,i,j] += fy2_ampx * cos(2 * π * fy2_kx * (x-xmid)/(xmax-xmin)) +
                          fy2_ampy * cos(2 * π * fy2_ky * (y-ymid)/(ymax-ymin))
        end
    end

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBranePert)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos

    # TODO: this guess works best for the conformal case. is there a better one?
    if id.xi0 == 0
        xi0 = (-a40)^0.25 - 1/AH_pos
    else
        xi0 = id.xi0
    end

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# PhiGaussian_u initial data

function analytic_phi(u, x, y, id::PhiGaussian_u)
    phi0   = id.phi0
    phi2   = id.phi2
    amp    = id.amp
    sigma  = id.sigma
    u0     = id.u0

    phi03  = phi0 * phi0 * phi0
    sigma2 = sigma * sigma

    myfunc = amp * exp( -(u-u0)*(u-u0) / (2*sigma2) )

    (phi2 + u * myfunc) / phi03
end

analytic_B1(u, x, y, id::PhiGaussian_u) = 0
analytic_B2(u, x, y, id::PhiGaussian_u) = 0
analytic_G(u, x, y, id::PhiGaussian_u)  = 0

function init_data!(ff::Boundary, sys::System, id::PhiGaussian_u)
    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::PhiGaussian_u)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    AH_pos  = id.AH_pos
    oophiM2 = id.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    # TODO: this guess works best for the conformal case. is there a better one?
    xi0 = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


#QNM in 1D initial data
analytic_B1(i, j, k, u, x, y,  id::QNM_1D, whichsystem)  = 3/2*0.1*u^8
analytic_B2(i, j, k, u, x, y,  id::QNM_1D, whichsystem)  = 1/2*0.1*u^8
analytic_G(i, j, k, u, x, y,  id::QNM_1D, whichsystem)   = 0
analytic_phi(i, j, k, u, x, y,  id::QNM_1D, whichsystem) = id.phi2/id.phi0^3

function init_data!(ff::Boundary, sys::System, id::QNM_1D)
    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::QNM_1D)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    AH_pos  = id.AH_pos
    oophiM2 = id.oophiM2
    phi04   = phi0 * phi0 * phi0 * phi0

    a40 = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    xi0 = 0

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


#Gaussian prob distributions for all the modes
analytic_phi(u, x, y, id::BlackBraneGaussPert) = id.phi2 / id.phi0^3
analytic_B1(u, x, y, id::BlackBraneGaussPert)  = 0
analytic_B2(u, x, y, id::BlackBraneGaussPert)  = 0
analytic_G(u, x, y, id::BlackBraneGaussPert)   = 0

function a4Perturbation(x::Real, y::Real, nx_max::Int, ny_max::Int, Lx::Real, Ly::Real,
                            a::Array{T,2}, b::Array{T,2}, c::Array{T,2}, d::Array{T,2}) where {T<:Real}

    sum  = 0.0
    for ny in 0:ny_max-1
        for nx in 0:nx_max-1
            sum += a[nx+1,ny+1] * cos(2 * π * nx * x / Lx) * cos(2 * π * ny * y / Ly) +
                    b[nx+1,ny+1] * cos(2 * π * nx * x / Lx) * sin(2 * π * ny * y / Ly) +
                    c[nx+1,ny+1] * sin(2 * π * nx * x / Lx) * cos(2 * π * ny * y / Ly) +
                    d[nx+1,ny+1] * sin(2 * π * nx * x / Lx) * sin(2 * π * ny * y / Ly)
        end
    end

    sum
end

function init_data!(ff::Boundary, sys::System{Inner}, id::BlackBraneGaussPert)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.oophiM2

    xmax = id.xmax
    xmin = id.xmin
    xmid = (xmax + xmin) / 2
    ymax = id.ymax
    ymin = id.ymin
    ymid = (ymax + ymin) / 2
    Lx   = xmax-xmin
    Ly   = ymax-ymin

    _, Nx, Ny = size(sys)
    xx        = sys.xcoord
    yy        = sys.ycoord

    phi04 = phi0 * phi0 * phi0 * phi0
    a40   = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    if id.recover == :no
        dist = Normal(0.0, id.sigma)
        a    = reshape(rand(dist, Nx*Ny), Nx, Ny)
        b    = reshape(rand(dist, Nx*Ny), Nx, Ny)
        c    = reshape(rand(dist, Nx*Ny), Nx, Ny)
        d    = reshape(rand(dist, Nx*Ny), Nx, Ny)
        Nkx  = Nx
        Nky  = Ny
    elseif id.recover == :yes
        ts         = OpenPMDTimeSeries(id.recover_dir, prefix="boundary_")
        a4old      = get_field(ts, it=ts.iterations[1], field="a4")[1][1,:,:]
        a, b, c, d = AdS5_3_1.Fourier_cos_sin(a4old)
        Nkx, Nky   = size(a)
        a[1,1]     = 0.
        if Nkx > Nx || Nky > Ny
            @warn "Too few points for such ammount of modes"
            return
        end
        a *= -3/4/epsilon
        b *= -3/4/epsilon
        c *= -3/4/epsilon
        d *= -3/4/epsilon
    else
        @warn "No such recover option"
        return
    end

    for j in 1:Ny
        for i in 1:Nx
            x = xx[i]
            y = yy[j]
            # I decided to change to δε/ε following a normal, so the normal pert has to be
            #multiplied by a factor to enter in the energy correctly.
            # Going to previous means 4/3*epsilon -> a40
            a4[1,i,j] += -4/3 * epsilon * a4Perturbation(x, y, Nkx, Nky, Lx, Ly, a, b, c, d)
        end
    end

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBraneGaussPert)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos

    # TODO: this guess works best for the conformal case. is there a better one?
    if id.xi0 == 0
        xi0 = (-a40)^0.25 - 1/AH_pos
    else
        xi0 = id.xi0
    end

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


#Same pert for all modes

function Perturbation(x::Real, y::Real, nx_max::Int, ny_max::Int, Lx::Real, Ly::Real)
    sum  = 0.0
    @inbounds for ny in 0:ny_max
        for nx in 0:nx_max
            if nx != 0 || ny != 0
                mode1 = cos(2 * π * nx * x / Lx) * cos(2 * π * ny * y / Ly)
                mode2 = cos(2 * π * nx * x / Lx) * sin(2 * π * ny * y / Ly)
                mode3 = sin(2 * π * nx * x / Lx) * cos(2 * π * ny * y / Ly)
                mode4 = sin(2 * π * nx * x / Lx) * sin(2 * π * ny * y / Ly)
                sum += mode1+mode2+mode3+mode4
            end
        end
    end
    sum
end

# analytic_phi(u, x, y, id::BlackBraneNoise) = id.phi2 / id.phi0^3 + id.pert * Perturbation(x, y, id.nx_max, id.ny_max, id.xmax-id.xmin, id.ymax-id.ymin)
# analytic_B1(u, x, y, id::BlackBraneNoise)  = id.pert * Perturbation(x, y, id.nx_max, id.ny_max, id.xmax-id.xmin, id.ymax-id.ymin)
# analytic_B2(u, x, y, id::BlackBraneNoise)  = id.pert * Perturbation(x, y, id.nx_max, id.ny_max, id.xmax-id.xmin, id.ymax-id.ymin)
# analytic_G(u, x, y, id::BlackBraneNoise)   = id.pert * Perturbation(x, y, id.nx_max, id.ny_max, id.xmax-id.xmin, id.ymax-id.ymin)
analytic_phi(u, x, y, id::BlackBraneNoise) = id.phi2 / id.phi0^3
analytic_B1(u, x, y, id::BlackBraneNoise)  = 0.0
analytic_B2(u, x, y, id::BlackBraneNoise)  = 0.0
analytic_G(u, x, y, id::BlackBraneNoise)   = 0.0

function init_data!(ff::Boundary, sys::System{Inner}, id::BlackBraneNoise)
    epsilon = id.energy_dens
    phi0    = id.phi0
    phi2    = id.phi2
    oophiM2 = id.oophiM2

    xmax = id.xmax
    xmin = id.xmin
    xmid = (xmax + xmin) / 2
    ymax = id.ymax
    ymin = id.ymin
    ymid = (ymax + ymin) / 2
    Lx   = xmax-xmin
    Ly   = ymax-ymin

    _, Nx, Ny = size(sys)
    xx        = sys.xcoord
    yy        = sys.ycoord

    phi04  = phi0 * phi0 * phi0 * phi0
    a40    = (-epsilon - phi0 * phi2 - phi04 * oophiM2 / 4 + 7 * phi04 / 36) / 0.75
    pert   = id.pert
    nx_max = id.nx_max
    ny_max = id.ny_max

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    for j in 1:Ny
        for i in 1:Nx
            x = xx[i]
            y = yy[j]
            a4[1,i,j] += -4/3*pert * Perturbation(x, y, nx_max, ny_max, Lx, Ly)
        end
    end

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBraneNoise)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos

    # TODO: this guess works best for the conformal case. is there a better one?
    if id.xi0 == 0
        xi0 = (-a40)^0.25 - 1/AH_pos
    else
        xi0 = id.xi0
    end

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# analytic BlackBrane with  numerical initial data for the potential 

analytic_B1(u, x, y, id::BlackBraneNumericalphi)  = 0
analytic_B2(u, x, y, id::BlackBraneNumericalphi)  = 0
analytic_G(u, x, y, id::BlackBraneNumericalphi)   = 0
#analytic_phi(u, x, y, id::BlackBrane) = 0
function analytic_phi(i, j, k, u, x, y, id::BlackBraneNumericalphi, whichsystem)
	uu = u
	initialphi=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/Initialphi_BB.h5")
	system_index = string(whichsystem+1)
	dset=initialphi[system_index]
	phi=read(dset)
	# here the indecex have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	phivalue = phi[k,j,i]
	#println("phi in u=$uu index: $i is $phivalue")
	phivalue
end

function init_data!(ff::Boundary, sys::System, id::BlackBraneNumericalphi)
    a40 = -id.energy_dens/0.75

    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, a40)
    fill!(fx2, 0)
    fill!(fy2, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBraneNumericalphi)
    a40     = -id.energy_dens/0.75
    AH_pos  = id.AH_pos
    xi0     = (-a40)^0.25 - 1/AH_pos

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# Numerical boosted Black Brane
analytic_phi(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem) = 0
#analytic_B1(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)  = 0
#analytic_B2(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)  = 0

function analytic_B1(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)
	uu = u
	dir = id.IDdir
	Bdirectory = dir*"InitialB_BBB.h5"
	initialB=h5open(Bdirectory)
	system_index = string(whichsystem+1)
	dset=initialB[system_index]
	B=read(dset)
	# here the indexes have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	Bvalue = B[k,j,i]
	
	
	Bvalue
end

function analytic_B2(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)
	uu = u
	dir = id.IDdir
	B2directory = dir*"InitialB2_BBB.h5"
	initialB2=h5open(B2directory)
	system_index = string(whichsystem+1)
	dset=initialB2[system_index]
	B2=read(dset)
	# here the indexes have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	B2value = B2[k,j,i]
	
	
	B2value
end

#function analytic_B2(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)
#	uu = u
#	dir = id.IDdir
#	initialB2=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/InitialB2_BBB.h5")
#	system_index = string(whichsystem+1)
#	dset=initialB2[system_index]
#	B2=read(dset)
	# here the indecex have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
#	Bvalue2 = B2[k,j,i]	
#	Bvalue2
#end

#function analytic_G(i, j, k, u, x, y, id::BoostedBBnumerical,whichsystem)  
	#uu = u
	#indexi = i
	#if j==1 
	#	if k==1
	#		#println("index i $indexi, u value:$uu")
	#	end
	#end
	#initialG=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/InitialG_BBB.h5")
	#system_index = string(whichsystem+1)
	#dset=initialG[system_index]
	#G=read(dset)
	# here the indecex have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be G[i,j,k]
	#Gvalue = G[k,j,i]
	#Gvalue
#end

function analytic_G(i, j, k, u, x, y, id::BoostedBBnumerical,whichsystem) 
	uu = u
	dir = id.IDdir
	Gdirectory = dir*"InitialG_BBB.h5"
	initialG=h5open(Gdirectory)
	system_index = string(whichsystem+1)
	dset=initialG[system_index]
	G=read(dset)
	# here the indexes have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	Gvalue = G[k,j,i]
		
	Gvalue
end



function init_data!(ff::Boundary, sys::System, id::BoostedBBnumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    
    a4  = geta4(ff)
    fx2 = getfx2(ff)
    fy2 = getfy2(ff)

    fill!(a4, 0)
    fill!(fx2, 0)
    fill!(fy2, 0)
    
    #a4data = h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/Initiala3_BBB.h5")
    #fx2data = h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/Initialfx_BBB.h5")
    #fy2data = h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/Initialfy_BBB.h5")
    dir = id.IDdir
    a4directory = dir*"Initiala3_BBB.h5"
    fx2directory = dir*"Initialfx_BBB.h5"
    fy2directory = dir*"Initialfy_BBB.h5"
    a4data = h5open(a4directory)
    fx2data = h5open(fx2directory)
    fy2data = h5open(fy2directory)
    a4read = read(a4data["a3"])
    fx2read = read(fx2data["fx"])
    fy2read = read(fy2data["fy"])
    for j in 1:Ny
        for i in 1:Nx      
                x = xx[i]
                y = yy[j]         
                a4[1,i,j] = a4read[j,i]
                fx2[1,i,j] = fx2read[j,i]
                fy2[1,i,j] = fy2read[j,i]
        end
    end
    ff
end

function init_data!(ff::Gauge, sys::System, id::BoostedBBnumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    AH_pos  = id.AH_pos
    
    xi  = getxi(ff)
    fill!(xi, 0)
    dir = id.IDdir
    xidirectory = dir*"Initialxi_BBB.h5"
    xidata = h5open(xidirectory)
    xiread = read(xidata["xi"])
    for j in 1:Ny
        for i in 1:Nx      
                x = xx[i]
                y = yy[j]         
                xi[1,i,j] = xiread[j,i]                
        end
    end
    
    #for j in 1:Ny
    #    for i in 1:Nx      
    #            x = xx[i]
    #            y = yy[j]         
    #            xi[1,i,j] = 1/4*(1-cos(2*π*x)*cos(2*π*x))
    #    end
    #end
    ff
end


