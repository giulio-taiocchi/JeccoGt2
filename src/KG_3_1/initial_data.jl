
Base.@kwdef struct Uniform2D{T} <: InitialData
    phi2  :: T  = 1.0
end

Base.@kwdef struct Sine2D{T} <: InitialData
    Lx :: T
    Ly :: T
    kx :: Int
    ky :: Int
end

Base.@kwdef struct Uniform{T} <: InitialData
    phi2  :: T  = 1.0
end

Base.@kwdef struct Polinomial{T} <: InitialData
    u0 :: T
    u1 :: T
    u2 :: T
    u3 :: T
end

function (id::InitialData)(bulkevols, boundary::Boundary, systems::SystemPartition)
    init_data!(boundary, systems[1], id)
    init_data!(bulkevols, systems, id)
    nothing
end

function init_data!(bulkevols, systems::SystemPartition, id::InitialData)
    # the Ref() makes its argument a scalar with respect to broadcast
    init_data!.(bulkevols, systems, Ref(id))
end

function init_data!(bulk::BulkEvolved, sys::System, id::InitialData)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    phi = getphi(bulk)

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u  = uu[a]
                x  = xx[i]
                y  = yy[j]
                phi[a,i,j] = analytic_phi(u, x, y, id)
            end
        end
    end

    bulk
end

analytic_phi(u, x, y, id::Uniform2D) = id.phi2

analytic_phi(u, x, y, id::Uniform) = id.phi2

analytic_phi(u, x, y, id::Polinomial) = id.u0 + u * id.u1 + u^2 * id.u2 + u^3 * id.u3

analytic_phi(u, x, y, id::Sine2D) =
   sin( 2*π * id.kx / id.Lx * x ) * sin( 2*π * id.ky / id.Ly * y )
   



function init_data!(boundary::Boundary, sys::System, id::InitialData)
    a4GF = geta4(boundary)

    fill!(a4GF, -1)
end
