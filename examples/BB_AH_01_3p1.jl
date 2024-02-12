
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  24,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  24,
    u_outer_min      =  0.1,
    u_outer_max      =  1.01,
    u_outer_domains  =  4,
    u_outer_nodes    =  36,
    u_inner_nodes    =  36,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

id   = BlackBrane(
   AH_pos = 1.0,
   #AH_pos = 0.9,
   energy_dens = 0.6,
)

evoleq = AffineNull(
    #phi0       = 0.0,
    #potential  = ZeroPotential(),
    #gaugecondition = ConstantAH(u_AH = 1.0),
    #gaugecondition = ConstantAH(u_AH = 0.5),
)

io = InOut(
    out_boundary_every  = 10,
    out_bulk_every      = 10,
    out_bulkconstrained_every = 1,
    #out_gauge_every     = 10,
    remove_existing     = true,
)

integration = Integration(
    dt              = 0.001,
    tmax            = 5.0,
    #ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
