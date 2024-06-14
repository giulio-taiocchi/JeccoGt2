
using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -1.0,
    x_max            =  1.0,
    x_nodes          =  16,
    y_min            = -1.0,
    y_max            =  1.0,
    y_nodes          =  16,
    u_outer_min      =  0.1,
    u_outer_max      =  1.1,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

id = BlackBranePert(
    energy_dens = 0.75,
    a4_ampx  = 0.05,
    #a4_kx    = 1,
    AH_pos   = 1.05,
    xmax     = grid.x_max,
    xmin     = grid.x_min,
    ymin     = grid.y_min,
    ymax     = grid.y_max,
)

evoleq = AffineNull(
    phi0           = 0.0,
    potential      = ZeroPotential(),
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    #out_boundary_every  = 200,
    out_bulkconstrained_every = 50,
    #out_bulk_every      = 1000,
    #out_gauge_every     = 200,
    remove_existing     = true,
)

integration = Integration(
    dt              = 0.001,
    tmax            = 20.0,
    ODE_method      = AdS5_3_1.AB4(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
