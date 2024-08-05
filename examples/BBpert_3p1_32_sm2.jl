using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.2,
    u_outer_domains  =  3,
    u_outer_nodes    =  28,
    u_inner_nodes    =  12,
    fd_order         =  2,
    sigma_diss       =  0.01,
)

id = BlackBranePert(
    #B_amp  = 0.1,
    energy_dens  = 2.0,
    #B_ny = 1.0,
    #a3_ampx = 0.1,
    #a3_kx  = 2,
    a3_ampy = 0.05,
    #a3_ky  = 10,
    xmax = grid.x_max,
    xmin = grid.x_min,
    ymax = grid.y_max,
    ymin = grid.y_min,
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    #out_boundary_every   = 1,
    #out_bulk_every        = 1,
    out_bulkconstrained_every = 40,
    #out_gauge_every      = 10,
    remove_existing     = true,
)

integration = Integration(
    #dt              = 0.001,
    tmax            = 10.0,
    ODE_method      = AdS4_3_1.VCABM3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
