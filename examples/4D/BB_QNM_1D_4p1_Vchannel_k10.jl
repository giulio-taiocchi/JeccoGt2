
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  64,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  64,
    u_outer_min      =  0.1,
    u_outer_max      =  1.01,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS4_3_1.BlackBranePert(
    energy_dens  = 2.0,
    fx1_ampy = 0.1,
    fx1_ky = 2,
    xmax = 5.0,
    xmin = -5.0,
    ymax = 5.0,
    ymin = -5.0,
     #AH_pos = 1.0,
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    out_boundary_every  = 200,
    #out_bulkconstrained_every = 100,
    #out_bulk_every      = 100,
    #out_gauge_every     = 10,
    remove_existing     = true,
)

integration = Integration(
    dt              = 0.003,
    tmax            = 250.0,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
