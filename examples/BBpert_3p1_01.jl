
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  64,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  64,
    u_outer_min      =  0.1,
    u_outer_max      =  1.003,
    u_outer_domains  =  2,
    u_outer_nodes    =  12,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

id = BlackBranePert(
    B_amp  = 0.1,
    energy_dens  = 1.5,
    B_ny = 1.0,
    #AH_pos = 0.9,
    #a3_ampx = 0.1,
    #a3_kx  = 2,
    #a3_ampy = 0.05,
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
    out_bulk_every        = 1,
    #out_constrained_every = 1,
    #out_gauge_every      = 10,
)

integration = Integration(
    dt              = 0.002,
    tmax            = 10.0,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
