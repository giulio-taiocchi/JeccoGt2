using Jecco.AdS4_3_1_s

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  96,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  96,
    u_outer_min      =  0.1,
    u_outer_max      =  1.01,
    u_outer_domains  =  8,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)



id = AdS4_3_1_s.BoostedBBnumerical(
    AH_pos = 1.0,
    a3_ampx	= 0.,
    a3_translx	= -5/2,
    A = 1.,
    B = 1.0,
    phase_a = 0.,
    phase_fx = 0.,
)


evoleq = AffineNull(
    source = GaussianSource(time = 0.0, sigmax=1,sigmay=1),
    gaugecondition = ConstantAH(u_AH = 1.0),
    #gaugecondition = ConstantAH(u_AH = 0.5),
)

io = InOut(
    out_boundary_every  = 1,
    out_bulkconstrained_every = 1,
    out_bulk_every      = 1,
    #out_gauge_every     = 10,
    remove_existing     = true,
)

integration = Integration(
    dt              = 0.001,
    tmax            = 50.0,
    #ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
