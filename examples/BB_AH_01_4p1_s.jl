using Jecco.AdS4_3_1_s

grid = SpecCartGrid3D(
    x_min            = -50.0,
    x_max            =  50.0,
    x_nodes          =  150,
    y_min            = -50.0,
    y_max            =  50.0,
    y_nodes          =  150,
    u_outer_min      =  0.1,
    u_outer_max      =  1.05,
    u_outer_domains  =  7,
    u_outer_nodes    =  12,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)



id = AdS4_3_1_s.BoostedBBnumerical(
    AH_pos = 1.0,
    a3_ampx	= 0.,
    a3_translx	= -1.375,
    A = 0.559017,
    B = 0.0,
    phase_a = 0.,
    phase_fx = 0.,
)


evoleq = AffineNull(
    source = GaussianSource(time = 0.0, sigmax=0.015625,sigmay=0.015625,Amp=1.0,t0=0.8),
    gaugecondition = ConstantAH(u_AH = 1.0),
    #gaugecondition = ConstantAH(u_AH = 0.5),
)

io = InOut(
    out_boundary_every  = 100,
    out_bulkconstrained_every = 25,
    out_bulk_every      = 500,
    #out_gauge_every     = 100,
    remove_existing     = true,
    checkpoint_every_walltime_hours = 1,
)

integration = Integration(
    dt              = 0.001,
    tmax            = 50.0,
    #ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
