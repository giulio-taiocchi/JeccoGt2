using Jecco.AdS4_3_1_s

grid = SpecCartGrid3D(
    x_min            = -100000.0,
    x_max            =  100000.0,
    x_nodes          =  128,
    y_min            = -100000.0,
    y_max            =  100000.0,
    y_nodes          =  128,
    u_outer_min      =  0.1,
    u_outer_max      =  1.03,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)



id = AdS4_3_1_s.BoostedBBnumerical(
    AH_pos = 0.9975136584870653,
    IDdir = "/home/giulio/University/PhD/Initial_data/sourced/",
)


evoleq = AffineNull(
    source = GaussianSource(time = 0.0, sigmax=0.4,sigmay=0.4, x0=0.0, y0=0.0,Amp=2.0,t0=-2.0, L=200000.0),
    #source = NoSource(),
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    #out_boundary_every  = 10,
    out_bulkconstrained_every = 200,
    #out_bulk_every      = 10,
    #out_gauge_every     = 1,
    remove_existing     = true,
    checkpoint_every_walltime_hours = 0.1,
    #recover                     = :yes,
    #recover_dir                 = "/home/giulio/University/PhD/JeccoGt2/examples/BB_AH_01_4p1_sA0d25XY300L500"
)

integration = Integration(
    #dt              = 0.001,
    tmax            = 20.0,
    #ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
