using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.01,
    u_outer_domains  =  3,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS5_3_1.BoostedBBnumerical(
    AH_pos = 1.0,
    #IDdir = "/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/",
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.00),
    
)

io = InOut(
    #out_boundary_every  = 10,
    #out_bulk_every      = 10,
    out_bulkconstrained_every = 1,
    #out_gauge_every     = 1,
    remove_existing     = true,
    checkpoint_every_walltime_hours = 0.1,
    #recover                     = :yes,
    #recover_dir                 = "/home/giulio/University/PhD/JeccoGt2/examples/FullBBB_3p1"
)

integration = Integration(

    dt              = 0.001,
    tmax            = 1,
    ODE_method      = AdS5_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
