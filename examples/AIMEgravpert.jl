using Jecco.AdS5_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.1,
    u_outer_max      =  1.1,
    u_outer_domains  =  1,
    u_outer_nodes    =  80,
    u_inner_nodes    =  24,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS5_3_1.BoostedBBnumerical(
    AH_pos = 1.0,
    
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.00),# this u_AH is a not used parameter
    
)

io = InOut(
    #out_boundary_every  = 3,
    out_bulk_every      = 20,
    #out_bulkconstrained_every = 50,
    #out_gauge_every     = 1,
    remove_existing     = true,
    #checkpoint_every_walltime_hours = 0.1,
    #recover                     = :yes,
    #recover_dir                 = "/home/giulio/University/PhD/JeccoGt2/examples/FullBBB_3p1"
)

integration = Integration(

    dt              = 0.001,
    tmax            = 50,
    ODE_method      = AdS5_3_1.VCABM3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
