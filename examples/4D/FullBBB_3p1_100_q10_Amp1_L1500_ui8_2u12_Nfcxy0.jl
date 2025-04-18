using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -750.0,
    x_max            =  750.0,
    x_nodes          =  100,
    y_min            = -750.0,
    y_max            =  750.0,
    y_nodes          =  100,
    u_outer_min      =  0.1,
    u_outer_max      =  1.2,
    u_outer_domains  =  2,
    u_outer_nodes    =  12,
    u_inner_nodes    =  8,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS4_3_1.BoostedBBnumerical(
    #AH_pos = 0.8338097372270935,
    IDdir = "/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/",
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.00),
    
)

io = InOut(
    out_boundary_every  = 200,
    #out_bulk_every      = 1,
    #out_bulkconstrained_every = 1,
    out_gauge_every     = 100,
    remove_existing     = true,
    checkpoint_every_walltime_hours = 0.1,
    #recover                     = :yes,
    #recover_dir                 = "/home/giulio/University/PhD/JeccoGt2/examples/FullBBB_3p1"
)

integration = Integration(

    #dt              = 0.0005,
    tmax            = 200,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
