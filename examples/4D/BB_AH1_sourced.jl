
using Jecco.AdS4_3_1_s

grid = SpecCartGrid3D(
    x_min            = -100000.0,
    x_max            =  100000.0,
    x_nodes          =  150,
    y_min            = -100000.0,
    y_max            =  100000.0,
    y_nodes          =  150,
    u_outer_min      =  0.3,
    u_outer_max      =  1.1,
    u_outer_domains  =  1,
    u_outer_nodes    =  12,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)

id = AdS4_3_1_s.BoostedBBnumerical(
    AH_pos = 1.0,
    IDdir = "/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/Sourced",
)

evoleq = AffineNull(
    #source = GaussianSource(time = 0.0, sigmax=0.8,sigmay=0.8, x0=0.0, y0=0.0,Amp=2.6,t0=-2, L=200000.0),
    source = NoSource(),
    gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    out_boundary_every  = 100,
    #out_bulk_every      = 10,
    remove_existing     = true,
    #out_gauge_every     = 1,
)

integration = Integration(
    #dt              = 0.002,
    tmax            = 30.0,
    #ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
