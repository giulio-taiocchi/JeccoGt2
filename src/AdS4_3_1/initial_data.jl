
Base.@kwdef struct BlackBrane_xi1{T} <: InitialData
    a30           :: T   = 1.0
    # guess for the AH position
    AH_pos        :: T   = 1.0
    xi_0          :: T   = 0.0
    xi_Ax         :: T   = 0.0
    xi_nx         :: Int = 1
    xi_Ay         :: T   = 0.0
    xi_ny         :: Int = 1
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    ahf           :: AHF = AHF()
end


abstract type ID_ConstantAH  <: InitialData end

Base.@kwdef struct BlackBrane{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    ahf           :: AHF = AHF()
end

Base.@kwdef struct BlackBranePert{T} <: ID_ConstantAH
    energy_dens   :: T   = 1.0
    AH_pos        :: T   = 1.0
    B_amp        :: T   = 0.0
    B_nx         :: Int = 1
    B_ny         :: Int = 2
    G_amp         :: T   = 0.0
    G_nx          :: Int = 1
    G_ny          :: Int = 2
    a3_ampx       :: T   = 0.0
    a3_ampy       :: T   = 0.0
    a3_kx         :: Int = 1
    a3_ky         :: Int = 1
    fx1_ampx      :: T   = 0.0
    fx1_ampy      :: T   = 0.0
    fx1_kx        :: Int = 1
    fx1_ky        :: Int = 1
    fy1_ampx      :: T   = 0.0
    fy1_ampy      :: T   = 0.0
    fy1_kx        :: Int = 1
    fy1_ky        :: Int = 1
    xi0           :: T   = 0.0
    xmax          :: T
    xmin          :: T
    ymax          :: T
    ymin          :: T
    ahf           :: AHF = AHF()
end


Base.@kwdef struct QNM_1D{T} <: InitialData
    energy_dens :: T   = 1.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
end

Base.@kwdef struct QNM_1DG{T} <: InitialData
    energy_dens :: T   = 1.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
end

Base.@kwdef struct vector_channel_QNM{T} <: InitialData
    energy_dens :: T   = 2.0
    AH_pos      :: T   = 1.0
    G_amp         :: T   = 0.0
    G_nx          :: Int = 1
    G_ny          :: Int = 1
    fx1_ampx    :: T   = 0.0
    fx1_ampy    :: T   = 0.0
    fx1_kx      :: Int = 1
    fx1_ky      :: Int = 1
    fy1_ampx    :: T   = 0.0
    fy1_ampy    :: T   = 0.0
    fy1_kx      :: Int = 1
    fy1_ky      :: Int = 1
    ahf         :: AHF = AHF()
end



Base.@kwdef struct BoostedBBSeminumerical{T} <: InitialData
    #energy_dens :: T   = 5.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
    a3_ampx	:: T = 0.0
    a3_translx	:: T = 0.0
    A		:: T = 0.0
    B		:: T = 0.0
    C		:: T = 0.0
    phase_a	:: T = 0.0
    phase_fx	:: T = 0.0
end

#does this should be InitialData ot ID_ConstantAh type

Base.@kwdef struct BoostedBBnumerical{T} <: InitialData
    #energy_dens :: T   = 5.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
    IDdir	:: AbstractString = "/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/"
end

Base.@kwdef struct BBnumerical{T} <: InitialData
    #energy_dens :: T   = 5.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
end

Base.@kwdef struct BB3Dnumerical{T} <: InitialData
    #energy_dens :: T   = 5.0
    AH_pos      :: T   = 1.0
    ahf         :: AHF = AHF()
end

function (id::InitialData)(bulkconstrains, bulkevols, bulkderivs, boundary::Boundary,
                           gauge::Gauge, horizoncache::HorizonCache, systems::SystemPartition,
                           evoleq::EvolutionEquations)
    _, Nx, Ny = size(systems[end])
    AH_pos    = id.AH_pos
    xi        = getxi(gauge)
    # function to solve the nested system
    nested = Nested(systems, bulkconstrains)

    init_data!(boundary, systems[1], id)

    init_data!(gauge, systems[end], id)
    init_data!(bulkevols, gauge, systems, id)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)
    # find the Apparent Horizon
    sigma = similar(gauge.xi)
    fill!(sigma, 1/AH_pos)  # initial guess
    find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
             horizoncache, systems[end], id.ahf)
    MaxAH=maximum(sigma)
    println("The horizon found is at rax=$MaxAH")
    println("AH_pos is $AH_pos")
    nothing
end



function (id::ID_ConstantAH)(bulkconstrains, bulkevols, bulkderivs, boundary::Boundary,
                           gauge::Gauge, horizoncache::HorizonCache, systems::SystemPartition,
                           evoleq::EvolutionEquations)
    _, Nx, Ny = size(systems[end])
    AH_pos    = id.AH_pos
    xi        = getxi(gauge)
    
    
    
    #printing u, added for numerical initial data
    #counting = 0
    #for nsys in systems
	    #println("domain $counting")
	    #actual_system = nsys
	    #Nu, Nx, Ny = size(actual_system)
	    #u_coordinates = actual_system.ucoord
	    #for a in 1:Nu
	    #	utoprint = u_coordinates[a]
	    #	println("u: $utoprint")
	    #end
	    #counting = counting +1
    #end
    
    

    # function to solve the nested system
    nested = Nested(systems, bulkconstrains)
    
   
    init_data!(boundary, systems[1], id)
    init_data!(gauge, systems[end], id)
    init_data!(bulkevols, gauge, systems, id)
    
    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)

    # find the Apparent Horizon
    sigma = similar(gauge.xi)
    fill!(sigma, 1/AH_pos)  # initial guess
    find_AH!(sigma, bulkconstrains[end], bulkevols[end], bulkderivs[end], gauge,
             horizoncache, systems[end], id.ahf)

    # assuming that the AH has been found, we now update xi and the bulk variables
    MaxAH=maximum(sigma)
    MaxOldxi= maximum(xi)
    println("The max horizon found is r=$MaxAH")
    println("old max xi is $MaxOldxi")
    AH_fixed_pos=evoleq.gaugecondition.u_AH
    for j in 1:Ny
        for i in 1:Nx
            xi[1,i,j] += -1 / AH_fixed_pos + sigma[1,i,j]
        end
    end

    init_data!(bulkevols, gauge, systems, id)

    # solve nested system for the constrained variables
    nested(bulkevols, boundary, gauge, evoleq)
    # TO DO printl(xi, and other)!!	
    # should AH_pos used up there be AH_pos fro the gauge_conditino structure!?
    # AH should now be at u = AH_pos

    nothing
end

function init_data!(bulkevols, gauge::Gauge, systems::SystemPartition,
                    id::InitialData)
    # the Ref() makes its argument a scalar with respect to broadcast
    #init_data!.(bulkevols, Ref(gauge), systems, Ref(id))
    counting = 0
    for k in systems
    	init_data!.(Ref(bulkevols.x[counting+1]), Ref(gauge), Ref(k), Ref(id),Ref(counting))
    	counting = counting + 1
    end
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Inner},
                    id::InitialData,counting)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord
    B  = getB(bulk)
    G   = getG(bulk)
    xi  = getxi(gauge)
    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u       = uu[a]
                x       = xx[i]
                y       = yy[j]
                xi_ij   = xi[1,i,j]
                aux     = 1 + xi_ij * u
                aux3    = aux * aux * aux
                u_old   = u / aux
                B_old   = analytic_B(a, i, j, u_old, x, y, id, counting)
                G_old   = analytic_G(a, i, j, u_old, x, y, id, counting)
                #B_old  = analytic_B(u_old, x, y, id)
                #G_old   = analytic_G(u_old, x, y, id)
                B[a,i,j]  = B_old / aux3
                G[a,i,j]  = G_old / aux3
                #Baij = B[a,i,j]
                #if j==5
	#		if i==5
	#			println("outside analytic  u= $u, aux= $aux, xi= $xi_ij, Binner= $Baij")
	#		end
	#	end
                
            end
        end
    end


    bulk
end

function init_data!(bulk::BulkEvolved, gauge::Gauge, sys::System{Outer},
                    id::InitialData,counting)
    Nu, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    uu = sys.ucoord

    B  = getB(bulk)
    G   = getG(bulk)
    xi  = getxi(gauge)

    for j in 1:Ny
        for i in 1:Nx
            for a in 1:Nu
                u         = uu[a]
                x         = xx[i]
                y         = yy[j]
                xi_ij     = xi[1,i,j]
                aux       = 1 + xi_ij * u
                aux3      = aux * aux * aux
                aux4      = aux * aux3
                u_old     = u / aux
                
                B_old     = analytic_B(a, i, j, u_old, x, y, id, counting)
                
                G_old     = analytic_G(a, i, j, u_old, x, y, id, counting)
                B_inner   = B_old / aux3
                G_inner   = G_old  / aux3
		#if j==5
		#	if i==5
		#		println("outside analytic  u= $u, aux= $aux, xi= $xi_ij, Binner= $B_inner, uold= $u_old")
		#	end
		#end
                B[a,i,j]  = u^3 * B_inner
                G[a,i,j]  = u^3 * G_inner
                #Baij = B[a,i,j]
                #if j==5
		#	if i==5
		#	println("post B is $Baij")
		#	end
		#end
            end
        end
    end

    

    bulk
end


# BlackBrane_xi1

analytic_B(i, j, k, u, x, y, id::BlackBrane_xi1, whichsystem)   = 0
analytic_G(i, j, k, u, x, y,  id::BlackBrane_xi1, whichsystem)    = 0

function init_data!(ff::Boundary, sys::System, id::BlackBrane_xi1)
    a30 = id.a30

    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBrane_xi1)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    xi  = getxi(ff)

    xmax  = id.xmax
    xmin  = id.xmin
    ymax  = id.ymax
    ymin  = id.ymin

    for j in 1:Ny
        for i in 1:Nx
            x         = xx[i]
            y         = yy[j]

            xi[1,i,j] = id.xi_0 +
                id.xi_Ax * sin( 2 * π * id.xi_nx * (xmax-x)/(xmax-xmin) ) +
                id.xi_Ay * sin( 2 * π * id.xi_ny * (ymax-y)/(ymax-ymin) )
        end
    end

    ff
end


# BlackBrane initial data

analytic_B(i, j, k, u, x, y, id::BlackBrane, whichsystem)  = 0
analytic_G(i, j, k, u, x, y, id::BlackBrane, whichsystem)   = 0

function init_data!(ff::Boundary, sys::System, id::BlackBrane)
    #a30 = -id.energy_dens/2
    a30 = -1
    

    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBrane)
    #a30     = -id.energy_dens/2
    a30 = -1
    AH_pos  = id.AH_pos
    #xi0     = (-a30)^0.25 - 1/AH_pos
    xi0 = 0

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end


# BlackBranePert initial data


function analytic_B(i, j, k,u, x, y, id::BlackBranePert, whichsystem)
    # add the perturbation on B
    pert_amp = id.B_amp
    xmax     = id.xmax
    xmin     = id.xmin
    ymax     = id.ymax
    ymin     = id.ymin
    # number of maxima in each direction
    nx       = id.B_nx
    ny       = id.B_ny

    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
        sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end


function analytic_G(i, j, k,u, x, y, id::BlackBranePert, whichsystem)
    # add the perturbation on G
    pert_amp = id.G_amp
    xmax     = id.xmax
    xmin     = id.xmin
    ymax     = id.ymax
    ymin     = id.ymin
    # number of maxima in each direction
    nx       = id.G_nx
    ny       = id.G_ny

    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
        sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end

function init_data!(ff::Boundary, sys::System{Inner}, id::BlackBranePert)
    epsilon = id.energy_dens

    # a3 perturbation amplitude
    ampx     = id.a3_ampx
    ampy     = id.a3_ampy
    fx1_ampx = id.fx1_ampx
    fx1_ampy = id.fx1_ampy
    fy1_ampx = id.fy1_ampx
    fy1_ampy = id.fy1_ampy
    # number of maxima
    kx     = id.a3_kx
    ky     = id.a3_ky
    fx1_kx = id.fx1_kx
    fx1_ky = id.fx1_ky
    fy1_kx = id.fy1_kx
    fy1_ky = id.fy1_ky
    
    xmax = id.xmax
    xmin = id.xmin
    xmid = (xmax + xmin) / 2
    ymax = id.ymax
    ymin = id.ymin
    ymid = (ymax + ymin) / 2

    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord

    a30   = (-epsilon) / 2

    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    for j in 1:Ny
        for i in 1:Nx
            x = xx[i]
            y = yy[j]
            a3[1,i,j]  += -a30 * ( ampx * cos(2 * π * kx * (x-xmid)/(xmax-xmin)) +
                                   ampy * cos(2 * π * ky * (y-ymid)/(ymax-ymin)) )

            fx1[1,i,j] += fx1_ampx * cos(2 * π * fx1_kx * (x-xmid)/(xmax-xmin)) +
                           fx1_ampy * cos(2 * π * fx1_ky * (y-ymid)/(ymax-ymin))

            fy1[1,i,j] += fy1_ampx * cos(2 * π * fy1_kx * (x-xmid)/(xmax-xmin)) +
                          fy1_ampy * cos(2 * π * fy1_ky * (y-ymid)/(ymax-ymin))
        end
    end

    ff
end

function init_data!(ff::Gauge, sys::System, id::BlackBranePert)
    a30     = -id.energy_dens/2
    AH_pos  = id.AH_pos

    # TODO: this guess works best for the conformal case. is there a better one?
    if id.xi0 == 0
        xi0 = (-a30)^0.25 - 1/AH_pos
    else
        xi0 = id.xi0
    end

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end




#QNM in 1D initial data
analytic_B(i, j, k,u, x, y, id::QNM_1D, whichsystem)  =  3/2*0.1 * u^6
analytic_G(i, j, k,u, x, y, id::QNM_1D, whichsystem)  = 0

function init_data!(ff::Boundary, sys::System, id::QNM_1D)
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    epsilon = id.energy_dens

    a30 = (-epsilon) / 2

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::QNM_1D)
    epsilon = id.energy_dens
    AH_pos  = id.AH_pos

    a30 = (-epsilon) / 2

    xi0 = 0

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

#QNM in 1D initial data for the tensor sector
analytic_B(i, j, k,u, x, y, id::QNM_1DG, whichsystem)  = 0
analytic_G(i, j, k,u, x, y, id::QNM_1DG, whichsystem)  = 3/2*0.1 * u^6

function init_data!(ff::Boundary, sys::System, id::QNM_1DG)
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    epsilon = id.energy_dens

    a30 = (-epsilon) / 2

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::QNM_1DG)
    epsilon = id.energy_dens
    AH_pos  = id.AH_pos

    a30 = (-epsilon) / 2

    #xi0 = 0
    
    xi0 = 0
    
    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

#QNM in 1D initial data for the vector sector
analytic_B(i, j, k,u, x, y, id::vector_channel_QNM, whichsystem)  = 0
function analytic_G(i, j, k,u, x, y, id::vector_channel_QNM, whichsystem)
	# add the perturbation on G
	    pert_amp = id.G_amp
	    # number of maxima in each direction
	    nx       = id.G_nx
	    ny       = id.G_ny

	    pert_amp * sin( 2 * π * nx * (xmax-x)/(xmax-xmin) ) *
		sin( -2 * π * ny * (ymax-y)/(ymax-ymin) )
end

function init_data!(ff::Boundary, sys::System, id::vector_channel_QNM)
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    epsilon = id.energy_dens

    a30 = (-epsilon) / 2

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::vector_channel_QNM)
    epsilon = id.energy_dens
    AH_pos  = id.AH_pos

    a30 = (-epsilon) / 2

    #xi0 = 0
    
    xi0 = 0
    
    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

#numerical Black Brane
function analytic_B(i, j, k, u, x, y, id::BBnumerical, whichsystem)
	uu = u
	uprec = precision(uu)
	initialB=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/InitialB_BB.h5")
	system_index = string(whichsystem+1)
	dset=initialB[system_index]
	B=read(dset)
	Bvalue = parse(Float64,B[i])
	Bprec = precision(Bvalue)
	if j==5
		if k==5
			@printf("%.40f\n", Bvalue)
			
			println("u=$uu , index: $i, B = $Bvalue")
		end
	end
	#println("B in u=$uu index: $i is $Bvalue")
	#println("THIS IS SYSTEM NUMBER $whichsystem")
	Bvalue
end
analytic_G(i, j, k, u, x, y, id::BBnumerical,whichsystem)  = 0

function init_data!(ff::Boundary, sys::System, id::BBnumerical)
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    #epsilon = id.energy_dens

    a30 = -1
    #a30 = -2

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!( ff::Gauge, sys::System, id::BBnumerical)
    #epsilon = id.energy_dens
    AH_pos  = id.AH_pos

    a30 = (-1)
    #a30 = -2

    #xi0 = 0.19931437035694333
    xi0 = 0

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

#numerical Black Brane 3D numerical initial data. Used for checking the 3D input system
function analytic_B(i, j, k, u, x, y, id::BB3Dnumerical, whichsystem)
	uu = u
	uprec = precision(uu)
	initialB=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/InitialB_BB.h5")
	system_index = string(whichsystem+1)
	dset=initialB[system_index]
	B=read(dset)
	Bvalue=B[k,j,i]
	#println("BID= $bkji , xGrid= $x")
	#println(" ")
	Bvalue
end
analytic_G(i, j, k, u, x, y, id::BB3Dnumerical,whichsystem)  = 0

function init_data!(ff::Boundary, sys::System, id::BB3Dnumerical)
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)

    #epsilon = id.energy_dens

    a30 = -1
    #a30 = -2

    fill!(a3, a30)
    fill!(fx1, 0)
    fill!(fy1, 0)

    ff
end

function init_data!(ff::Gauge, sys::System, id::BB3Dnumerical)
    #epsilon = id.energy_dens
    AH_pos  = id.AH_pos

    a30 = -1
    #a30 = -2

    #xi0 = 0.19931437035694333
    xi0 = 0

    xi  = getxi(ff)

    fill!(xi, xi0)

    ff
end

#Semi numerical boosted Black Brane
function analytic_B(i, j, k, u, x, y, id::BoostedBBSeminumerical, whichsystem)
	uu = u
	initialB=h5open("/home/giulio/University/PhD/JeccoNewTest/Jecco_G/examples/InitialB_BBB.h5")
	system_index = string(whichsystem+1)
	dset=initialB[system_index]
	B=read(dset)
	# here the indecex have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	Bvalue = B[k,j,i]
	
	
	Bvalue
end
analytic_G(i, j, k, u, x, y, id::BoostedBBSeminumerical,whichsystem)  = 0



function init_data!(ff::Boundary, sys::System, id::BoostedBBSeminumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)
    ampx = id.a3_ampx
    transl = id.a3_translx
    phfx = id.phase_fx
    pha = id.phase_a
    AA = id.A
    BB = id.B
    CC = id.C


    fill!(a3, 0)
    fill!(fx1, 0)
    fill!(fy1, 0)
    for j in 1:Ny
        for i in 1:Nx      
                x = xx[i]
                y = yy[j]         
                a3[1,i,j] = transl + ampx * cos(pha*π*y)
                fx1[1,i,j] = AA * cos(phfx*π*y) * sqrt(BB+ CC * cos(phfx*π*y) *cos(phfx*π*y))
        end
    end
    ff
end

function init_data!(ff::Gauge, sys::System, id::BoostedBBSeminumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    AH_pos  = id.AH_pos
    
    xi  = getxi(ff)
    fill!(xi, 0)
    
    
    #for j in 1:Ny
    #    for i in 1:Nx      
    #            x = xx[i]
    #            y = yy[j]         
    #            xi[1,i,j] = 1/4*(1-cos(2*π*x)*cos(2*π*x))
    #    end
    #end
    ff
end


# Numerical boosted Black Brane
function analytic_B(i, j, k, u, x, y, id::BoostedBBnumerical, whichsystem)
	uu = u
	dir = id.IDdir
	Bdirectory = dir*"InitialB_BBB.h5"
	initialB=h5open(Bdirectory)
	system_index = string(whichsystem+1)
	dset=initialB[system_index]
	B=read(dset)
	# here the indexes have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	Bvalue = B[k,j,i]
	
	
	Bvalue
end

function analytic_G(i, j, k, u, x, y, id::BoostedBBnumerical,whichsystem) 
	uu = u
	dir = id.IDdir
	Gdirectory = dir*"InitialG_BBB.h5"
	initialG=h5open(Gdirectory)
	system_index = string(whichsystem+1)
	dset=initialG[system_index]
	G=read(dset)
	# here the indexes have to be inverted since julia and mathematica input and output mechanism is the opposite
	# should be B[i,j,k]
	Gvalue = G[k,j,i]
		
	Gvalue
end



function init_data!(ff::Boundary, sys::System, id::BoostedBBnumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    
    a3  = geta3(ff)
    fx1 = getfx1(ff)
    fy1 = getfy1(ff)
    


    fill!(a3, 0)
    fill!(fx1, 0)
    fill!(fy1, 0)
    dir = id.IDdir
    a3directory = dir*"Initiala3_BBB.h5"
    fx1directory = dir*"Initialfx_BBB.h5"
    fy1directory = dir*"Initialfy_BBB.h5"
    a3data = h5open(a3directory)
    fxdata = h5open(fx1directory)
    fydata = h5open(fy1directory)
    a3read = read(a3data["a3"])
    fxread = read(fxdata["fx"])
    fyread = read(fydata["fy"])
    for j in 1:Ny
        for i in 1:Nx      
                x = xx[i]
                y = yy[j]         
                a3[1,i,j] = a3read[j,i]
                fx1[1,i,j] = fxread[j,i]
                fy1[1,i,j] = fyread[j,i]
        end
    end
    ff
end

function init_data!(ff::Gauge, sys::System, id::BoostedBBnumerical)
    _, Nx, Ny = size(sys)
    xx = sys.xcoord
    yy = sys.ycoord
    AH_pos  = id.AH_pos
    
    
    xi  = getxi(ff)
    fill!(xi, 0)
    dir = id.IDdir
    xidirectory = dir*"Initialxi_BBB.h5"
    xidata = h5open(xidirectory)
    xiread = read(xidata["xi"])
    for j in 1:Ny
        for i in 1:Nx      
                x = xx[i]
                y = yy[j]         
                xi[1,i,j] = xiread[j,i]                
        end
    end
    #for j in 1:Ny
    #    for i in 1:Nx      
    #            x = xx[i]
    #            y = yy[j]         
    #            xi[1,i,j] = 1/4*(1-cos(2*π*x)*cos(2*π*x))
    #    end
    #end
    ff
end



