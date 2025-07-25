
function setup_rhs(tmp::EvolVars, bulkconstrains::BulkPartition{Nsys},
                   bulkderivs::BulkPartition{Nsys},
                   cache::HorizonCache, systems::SystemPartition,
                   integration::Integration,tzero) where {Nsys}
    # function to solve the nested system
    nested = Nested(systems, bulkconstrains, bulkderivs)

    diss_cache = similar(tmp)

    bulkevols_cache   = getbulkevolvedpartition(diss_cache)
    boundary_cache    = getboundary(diss_cache)
    gauge_cache       = getgauge(diss_cache)

    function (ff_t::EvolVars, ff::EvolVars, evoleq::EvolutionEquations, t)
        #println("t in setup_RHS is $t")
	evoleq.source.time = t+tzero
	evoleq.source.step += 1
	timing = getSourcetime(evoleq.source)
	#println("source time in setup_RHS is $timing")
        bulkevols_t = getbulkevolvedpartition(ff_t)
        boundary_t  = getboundary(ff_t)
        gauge_t     = getgauge(ff_t)

        bulkevols   = getbulkevolvedpartition(ff)
        boundary    = getboundary(ff)
        gauge       = getgauge(ff)


        # filter state vector after each integration (sub)step.
        vprint("INFO: filtering")

        # Note: one problem with this approach is that it always filters the
        # state vector upon calling, ie, if this function is called twice with
        # the *same* state vector (say, for error estimate purposes), the state
        # vector will be filtered *twice*.
        if t > 0 && integration.filter_poststep
            @inbounds for aa in 1:Nsys
                sys            = systems[aa]
                bulkevol       = bulkevols[aa]
                bulkevol_cache = bulkevols_cache[aa]

                apply_dissipation!(bulkevol, bulkevol_cache, sys)

                # exponential filter
                sys.filters(bulkevols[aa])
            end
            apply_dissipation!(boundary, boundary_cache, systems[1])
            apply_dissipation!(gauge, gauge_cache,  systems[Nsys])
        end


        vprint("INFO: compute_boundary_t")
        compute_boundary_t!(boundary_t, bulkevols[1], boundary, gauge, systems[1], evoleq)

        # solve nested system for the constrained variables
        vprint("INFO: nested system")
        #println("enter nested from setup rhs")
        nested(bulkevols, boundary, gauge, evoleq)
	#println("exit nested from setup rhs")
        vprint("INFO: compute_xi_t")
        compute_xi_t!(gauge_t, bulkconstrains[Nsys], bulkevols[Nsys], bulkderivs[Nsys],
                      gauge, cache, systems[Nsys], evoleq.gaugecondition,evoleq)

        vprint("INFO: bulkevolved_t")
        # TODO: check if this loop is thread-safe
        # @inbounds @threads for aa in 1:Nsys
        @inbounds for aa in 1:Nsys
            sys           = systems[aa]
            bulkevol_t    = bulkevols_t[aa]
            bulkevol      = bulkevols[aa]
            bulkconstrain = bulkconstrains[aa]

            compute_bulkevolved_t!(bulkevol_t, bulkconstrain, gauge_t, bulkevol,
                                   boundary, gauge, sys, evoleq)
        end
        sync_bulkevolved!(bulkevols_t, bulkconstrains, gauge_t, systems, evoleq)

	
	
        nothing
    end
end
