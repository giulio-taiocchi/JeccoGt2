
@inline S_inner_to_outer(S_in, u, xi,S0,S0_t) =  S0/u + S0*xi + S0_t + u*u*S_in#1/u + xi + u*u*S_in
    #

@inline Fx_inner_to_outer(F_in, u, S0, S0_x, S0_t, S0_tx) =  F_in*u + (S0*S0_tx - S0_t*S0_x)/(S0*S0)
@inline Fy_inner_to_outer(F_in, u, S0, S0_y, S0_t, S0_ty) =  F_in*u + (S0*S0_ty - S0_t*S0_y)/(S0*S0) 


@inline A_inner_to_outer(A_in, u, xi, S0, S0_t, S0_tt, S0_x, S0_xx, S0_y, S0_yy) = 1/(u*u)+ (2*xi)/u + u*A_in + (S0 * S0*S0_t * S0_t - 2*S0*S0*S0*S0_tt + S0_x*S0_x + S0_y*S0_y- S0*(S0_xx + S0_yy) + S0*S0*S0*S0*(xi*xi))/(S0*S0*S0*S0)#1/(u*u) + 2*xi/u + xi*xi + u * A_in #
	


@inline S_u_inner_to_outer(S_u_in, S_in, u, xi,S0) =
    -(S0/u^2) + 2*u*(S_in) + u^2*(S_u_in)
    #   -1/(u*u)+ 2 *u * S_in + u*u * S_u_in

@inline F_u_inner_to_outer(F_u_in, F_in, u) = F_in + u* F_u_in


@inline A_u_inner_to_outer(A_u_in, A_in, u, xi) =
    #-2/(u*u*u) - 2*xi/(u*u) + A_in + u* A_u_in
    -2/(u*u*u) - 2*xi/(u*u) + A_in + u* A_u_in


@inline Sd_inner_to_outer(Sd_in, u, xi, S0, S0_t, S0_x, S0_xx, S0_y, S0_yy) = S0/(2*u*u) + (S0*xi + (S0_t))/u + (S0*S0*S0*S0*xi*xi + 2*S0*S0*S0*xi*(S0_t) + S0*S0*S0_t*S0_t + S0_x*S0_x - S0*S0_xx + S0_y*S0_y - S0*S0_yy)/(2*S0*S0*S0) + u*Sd_in 
#1/(2*u*u) + xi/u + xi*xi/2  + u* Sd_in # sourceless	
    

@inline Bd_inner_to_outer(Bd_in, u) = u*u*Bd_in


