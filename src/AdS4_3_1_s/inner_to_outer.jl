
@inline S_inner_to_outer(S_in, u, xi,S0,S0_t) =
    #S0/u + S_in*u^2 + S0*xi + (S0_t)
    S0/u + S0*xi + (S0_t) + u^2*S_in

@inline Fx_inner_to_outer(F_in, u, S0, S0_x, S0_t, S0_tx) =  F_in*u + (S0*S0_tx - S0_t*S0_x)/(S0*S0)
@inline Fy_inner_to_outer(F_in, u, S0, S0_y, S0_t, S0_ty) =  F_in*u + (S0*S0_tx - S0_t*S0_x)/(S0*S0) 


@inline A_inner_to_outer(A_in, u, xi, S0, S0_t, S0_tt, S0_x, S0_xx, S0_y, S0_yy) = 1/(u*u) + A_in*u + (2*xi)/u +(S0*S0*S0_t*S0_t - 2*S0*S0*S0*S0_tt + S0_x*S0_x + S0_y*S0_y - S0*(S0_xx + S0_yy) + S0*S0*S0*S0*(xi*xi - 2*xi_t))/(S0*S0*S0*S0)
	


@inline S_u_inner_to_outer(S_u_in, S_in, u, xi,S0) =-1/(u*u)+ 2 *u * S_in + u*u * S_u_in
   # -(S0/u^2) + 2*u*(S_in) + u^2*(S_u_in)
       

@inline F_u_inner_to_outer(F_u_in, F_in, u) = F_in + u* F_u_in


@inline A_u_inner_to_outer(A_u_in, A_in, u, xi) =
    #-2/(u*u*u) - 2*xi/(u*u) + A_in + u* A_u_in
    -2/(u*u*u) - 2*xi/(u*u) + A_in + u* A_u_in


@inline Sd_inner_to_outer(Sd_in, u, xi, S0, S0_t, S0_x, S0_xx, S0_y, S0_yy) =
	####S0/(2*u^2) + u*Sd_in + (S0*xi +S0_t)/u +(S0*xi^2/2 + 2*xi*S0_t/2 + S0_t^2/(2*S0) + S0_x^2/(2*S0^3) - S0_xx/(2*S0^2) + S0_y^2/(2*S0^3) - S0_yy/(2*S0^2))###
	#(4*S0^4*xi^2 + 8*S0^3*xi*(S0_t) + 4*S0^2*(S0_t)^2 + 4*(S0_x)^2 - 4*S0*(S0_xx) + 4*(S0_y)^2 - 4*S0*(S0_yy))/(8*S0^3)
    S0/(2*u*u) + (S0*xi + (S0_t))/u + (S0*S0*S0*S0*xi*xi + 2*S0*S0*S0*xi*(S0_t) + S0*S0*S0_t*S0_t + S0_x*S0_x - S0*S0_xx + S0_y*S0_y - S0*S0_yy)/(2*S0*S0*S0) + u*Sd_in 

@inline Bd_inner_to_outer(Bd_in, u) = u*u*Bd_in


