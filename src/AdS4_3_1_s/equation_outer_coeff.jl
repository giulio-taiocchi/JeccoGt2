
#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_outer("B1")

should expand to

  B1t = B1_x -  (Fx + xi_x) * B1p

etc.

=#
macro tilde_outer(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx + xi_x) * $fp) )
end
macro hat_outer(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy + xi_y)  * $fp) )
end
macro bar_outer(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx + xi_x)  * ( -2*($fp_x) + (Fx + xi_x) * ($fpp) )) )
end
macro star_outer(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy + xi_y)  * ( -2*($fp_y) + (Fy + xi_y) * ($fpp) )) )
end
macro cross_outer(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx + xi_x) * ($fp_y) -
                 (Fy + xi_y)  * ( $fp_x - (Fx + xi_x)  * ($fpp) ) ) )
end


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (S0, S0_t, u, xi, B, Bp, G, Gp) = vars
	
	ABCS[1] = 4*u^4

	ABCS[2] = 8*u^3

	ABCS[3] = Gp^2 + Bp^2*cosh(G)^2

	ABCS[4] =0


    nothing
end

# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
       S0, S0_x, S0_y, S0_t, S0_tx, S0_ty, u, xi, xi_x, xi_y,
        B     ,        G      ,        S      ,
        Bp    ,        Gp     ,        Sp     ,
        Bpp   ,        Gpp    ,        Spp    ,
        B_x   ,        G_x    ,        S_x    ,
        B_y   ,        G_y    ,        S_y    ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,
        Bp_y  ,        Gp_y   ,        Sp_y
    ) = vars

 



x0 = S ^ 2
x1 = 2 * x0
x2 = u ^ 4 * x1
x3 = u ^ 2
x4 = 2 * G
x5 = cosh(x4)
x6 = Bp * x5
x7 = sinh(G)
x8 = cosh(G)
x9 = Bp * x7 * x8
x10 = exp(B)
x11 = x1 * x10
x12 = x11 * x3
x13 = Bp ^ 2
x14 = 2 * Bp
x15 = Bpp * S + Sp * x14
x16 = S * x13 - x15
x17 = x8 ^ 2
x18 = 2 * S
x19 = x17 * x18
x20 = Gp * x1
x21 = sinh(x4)
x22 = Bp * x21
x23 = x20 * x22
x24 = Sp ^ 2
x25 = Gp ^ 2 * x0
x26 = 4 * S
x27 = Spp * x26 - 4 * x24 + 2 * x25
x28 = S * x13 + x15
x29 = Gp * S * x14
x30 = Gp * x18 * x6
x31 = 4 * Sp
x32 = Gp * x31 + Gpp * x18
x33 = S * x10 * (x21 * x28 - x29 + x30 - x32)
x34 = x10 * x26
x35 = S_y * x34
x36 = B_y * x10
x37 = Bp * x17
x38 = x1 * x17
x39 = x0 * x21
x40 = Bp_y * x10
x41 = Bp * x39
x42 = G_y * x11
x43 = S * x17
x44 = Gp * x41
x45 = Spp * x18 - 2 * x24 + x25
x46 = x16 * x21 + x29 - x30 - x32
AA[1,1] = x2
AA[1,2] = 0
BB[1,1] = x0 * x3 * (Bp + 4 * u + x6)
BB[1,2] = x12 * (Gp - x9)
CC[1,1] = x16 * x19 - x23 + x27
CC[1,2] = x33
SS[1] = -B_x * x1 * x37 + Bp_x * x38 + G_x * x1 * x22 - G_x * x20 + Gp * x35 + Gp_y * x11 + S_x * x26 * x37 + S_x * x31 - Sp_x * x26 + x20 * x36 + x33 * xi_y - x35 * x9 - x36 * x41 - x39 * x40 - x42 * x6 + 2 * xi_x * (x16 * x43 - x44 + x45)
AA[2,1] = 0
AA[2,2] = x10 * x2
BB[2,1] = x1 * x3 * (Gp + x9)
BB[2,2] = x12 * (2 * u - x37)
CC[2,1] = S * x46
CC[2,2] = x10 * (x19 * x28 + x23 + x27)
SS[2] = -B_x * x20 - B_x * x41 - B_y * x11 * x37 + 2 * Bp * G_x * x0 * x5 + 4 * Bp * S * S_x * x7 * x8 + Bp_x * x0 * x21 - G_y * x10 * x20 + 4 * Gp * S * S_x + 2 * Gp_x * x0 + S * x46 * xi_x + 4 * S_y * Sp * x10 - Sp_y * x34 + 2 * x10 * xi_y * (x28 * x43 + x44 + x45) - x22 * x42 - x35 * x37 - x38 * x40
    
    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
         S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy,u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,        G      ,        S      ,    Fx     ,    Fy     ,
        Bp    ,        Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,        Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,        G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,        G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,        Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,        G_xx   ,        S_xx   ,
        B_yy  ,        G_yy   ,        S_yy   ,
                        G_xy   ,        S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Bp")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Bp")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")

    
   
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    x0 = exp(B)
x1 = 8 * x0
x2 = S ^ 2
x3 = cosh(G)
x4 = 4 * x3
x5 = sinh(G)
x6 = 2 * Gt
x7 = 2 * x0 * x2
x8 = 4 * St
x9 = Fxt + xi_xx
x10 = 4 * Sp
x11 = Fx + xi_x
x12 = 4 * x11
x13 = x11 ^ 2
x14 = 4 * Spp
x15 = Spp * x11
x16 = Spt + x15
x17 = 4 * Sh
x18 = Fy + xi_y
x19 = Fxh + Fyt + 2 * xi_xy
x20 = 2 * x18
x21 = -Fyp
x22 = Gpp * x11
x23 = 2 * Fy + 2 * xi_y
x24 = Gpt + x22
x25 = 2 * x9
x26 = 2 * x11
x27 = 2 * x13
x28 = 2 * Bt
x29 = 2 * Fx + 2 * xi_x
x30 = Fyh + xi_yy
x31 = 4 * x18
x32 = x18 ^ 2
x33 = 2 * x30
x34 = 2 * x32
x35 = 2 * Bh
x36 = 2 * Bph
ABCS[1] = 0
ABCS[2] = -S ^ 3 * u ^ 2 * x1
ABCS[3] = Sp * x1 * x2
ABCS[4] = -12 * S ^ 4 * x0 + S * x0 * (x3 * (-Gh * x8 - Gt * x17) + x5 * (-8 * Sc + 4 * Sp * x19 + 8 * Spt * x18 - 8 * x15 * x18 - x16 * (8 * Fy + 8 * xi_y))) + S * (Gt * x5 * x8 - x3 * (Bt * x8 - 4 * Sb + Spt * x12 + x10 * x9 - x12 * x16 + x13 * x14)) + Sh * St * x1 * x5 - St ^ 2 * x4 + x2 * (x3 * (-2 * Bb + Bp * x25 + Bpp * x27 + Bpt * x26 + 2 * Bt ^ 2 + Fxp ^ 2 + Fxp * x28 - 2 * Fxpt + 2 * Gt ^ 2 - x29 * (Bpp * x11 + Bpt)) + x5 * (2 * Gb - Gp * x25 - Gpp * x27 - Gpt * x26 + x24 * x29 - x6 * (Fxp + x28))) + x3 * x7 * (-2 * Gc + Gh * (Bt + Fxp) + Gp * x19 + Gpt * x20 + Gt * (-Bh - x21) - x20 * x22 - x23 * x24) + x5 * x7 * (-Fxp * Fyp + Fxph + Fypt - Gh * x6) + (S * (Gh * x17 * x5 + x3 * (Bh * x17 - Sph * x31 + 4 * Ss - x10 * x30 + x14 * x32 + x31 * (Sph + Spp * x18))) - Sh ^ 2 * x4 + x2 * (x3 * (2 * Bh ^ 2 - Bp * x33 + Bpp * x34 + 2 * Bs + Fyp ^ 2 - Fyp * x35 - 2 * Fyph + 2 * Gh ^ 2 - x18 * x36 + x18 * (Bpp * x20 + x36)) + x5 * (2 * Gh * (x21 + x35) - Gp * x33 - Gph * x20 + Gpp * x34 + 2 * Gs + x23 * (Gph + Gpp * x18)))) * exp(2 * B)

    nothing
end


# this is another coupled equation, for Bd and Gd. the notation used is
#
# ( A11 d_uu Bd + A12 d_uu Gd + B11 d_u Bd + B12 d_u Gd + C11 Bd + C12 Gd ) = -S1
# ( A21 d_uu Bd + A22 d_uu Gd + B21 d_u Bd + B22 d_u Gd + C21 Bd + C22 Gd ) = -S2

function BdGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
       S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy,u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,      G      ,        S      ,    Fx     ,    Fy     ,  Sd,
        Bp    ,      Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,      Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,      G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,      G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,      Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,      Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,      G_xx   ,        S_xx   ,
        B_yy  ,      G_yy   ,        S_yy   ,
                      G_xy   ,        S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")


    
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    tanhG   = tanh(G)
    sinhG   = sinh(G)
    sechG   = sech(G)




x0 = exp(B)
x1 = 8 * x0
x2 = S ^ 3
x3 = u ^ 2 * x2
x4 = tanh(G)
x5 = S ^ 2
x6 = x1 * x5
x7 = Bp * x2
x8 = 4 * S
x9 = Fxp ^ 2
x10 = exp(2 * B)
x11 = Fyp * Sh
x12 = Fyp ^ 2 * S
x13 = 4 * x0
x14 = sinh(G)
x15 = cosh(G)
x16 = 2 * S
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x1 * x3
BB[1,2] = 0
CC[1,1] = x6 * (Gp * S * x4 + Sp)
CC[1,2] = x1 * x4 * x7
SS[1] = Bp * Sd * x6 + (8 * Fxp * St - Fxpt * x8 + 2 * S * x9 - x0 * x8 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp) - x10 * (-Fyph * x8 + 8 * x11 + 2 * x12)) * sech(G)
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x13 * x3
CC[2,1] = -2 * x0 * x7 * sinh(2 * G)
CC[2,2] = Sp * x13 * x5
SS[2] = -4 * Fxp * St * x14 + 4 * Gp * Sd * x0 * x5 - S * x14 * (-2 * Fxpt + x9) - x0 * x15 * x16 * (Bp * (Fxh - Fyt) + Bt * Fyp - Fxp * (Bh + Fyp) + Fxph + Fypt) + x0 * x15 * (4 * Fxp * Sh + 4 * Fyp * St) - x10 * x14 * (-Fyph * x16 + 4 * x11 + x12)
    nothing
end


function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
         S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd, Bd, Gd, 
        Bp  ,  Gp  ,  Sp   , Fxp   , Fyp   ,
        Bpp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,
        B_x ,  G_x ,  S_x  , Fx_x  , Fy_x  ,
        B_y ,  G_y ,  S_y  , Fx_y  , Fy_y  ,
        Bp_x,  Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        Bp_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,
        B_xx,  G_xx,  S_xx ,
        B_yy,  G_yy,  S_yy ,
                G_xy,  S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")


    #expB1   = exp(B1)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


x0 = exp(B)
x1 = S ^ 4
x2 = 2 * x1
x3 = cosh(G)
x4 = 4 * x3
x5 = sinh(G)
x6 = St * x5
x7 = 4 * Bt
x8 = Fxt + xi_xx
x9 = 4 * Sp
x10 = Fx + xi_x
x11 = 4 * x10
x12 = x10 ^ 2
x13 = 4 * Spp
x14 = Spp * x10
x15 = Spt + x14
x16 = S ^ 2
x17 = 2 * x8
x18 = 2 * x10
x19 = 2 * x12
x20 = 2 * Fx + 2 * xi_x
x21 = Gpp * x10
x22 = Gpt + x21
x23 = Fy + xi_y
x24 = 2 * x23
x25 = Fxh + Fyt + 2 * xi_xy
x26 = 2 * Fy + 2 * xi_y
x27 = 4 * x23
x28 = 4 * Bh
x29 = Fyh + xi_yy
x30 = x23 ^ 2
x31 = 2 * x30
ABCS[1] = u ^ 4 * x0 * x2
ABCS[2] = 4 * u ^ 3 * x0 * x1
ABCS[3] = 0
ABCS[4] = S * (-4 * Gt * x6 + x3 * (-4 * Sb + Spt * x11 + St * x7 - x11 * x15 + x12 * x13 + x8 * x9)) + St ^ 2 * x4 + x0 * (4 * S * (x3 * (Gh * St + Gt * Sh) + x5 * (2 * Sc - Sp * x25 - Spt * x24 + x14 * x24 + x15 * x26)) - 8 * Sh * x6 + x16 * (-8 * Sd * Sp + x3 * (2 * Bh * Gt - 2 * Bt * Gh + 4 * Gc - 2 * Gp * x25 - Gpt * x27 + x21 * x27 + 2 * x22 * x26) + x5 * (-2 * Fxp * Fyp + 4 * Gh * Gt)) + x2 * (Bd * Bp * x3 ^ 2 + Gd * Gp)) + x16 * (x3 * (2 * Bb - Bp * x17 - Bpp * x19 - Bpt * x18 - 2 * Bt ^ 2 + Fxp ^ 2 - 2 * Gt ^ 2 + x20 * (Bpp * x10 + Bpt)) + x5 * (-2 * Gb + Gp * x17 + Gpp * x19 + Gpt * x18 + Gt * x7 - x20 * x22)) + (S * (-4 * Gh * Sh * x5 - x3 * (Sh * x28 - Sph * x27 + 4 * Ss + x13 * x30 + x27 * (Sph + Spp * x23) - x29 * x9)) + Sh ^ 2 * x4 + x16 * (x3 * (-2 * Bh ^ 2 + 2 * Bp * x29 + 2 * Bph * x23 - Bpp * x31 - 2 * Bs + Fyp ^ 2 - 2 * Gh ^ 2 - x23 * (2 * Bph + Bpp * x24)) + x5 * (-Gh * x28 + 2 * Gp * x29 + 2 * Gph * x23 - Gpp * x31 - 2 * Gs - x26 * (Gph + Gpp * x23)))) * exp(2 * B)
    nothing
end

function xi_t_eq_coeff(vars::Tuple, ::Outer)
    (
       S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,  kappa, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd ,  Bd  , Gd,  A   ,
        Bp  ,  Gp  ,  Sp   , Fxp   , Fyp   , Sdp,  Bdp , Gdp, Ap  ,
        Bpp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,                   App ,
        B_x ,   G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x, Bd_x, Gd_x,A_x ,
	B_y ,  G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y, Bd_y,Gd_y, A_y ,
        Bp_x,  Gp_x,  Sp_x , Fxp_x , Fyp_x ,                   Ap_x,
        Bp_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,                   Ap_y,
                                 Fy_xx ,       A_xx,
                                          Fx_yy ,               A_yy,
                                          Fx_xy , Fy_xy ,       A_xy,Fx_xx,Fy_yy,
                                          B_xx,B_yy,B_xy,G_xx,G_yy,G_xy,S_xx,S_yy,S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")
    @tilde_outer("Bd")
    @tilde_outer("Gd")
    @tilde_outer("A")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")
    @hat_outer("Bd")
    @hat_outer("Gd")
    @hat_outer("A")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")
    @bar_outer("Fx")
    @bar_outer("Fy")
    #@bar_outer("Sd")
    #@bar_outer("Bd")
    #@bar_outer("Gd")
    @bar_outer("A")
    
    @tilde_outer("Sd")
    @tilde_outer("Bd")
    @tilde_outer("Gd")
    

    @star_outer("A")
    @star_outer("B")
    @star_outer("G")
    @star_outer("S")
    @star_outer("Fx")
    @star_outer("Fy")
    #@star_outer("Sd")
    #@star_outer("Bd")
    #@star_outer("Gd")
   

    @tilde_outer("Bp")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")
    @tilde_outer("Ap")

    @hat_outer("Bp")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")
    @hat_outer("Ap")

    @cross_outer("A")
    @cross_outer("Fx")
    @cross_outer("Fy")
    @cross_outer("S")
    @cross_outer("G")
    

  x0 = S ^ 4
x1 = exp(B)
x2 = x0 * x1
x3 = cosh(G)
x4 = 8 * x3
x5 = x2 * x4
x6 = exp(3 * B)
x7 = x0 * x6
x8 = x4 * x7
x9 = exp(2 * B)
x10 = x0 * x9
x11 = sinh(G)
x12 = 16 * x11
x13 = sech(G)
x14 = x3 ^ 2
x15 = 8 * x14
x16 = Fxp * x2
x17 = x10 * x11
x18 = Fx + xi_x
x19 = Bp * x18
x20 = x15 * x2
x21 = 16 * x14
x22 = Sp * x18
x23 = S ^ 3
x24 = x1 * x23
x25 = x22 * x24
x26 = Gp * x18
x27 = Fy + xi_y
x28 = Gp * x27
x29 = x10 * x15
x30 = Sp * x27
x31 = 16 * x3
x32 = x31 * x9
x33 = x11 * x23
x34 = x32 * x33
x35 = Bt + x19
x36 = Gt + x26
x37 = x11 * x36
x38 = Gh + x28
x39 = x15 * x7
x40 = Bp * x27
x41 = Bh + x40
x42 = x11 * x38
x43 = S * x31
x44 = St * x43
x45 = 16 * Fx
x46 = St * x45
x47 = Sp * x3
x48 = 16 * xi_x
x49 = St * x48
x50 = S ^ 2
x51 = x4 * x50
x52 = 32 * x1
x53 = S * x11
x54 = S * x3
x55 = Bp * x54
x56 = Sp * x54
x57 = Bt * x56
x58 = Fxp * x56
x59 = Bp * Sd
x60 = x23 * x59
x61 = x3 * x50
x62 = Bpt * x61
x63 = 24 * xi_x
x64 = Fxp * x51
x65 = x12 * x50
x66 = x1 * x50
x67 = Sd * Sp
x68 = Fx * Sph
x69 = S * x1
x70 = x12 * x69
x71 = Fy * Spt
x72 = Gh * x1
x73 = Gt * x1
x74 = x43 * x73
x75 = x1 * x53
x76 = Bp * Fx
x77 = Sp * xi_x
x78 = 32 * x54
x79 = Fxp ^ 2
x80 = x11 * x50
x81 = 8 * x1 * x80
x82 = x31 * x66
x83 = 8 * Fx
x84 = Bp * Fxp * x61
x85 = Bp * x80
x86 = Gt * x85
x87 = Fx ^ 2
x88 = Bp * Sp
x89 = x43 * x88
x90 = 8 * xi_x
x91 = xi_x ^ 2
x92 = Bpp * Fx
x93 = Fx * Gp
x94 = 24 * x93
x95 = Bt * x80
x96 = Gp * xi_x
x97 = 24 * x96
x98 = 64 * Fy
x99 = Fx * Spp
x100 = x75 * x99
x101 = x1 * x54
x102 = Gh * x101
x103 = Gp * Sh * x101
x104 = 64 * xi_y
x105 = Fy * Gp
x106 = x1 * x44
x107 = Fy * Sp
x108 = x98 * xi_x
x109 = Spp * x75
x110 = Gp * xi_y
x111 = Sp * xi_y
x112 = x104 * xi_x
x113 = x31 * x50
x114 = Bpp * x91
x115 = Bp * x66
x116 = x115 * x12
x117 = x1 * x80
x118 = Fypp * x117
x119 = Fx * Gph
x120 = x1 * x3
x121 = x120 * x50
x122 = 24 * x121
x123 = x1 * x51
x124 = Gp * x123
x125 = Fxp * Fyp
x126 = Fxpp * x81
x127 = Fy * Gpt
x128 = Sh * x32
x129 = Bp * xi_x
x130 = Bp * Gp
x131 = Bp ^ 2
x132 = Gp ^ 2
x133 = x1 * x85
x134 = Fyp * x133
x135 = Bp * Gh * x121
x136 = Bp * Fy
x137 = Fxp * x66
x138 = x12 * x137
x139 = Bp * xi_y
x140 = Gt * x122
x141 = Fx * Gpp
x142 = x121 * x141
x143 = Fyp * x121
x144 = x72 * x80
x145 = Fxp * x122
x146 = S * x32
x147 = Fyp * x146
x148 = 24 * x73
x149 = Gpp * x121
x150 = Fy ^ 2
x151 = xi_y ^ 2
x152 = 32 * x121
x153 = x152 * x93
x154 = x152 * x96
x155 = Fyp ^ 2
x156 = 3 * G
x157 = cosh(x156)
x158 = x146 * x88
x159 = 2 * Fy
x160 = Fx * x117
x161 = x131 * x160
x162 = 2 * xi_y
x163 = x117 * xi_x
x164 = x131 * x163
x165 = 40 * x132
x166 = x160 * x165
x167 = x163 * x165
x168 = x131 * x50
x169 = x1 * x168 * sinh(x156)
x170 = Fx * x169
x171 = x169 * xi_x
x172 = 2 * x120
x173 = Sp * x120
x174 = 3 * x75
x175 = x101 * x40
x176 = S ^ 6 * x9
x177 = S ^ 5
x178 = 4 * x176
x179 = 8 * x17
x180 = 8 * x9
x181 = A * x180
x182 = A * x10
x183 = 4 * x17
x184 = 12 * x3
x185 = x18 ^ 2
x186 = x185 * x2
x187 = x18 ^ 4
x188 = x3 ^ 4
x189 = 4 * Sp
x190 = Bp * S
x191 = x189 * x190
x192 = x18 ^ 3
x193 = 2 * A
x194 = x17 * x193
x195 = Bpp * x50
x196 = 2 * x187
x197 = 4 * x3
x198 = A * x197
x199 = x16 * x19
x200 = A * x120
x201 = 4 * x23
x202 = x200 * x201
x203 = x27 ^ 2
x204 = x203 * x7
x205 = 2 * x188
x206 = Bp * x192 * x50
x207 = x205 * x206
x208 = 2 * G
x209 = sinh(x208)
x210 = sinh(4 * G)
x211 = x209 ^ 2
x212 = A * x3
x213 = 2 * x212
x214 = Bpp * x186 * x213
x215 = x10 * x197
x216 = Ah * x215
x217 = x23 * x6
x218 = x217 * x30
x219 = x186 * x197
x220 = Fxp * x27
x221 = Ap * x183
x222 = Fyp * x18
x223 = x173 * x185 * x201
x224 = x3 ^ 3
x225 = x19 * x224
x226 = At * x215
x227 = x12 * x23 * x9
x228 = Sd * x227
x229 = 2 * x192
x230 = x209 * x69
x231 = Gd * Gp
x232 = Fxpp * x18 + Fxpt
x233 = x197 * x2
x234 = x233 * x35
x235 = x11 ^ 2
x236 = x15 * x9
x237 = x27 * x50
x238 = x180 * x235
x239 = x18 * x50
x240 = x209 * x66
x241 = x11 * x186
x242 = x130 * x193 * x241
x243 = A * Bp
x244 = x223 * x243
x245 = x10 * x213
x246 = Fxp * x28
x247 = Fyp * x26
x248 = x11 * x185
x249 = x189 * x248
x250 = Fxp * x18
x251 = Fxt + x250 + xi_xx
x252 = 8 * x2
x253 = x11 * x251
x254 = Sd * x4
x255 = x203 * x217
x256 = x255 * x47
x257 = App * x27
x258 = x179 * x18
x259 = x186 * x224
x260 = 4 * Bp
x261 = Bd * x260
x262 = x11 ^ 3
x263 = 8 * Gd
x264 = Bp * x263
x265 = x11 * x204
x266 = 4 * x50
x267 = x11 * x266
x268 = x224 * x267
x269 = exp(4 * B)
x270 = x27 ^ 4
x271 = x269 * x270
x272 = x1 * x185
x273 = 8 * x224
x274 = Fxp * S
x275 = x27 ^ 3
x276 = x275 * x6
x277 = 2 * x276
x278 = x14 * x192
x279 = x278 * x66
x280 = 2 * Gp
x281 = Fxh + x220 + xi_xy
x282 = Fyt + x222 + xi_xy
x283 = Spp * x18 + Spt
x284 = St + x22
x285 = x278 * x284
x286 = x284 ^ 2
x287 = Fxph + Fxpp * x27
x288 = Fypp * x27
x289 = Fyph + x288
x290 = x35 ^ 2
x291 = Fypp * x18 + Fypt
x292 = x36 ^ 2
x293 = 2 * x14
x294 = Sdp * x18 + Sdt
x295 = App * x18 + Apt
x296 = Bdp * x18 + Bdt
x297 = A * x131
x298 = x209 * x50
x299 = x211 * x271
x300 = x298 * x6
x301 = x18 * x27
x302 = x193 * x224 * x7
x303 = x189 * x255
x304 = x10 * x4
x305 = x227 * x27
x306 = 6 * x209
x307 = x306 * x69
x308 = x18 * x23 * x32
x309 = Fyp * x27
x310 = Fyh + x309 + xi_yy
x311 = x197 * x7
x312 = x310 * x7
x313 = Bd * x31
x314 = 4 * S
x315 = Bp * x314
x316 = x14 * x185
x317 = 8 * x316
x318 = x274 * x284
x319 = x189 * x278
x320 = x217 * x254
x321 = kappa * x8
x322 = Bpp * x203
x323 = Fxp * x38
x324 = x245 * x36
x325 = Gpp * x18 + Gpt
x326 = x204 * x224
x327 = x204 * x262
x328 = 8 * Bd * Gp
x329 = Bd * Sp
x330 = x269 * x50
x331 = x210 / 2
x332 = Bpp * x27
x333 = x14 * x50
x334 = x276 * x333
x335 = 2 * x334
x336 = x18 * x24
x337 = 8 * Sd
x338 = Bd * Gp
x339 = x11 * x2 * x317
x340 = cosh(x208)
x341 = x235 * x4
x342 = Bpp * x18
x343 = Bpt + x342
x344 = Sh + x30
x345 = x344 ^ 2
x346 = x266 * x6
x347 = x41 ^ 2
x348 = x38 ^ 2
x349 = x180 * x281
x350 = x220 * x349 * x50
x351 = 3 * x9
x352 = x282 * x50
x353 = x316 * x50
x354 = Sdh + Sdp * x27
x355 = x18 * x227
x356 = Aph + x257
x357 = Bdh + Bdp * x27
x358 = x27 * x357
x359 = A * x183
x360 = Bp * x281
x361 = Gp * x282
x362 = Sph + Spp * x27
x363 = 4 * x11
x364 = A * x2
x365 = 2 * x364
x366 = x194 * x301
x367 = x1 * x266
x368 = x235 * x367
x369 = Bp * x192 * x368
x370 = x28 * x369
x371 = x266 * x9
x372 = x185 * x368
x373 = x203 * x6
x374 = x306 * x373
x375 = x14 * x194
x376 = Fxp * x40
x377 = x168 * x211
x378 = x198 * x217
x379 = Gph + Gpp * x27
x380 = 2 * x379
x381 = x10 * x212
x382 = x18 * x343
x383 = 2 * x325
x384 = x27 * x383
x385 = x27 * x7
x386 = Bd * x38
x387 = x211 * x275
x388 = Bp * x387
x389 = Fyp / 2
x390 = x192 * x3
x391 = 2 * x3
x392 = x185 * x391
x393 = Fxp * x36 * x80
x394 = 16 * x9
x395 = x27 * x284
x396 = x23 * x395
x397 = x18 * x41
x398 = Gp * x192
x399 = x391 * x80
x400 = x35 * x399
x401 = x27 * x6
x402 = x27 * x41
x403 = x217 * x27
x404 = x42 * x7
x405 = 8 * x404
x406 = x120 * x266
x407 = x262 * x406
x408 = 2 * x353
x409 = Bp * x251
x410 = x269 * x275
x411 = x172 * x80
x412 = x185 * x411
x413 = 4 * A
x414 = 6 * x38
x415 = x324 * x40
x416 = x194 * x38
x417 = x11 * x344
x418 = x203 * x9
x419 = x19 * x50
x420 = x205 * x418 * x419
x421 = Fyp * x203
x422 = x235 * x346
x423 = x421 * x422
x424 = x14 * x203
x425 = x424 * x50
x426 = 2 * x425
x427 = x426 * x6
x428 = x101 * x249
x429 = x276 * x53
x430 = x284 * x35
x431 = Bph + x332
x432 = 4 * x37
x433 = Fyp * x316
x434 = S * x344
x435 = x269 * x421
x436 = 2 * x36
x437 = x276 * x36
x438 = x185 * x27
x439 = x11 * x18
x440 = Bd * x35
x441 = Bd * x209
x442 = x11 * x27
x443 = x131 * x224
x444 = x275 * x346
x445 = x185 * x35
x446 = x120 * x267
x447 = x301 * x446
x448 = x19 * x209
x449 = 8 * x101
x450 = x330 * x431
x451 = x205 * x275
x452 = 6 * x9
x453 = Gt + x93 + x96
x454 = x236 * x50
x455 = 6 * Gp
x456 = x316 * x66
x457 = x343 * x80
x458 = x373 * x399
x459 = x11 * x209 * x7
x460 = 4 * x459
x461 = x209 * x27
x462 = A * x459
x463 = x185 * x9
x464 = Sh + x107 + x111
x465 = x224 * x260
x466 = x344 * x9
x467 = S * x18
x468 = 16 * x235 * x466 * x467
x469 = S * x394
x470 = x395 * x469
x471 = x185 * x36
x472 = 4 * x316
x473 = x11 * x54
x474 = x189 * x373 * x473
x475 = x364 * x391
x476 = x198 * x23
x477 = x36 * x38
x478 = x330 * x41
x479 = x408 * x9
x480 = Fyp * x41
x481 = 2 * x330 * x424
x482 = x185 * x54
x483 = x235 * x301
x484 = x371 * x483
x485 = Bph + Bpp * Fy + Bpp * xi_y
x486 = x117 * x224
x487 = x203 * x284
x488 = x30 * x344
x489 = x83 + x90
x490 = Bd * x27
x491 = x392 * x80
x492 = x491 * x9
x493 = x11 * x3
x494 = x346 * x493
x495 = x18 * x28
x496 = x40 * x41
x497 = x426 * x9
x498 = x235 * x281
x499 = x203 * x346
x500 = x424 * x6
x501 = x224 * x277
x502 = x4 * x417
x503 = x314 * x424
x504 = x503 * x9
x505 = Fx * Sp + St + x77
x506 = Bh + x136 + x139
x507 = 8 * x373
x508 = x38 * x507
x509 = x209 * x40
x510 = x19 * x41
x511 = x19 * x507 * x54
x512 = x224 * x40
x513 = x3 * x42
x514 = x3 * x37 * x371
x515 = x205 * x463
x516 = x203 * x422
x517 = Bpp * xi_x + Bpt + x92
x518 = x487 * x9
x519 = x19 * x27
x520 = x373 * x51
x521 = S * x21
x522 = Gh + x105 + x110
x523 = Bt + x129 + x76
x524 = x373 * x80
x525 = x11 * x301
x526 = S * x472
x527 = x526 * x9
x528 = x344 * x41
x529 = x269 * x503
x530 = x35 * x38
x531 = x203 * x269
x532 = x185 * x40
x533 = x523 * x532
x534 = x18 * x485
x535 = x262 * x3 * x499
x536 = x180 * x301
x537 = x41 * x513
x538 = x51 * x9
x539 = x159 + x162
x540 = Spp * x203 + Ss + x362 * x539
x541 = 2 * Spt
x542 = Fx * (x541 + x99) + Sb + Spp * x91 + xi_x * (x541 + 2 * x99)
x543 = Bs + x322 + x431 * x539
x544 = Gpp * x203 + Gs + x379 * x539
x545 = 2 * Bpt
x546 = Bb + Fx * (x545 + x92) + x114 + xi_x * (x545 + 2 * x92)
x547 = 2 * Gpt
x548 = Fx * (x141 + x547) + Gb + Gpp * x91 + xi_x * (2 * x141 + x547)
x549 = x203 * x399
x550 = 3 * x99
x551 = 3 * Spp
x552 = Fy * x550 + Sc + x68 + x71 + xi_x * (Fy * x551 + Sph + x551 * xi_y) + xi_y * (Spt + x550)
x553 = 3 * x141
x554 = 3 * Gpp
x555 = Fy * x553 + Gc + x119 + x127 + xi_x * (Fy * x554 + Gph + x554 * xi_y) + xi_y * (Gpt + x553)
axx = -x5
ayy = -x8
axy = x10 * x12
bx = x13 * (Fyp * x17 * x4 - x11 * x26 * x5 - x15 * x16 + x19 * x20 + x20 * x35 + x21 * x25 + x28 * x29 + x29 * x38 - x30 * x34 - x37 * x5)
by = x13 * (8 * Fxp * x0 * x11 * x3 * x9 - Fyp * x39 + 8 * Gp * x0 * x14 * x18 * x9 + 16 * Sp * x14 * x23 * x27 * x6 + 8 * x0 * x14 * x36 * x9 - x11 * x28 * x8 - x22 * x34 - x39 * x40 - x39 * x41 - x42 * x8)
cc = -x66 * (-Bb * x51 + 8 * Bh ^ 2 * x3 * x50 * x9 + 8 * Bh * S * x1 * (Fyp * x101 + 2 * Gh * x75 - Gt * x54 + Sh * x172 + x105 * x174 + x110 * x174 + x159 * x173 + x162 * x173 + 5 * x175 - x54 * x93 - x54 * x96) + 8 * Bp * Bt * Fx * x3 * x50 + 8 * Bp * Bt * x3 * x50 * xi_x + 32 * Bp * Fx * Fy * S * Sp * x1 * x11 + 32 * Bp * Fx * S * Sp * x1 * x11 * xi_y + 8 * Bp * Fxt * x3 * x50 + 40 * Bp * Fy * Fyp * x3 * x50 * x9 + 40 * Bp * Fy * Gh * x11 * x50 * x9 + 96 * Bp * Fy * Gp * x11 * x50 * x9 * xi_y + 16 * Bp * Fy * S * Sh * x3 * x9 + 32 * Bp * Fy * S * Sp * x1 * x11 * xi_x + 24 * Bp * Fyh * x3 * x50 * x9 + 40 * Bp * Fyp * x3 * x50 * x9 * xi_y + 40 * Bp * Gh * x11 * x50 * x9 * xi_y + 48 * Bp * Gp * x11 * x150 * x50 * x9 + 48 * Bp * Gp * x11 * x151 * x50 * x9 + 16 * Bp * S * Sh * x3 * x9 * xi_y + 32 * Bp * S * Sp * x1 * x11 * xi_x * xi_y + 24 * Bp * x3 * x50 * x9 * xi_yy + 8 * Bp * x3 * x50 * xi_xx + 24 * Bph * Fy * x3 * x50 * x9 + 24 * Bph * x3 * x50 * x9 * xi_y + 64 * Bpp * Fy * x3 * x50 * x9 * xi_y - Bpp * x113 * x87 + 32 * Bpp * x150 * x3 * x50 * x9 + 32 * Bpp * x151 * x3 * x50 * x9 + 8 * Bs * x3 * x50 * x9 + 8 * Bt ^ 2 * x3 * x50 + 8 * Bt * Fy * Gp * x1 * x3 * x50 + 8 * Bt * Gh * x1 * x3 * x50 + 8 * Bt * Gp * x1 * x3 * x50 * xi_y - Bt * Gt * x65 - Bt * x44 - Bt * x64 + 24 * Fx * Fxp * Gp * x11 * x50 + 8 * Fx * Fxpp * x3 * x50 + 16 * Fx * Fyp * S * Sp * x1 * x11 + 24 * Fx * Gp * Gt * x3 * x50 + 16 * Fx * Gp * S * St * x11 + 32 * Fx * Gpp * x11 * x50 * xi_x + 24 * Fx * Gpt * x11 * x50 + 16 * Fx * Gt * S * Sp * x11 + 16 * Fx * S * Spt * x3 + 16 * Fx * Sh * Sp * x1 * x11 + 2 * Fx * x131 * x157 * x50 * xi_x + 6 * Fx * x131 * x3 * x50 * xi_x + 40 * Fx * x132 * x3 * x50 * xi_x - 24 * Fx * x62 - Fxh * x116 - Fxh * x124 + 16 * Fxp * Fy * S * Sp * x1 * x11 + 24 * Fxp * Gp * x11 * x50 * xi_x + 8 * Fxp * Gt * x11 * x50 + 16 * Fxp * S * Sp * x1 * x11 * xi_y - Fxph * x81 + 8 * Fxpp * x3 * x50 * xi_x + 8 * Fxpt * x3 * x50 + 8 * Fxt * Gp * x11 * x50 + 24 * Fy * Fyp * Gp * x11 * x50 * x9 + 8 * Fy * Fypp * x3 * x50 * x9 + 24 * Fy * Gh * Gp * x3 * x50 * x9 + 16 * Fy * Gh * S * Sp * x11 * x9 + 16 * Fy * Gp * S * Sh * x11 * x9 + 24 * Fy * Gph * x11 * x50 * x9 + 64 * Fy * Gpp * x11 * x50 * x9 * xi_y + 16 * Fy * S * Sph * x3 * x9 + 64 * Fy * S * Spp * x3 * x9 * xi_y + 16 * Fy * Sp * St * x1 * x11 - Fy * x126 + 2 * Fy * x131 * x157 * x50 * x9 * xi_y + 70 * Fy * x131 * x3 * x50 * x9 * xi_y + 40 * Fy * x132 * x3 * x50 * x9 * xi_y - Fy * x166 - Fy * x167 + 8 * Fyh * Gp * x11 * x50 * x9 + 8 * Fyp * Gh * x11 * x50 * x9 + 24 * Fyp * Gp * x11 * x50 * x9 * xi_y + 16 * Fyp * S * Sp * x1 * x11 * xi_x - Fyp * x51 * x73 + 8 * Fyph * x3 * x50 * x9 + 8 * Fypp * x3 * x50 * x9 * xi_y - Fypt * x81 - Fyt * x116 - Fyt * x124 + 8 * Gb * x11 * x50 - Gc * x82 + 8 * Gh ^ 2 * x3 * x50 * x9 + 24 * Gh * Gp * x3 * x50 * x9 * xi_y - Gh * Gt * x12 * x66 + 16 * Gh * S * Sh * x11 * x9 + 16 * Gh * S * Sp * x11 * x9 * xi_y + 24 * Gp * Gt * x3 * x50 * xi_x + 16 * Gp * S * Sh * x11 * x9 * xi_y + 16 * Gp * S * St * x11 * xi_x + 8 * Gp * x11 * x50 * x9 * xi_yy + 8 * Gp * x11 * x50 * xi_xx - Gp * x82 * xi_xy + 24 * Gph * x11 * x50 * x9 * xi_y - Gph * x121 * x63 + 32 * Gpp * x11 * x150 * x50 * x9 + 32 * Gpp * x11 * x151 * x50 * x9 + 16 * Gpp * x11 * x50 * x87 + 16 * Gpp * x11 * x50 * x91 + 24 * Gpt * x11 * x50 * xi_x - Gpt * x122 * xi_y + 8 * Gs * x11 * x50 * x9 + 8 * Gt ^ 2 * x3 * x50 + 16 * Gt * S * Sp * x11 * xi_x + 16 * Gt * S * St * x11 + 16 * S * Sb * x3 + 16 * S * Sph * x3 * x9 * xi_y + 32 * S * Spp * x150 * x3 * x9 + 32 * S * Spp * x151 * x3 * x9 + 16 * S * Spt * x3 * xi_x + 16 * S * Ss * x3 * x9 - Sc * x52 * x53 - Sh ^ 2 * x32 + 16 * Sh * Sp * x1 * x11 * xi_x + 32 * Sh * St * x1 * x11 - Sh * x74 + 16 * Sp * St * x1 * x11 * xi_y - Sp * x102 * x45 - Sp * x102 * x48 - Sph * x48 * x75 - Spt * x70 * xi_y - St ^ 2 * x31 - x100 * x104 - x100 * x98 - x103 * x45 - x103 * x48 - x104 * x142 - x105 * x106 - x105 * x145 - x105 * x148 * x80 - x106 * x110 - x107 * x128 - x107 * x139 * x78 * x9 - x107 * x147 - x107 * x74 - x108 * x109 - x108 * x149 - x109 * x112 - x110 * x145 - x110 * x148 * x80 - x111 * x128 - x111 * x147 - x111 * x74 - x112 * x149 - x113 * x114 - x118 * x83 - x118 * x90 - x119 * x122 - x122 * x127 - x125 * x81 - x126 * xi_y - 32 * x129 * x80 * x93 - x130 * x65 * x87 - x130 * x65 * x91 + x131 * x150 * x157 * x50 * x9 + 35 * x131 * x150 * x3 * x50 * x9 + x131 * x151 * x157 * x50 * x9 + 35 * x131 * x151 * x3 * x50 * x9 + x131 * x157 * x50 * x87 + x131 * x157 * x50 * x91 + 3 * x131 * x3 * x50 * x87 + 3 * x131 * x3 * x50 * x91 + 20 * x132 * x150 * x3 * x50 * x9 + 20 * x132 * x151 * x3 * x50 * x9 + 20 * x132 * x3 * x50 * x87 + 20 * x132 * x3 * x50 * x91 - x134 * x45 - x134 * x48 - x135 * x83 - x135 * x90 - x136 * x138 - x136 * x140 - x136 * x153 - x136 * x154 - x138 * x139 - x139 * x140 - x139 * x153 - x139 * x154 - x142 * x98 - x143 * x94 - x143 * x97 - x144 * x94 - x144 * x97 - x150 * x158 - x151 * x158 + 4 * x155 * x3 * x50 * x9 - x159 * x161 - x159 * x164 - x159 * x170 - x159 * x171 - x161 * x162 - x162 * x164 - x162 * x170 - x162 * x171 - x166 * xi_y - x167 * xi_y - 48 * x2 + 4 * x3 * x50 * x79 - x44 * x72 - x45 * x57 - x45 * x58 - x46 * x47 - x46 * x55 - x47 * x49 - x48 * x57 - x48 * x58 - x49 * x55 - x52 * x60 - x52 * x85 * xi_xy - 32 * x61 * x92 * xi_x - x62 * x63 - x64 * x72 - 32 * x66 * x67 - x68 * x70 - x70 * x71 - x76 * x77 * x78 - x83 * x84 - x83 * x86 - x84 * x90 - x86 * x90 - x87 * x89 - x89 * x91 - x94 * x95 - x95 * x97) / 2
SS = 2 * A * Bp * Fxp * x0 * x1 * x18 * x224 + A * Bp * Fxp * x0 * x209 * x27 * x3 * x9 + A * Bp * Fyp * x0 * x11 * x209 * x27 * x6 + 4 * A * Bp * Fyp * x0 * x27 * x3 * x6 + 4 * A * Bp * Gp * x0 * x1 * x11 * x14 * x185 + 2 * A * Bp * Gp * x0 * x11 * x203 * x340 * x6 + 2 * A * Bp * Gp * x0 * x11 * x203 * x6 + 4 * A * Bp * Sp * x1 * x185 * x224 * x23 + 8 * A * Bp * Sp * x11 * x18 * x23 * x27 * x9 + 4 * A * Bp * Sp * x203 * x23 * x235 * x3 * x6 + 2 * A * Bp * x0 * x1 * x11 * x18 * x340 * x36 + 4 * A * Bp * x0 * x1 * x11 * x18 * x36 + 2 * A * Bp * x0 * x1 * x18 * x224 * x35 + 2 * A * Bp * x0 * x1 * x251 * x3 + 2 * A * Bp * x0 * x11 * x209 * x27 * x36 * x9 + 4 * A * Bp * x0 * x11 * x27 * x38 * x6 + A * Bp * x0 * x18 * x209 * x3 * x41 * x9 + 2 * A * Bp * x0 * x18 * x3 * x340 * x38 * x9 + A * Bp * x0 * x209 * x27 * x3 * x35 * x9 + 2 * A * Bp * x0 * x209 * x27 * x3 * x38 * x6 + 2 * A * Bp * x0 * x224 * x27 * x41 * x6 + 4 * A * Bp * x0 * x27 * x3 * x41 * x6 + 6 * A * Bp * x0 * x3 * x310 * x6 + 2 * A * Bp * x1 * x11 * x18 * x209 * x23 * x284 + 4 * A * Bp * x1 * x18 * x23 * x284 * x3 + 4 * A * Bp * x11 * x14 * x23 * x27 * x284 * x9 + 2 * A * Bp * x18 * x209 * x23 * x3 * x344 * x9 + 4 * A * Bp * x224 * x23 * x27 * x344 * x6 + 2 * A * Bpp * x0 * x1 * x185 * x224 + A * Bpp * x0 * x11 * x203 * x209 * x6 + 2 * A * Bpp * x0 * x203 * x3 * x6 + 2 * A * Fxp * Gp * x0 * x1 * x11 * x18 + 4 * A * Fxp * Sp * x11 * x23 * x27 * x9 + 2 * A * Fxp * x0 * x1 * x11 * x36 + 2 * A * Fyp * Gp * x0 * x11 * x27 * x6 + 4 * A * Fyp * Sp * x11 * x18 * x23 * x9 + 2 * A * Fyp * x0 * x11 * x38 * x6 + 2 * A * Fyp * x0 * x3 * x41 * x6 + 8 * A * Gp * Sp * x18 * x23 * x27 * x3 * x9 + 2 * A * Gp * x0 * x1 * x11 * x251 + 2 * A * Gp * x0 * x1 * x18 * x3 * x36 + 2 * A * Gp * x0 * x11 * x27 * x41 * x6 + 2 * A * Gp * x0 * x11 * x310 * x6 + 2 * A * Gp * x0 * x27 * x3 * x38 * x6 - A * Gp * x11 * x303 - A * Gp * x24 * x249 + 4 * A * Sp * x1 * x18 * x284 * x3 * x50 + 4 * A * Sp * x27 * x3 * x344 * x50 * x6 + A * x0 * x1 * x11 * x131 * x185 * x209 + A * x0 * x1 * x11 * x18 * x209 * x343 + 2 * A * x0 * x1 * x11 * x18 * x325 + 2 * A * x0 * x1 * x11 * x548 + 2 * A * x0 * x1 * x131 * x185 * x3 + A * x0 * x1 * x132 * x185 * x3 + 2 * A * x0 * x1 * x232 * x3 + 2 * A * x0 * x1 * x290 * x3 + 2 * A * x0 * x1 * x292 * x3 + A * x0 * x1 * x3 * x79 + A * x0 * x11 * x131 * x203 * x209 * x6 + 2 * A * x0 * x11 * x14 * x27 * x343 * x9 + 2 * A * x0 * x11 * x27 * x379 * x6 + 4 * A * x0 * x11 * x38 * x41 * x6 + 2 * A * x0 * x11 * x544 * x6 + 2 * A * x0 * x131 * x203 * x3 * x6 + A * x0 * x132 * x203 * x3 * x6 + A * x0 * x155 * x3 * x6 + A * x0 * x18 * x209 * x3 * x431 * x9 + 2 * A * x0 * x224 * x27 * x431 * x6 + 2 * A * x0 * x289 * x3 * x6 + 2 * A * x0 * x3 * x347 * x6 + 2 * A * x0 * x3 * x348 * x6 + 2 * A * x0 * x3 * x35 * x38 * x9 + 2 * A * x0 * x3 * x543 * x6 + 4 * A * x1 * x11 * x23 * x284 * x36 + 4 * A * x1 * x23 * x3 * x542 + 4 * A * x11 * x18 * x23 * x362 * x9 + 4 * A * x11 * x23 * x27 * x283 * x9 + 4 * A * x11 * x23 * x344 * x38 * x6 - A * x11 * x284 * x30 * x371 + 8 * A * x11 * x284 * x344 * x50 * x9 - 12 * A * x176 - A * x19 * x234 - A * x22 * x371 * x417 + 4 * A * x23 * x3 * x344 * x41 * x6 + 4 * A * x23 * x3 * x540 * x6 - A * x256 * x260 + 4 * Ab * x0 * x1 * x3 - Ac * x179 + 4 * Ah * Bp * x0 * x224 * x27 * x6 - Ah * Fxp * x183 + 4 * Ah * Fyp * x0 * x3 * x6 + 4 * Ah * Gp * x0 * x11 * x27 * x6 + 8 * Ah * Sp * x11 * x18 * x23 * x9 + 4 * Ah * x0 * x11 * x38 * x6 + 4 * Ah * x0 * x3 * x41 * x6 - Ah * x218 * x4 - Ah * x235 * x311 * x40 + 4 * Ap * Bp * x0 * x203 * x3 * x6 - Ap * Bp * x219 + 4 * Ap * Fxp * x0 * x1 * x18 * x3 + 4 * Ap * Fyp * x0 * x27 * x3 * x6 + 4 * Ap * Gp * x0 * x1 * x11 * x185 + 4 * Ap * Gp * x0 * x11 * x203 * x6 + 8 * Ap * Sd * x177 * x9 + 8 * Ap * Sp * x11 * x18 * x23 * x27 * x9 + 4 * Ap * x0 * x11 * x281 * x9 + 4 * Ap * x0 * x11 * x282 * x9 - Ap * x18 * x28 * x304 - Ap * x223 - Ap * x233 * x251 - 4 * Ap * x256 - Ap * x310 * x311 + 4 * App * x0 * x203 * x3 * x6 - App * x186 * x197 + 4 * As * x0 * x3 * x6 + 4 * At * Bp * x0 * x1 * x18 * x235 * x3 + 4 * At * Fxp * x0 * x1 * x3 - At * Fyp * x183 + 4 * At * Gp * x0 * x1 * x11 * x18 + 8 * At * Sp * x11 * x23 * x27 * x9 + 4 * At * x0 * x1 * x11 * x36 - 4 * At * x2 * x225 - At * x234 - At * x25 * x4 - Bd ^ 2 * x14 * x178 + 8 * Bd * Bp * x0 * x11 * x14 * x18 * x27 * x9 + 8 * Bd * Gp * x0 * x1 * x185 * x209 * x3 + 8 * Bd * Gp * x0 * x1 * x185 * x262 + 8 * Bd * Gp * x0 * x11 * x14 * x203 * x6 + 8 * Bd * Gp * x0 * x11 * x203 * x340 * x6 + 16 * Bd * Sd * x177 * x9 + 8 * Bd * Sp * x1 * x185 * x224 * x23 + 8 * Bd * Sp * x203 * x23 * x235 * x3 * x6 + 16 * Bd * Sp * x203 * x23 * x3 * x6 + 8 * Bd * x0 * x1 * x11 * x18 * x340 * x36 + 8 * Bd * x0 * x1 * x18 * x224 * x35 + 8 * Bd * x0 * x11 * x209 * x27 * x36 * x9 + 8 * Bd * x0 * x11 * x281 * x9 + 8 * Bd * x0 * x11 * x282 * x9 + 4 * Bd * x0 * x18 * x209 * x3 * x41 * x9 + 8 * Bd * x0 * x18 * x3 * x340 * x38 * x9 + 8 * Bd * x0 * x18 * x3 * x38 * x9 + 4 * Bd * x0 * x209 * x27 * x3 * x35 * x9 + 8 * Bd * x0 * x209 * x27 * x3 * x38 * x6 + 8 * Bd * x0 * x224 * x27 * x41 * x6 + 8 * Bd * x0 * x27 * x3 * x36 * x9 + 8 * Bd * x1 * x11 * x18 * x209 * x23 * x284 + 16 * Bd * x11 * x14 * x23 * x27 * x284 * x9 - Bd * x11 * x29 * x397 - Bd * x14 * x344 * x355 + 8 * Bd * x18 * x209 * x23 * x3 * x344 * x9 - Bd * x22 * x305 + 16 * Bd * x224 * x23 * x27 * x344 * x6 - 16 * Bd * x224 * x284 * x336 - Bd * x402 * x460 + 4 * Bp * Fxp * x1 * x11 * x185 * x224 * x27 * x50 + Bp * Fxp * x11 * x18 * x203 * x209 * x3 * x50 * x9 + Bp * Fxp * x192 * x211 * x50 / 2 + 4 * Bp * Fyp * x18 * x203 * x262 * x3 * x50 * x6 + 2 * Bp * Fyp * x185 * x188 * x27 * x50 * x9 + 2 * Bp * Fyp * x188 * x269 * x275 * x50 + 8 * Bp * Gd * x0 * x1 * x11 * x185 + 8 * Bp * Gd * x0 * x1 * x185 * x262 + 8 * Bp * Gd * x0 * x11 * x14 * x203 * x6 - Bp * Gd * x339 + 2 * Bp * Gp * x1 * x192 * x211 * x27 * x50 + 4 * Bp * Gp * x11 * x224 * x269 * x270 * x50 + 4 * Bp * Gp * x18 * x235 * x275 * x340 * x50 * x6 + 4 * Bp * Gp * x18 * x235 * x275 * x50 * x6 + Bp * Gp * x187 * x209 * x50 + Bp * Gp * x187 * x210 * x50 / 2 + 8 * Bp * S * Sp * x1 * x11 * x192 * x224 * x27 + 4 * Bp * S * Sp * x14 * x187 + 4 * Bp * S * Sp * x18 * x209 * x275 * x6 + 8 * Bp * S * Sp * x18 * x262 * x275 * x3 * x6 + Bp * S * Sp * x187 * x211 + 4 * Bp * S * Sp * x188 * x269 * x270 + 4 * Bp * S * x1 * x11 * x185 * x224 * x27 * x505 + 8 * Bp * S * x1 * x11 * x185 * x27 * x284 * x3 + 4 * Bp * S * x1 * x11 * x192 * x224 * x344 + 8 * Bp * S * x1 * x185 * x262 * x27 * x3 * x505 + 12 * Bp * S * x11 * x18 * x203 * x224 * x344 * x6 + 4 * Bp * S * x11 * x224 * x275 * x505 * x6 + 4 * Bp * S * x14 * x185 * x27 * x344 * x9 + 4 * Bp * S * x14 * x269 * x275 * x344 + 4 * Bp * S * x18 * x188 * x203 * x284 * x9 + Bp * S * x185 * x211 * x27 * x344 * x9 + 4 * Bp * S * x188 * x192 * x284 + Bp * S * x211 * x269 * x275 * x344 + 8 * Bp * Sd * x1 * x185 * x23 * x235 * x3 + 8 * Bp * Sd * x203 * x224 * x23 * x6 + 4 * Bp * x1 * x11 * x18 * x251 * x27 * x3 * x50 + 6 * Bp * x1 * x11 * x185 * x224 * x27 * x35 * x50 + 4 * Bp * x1 * x11 * x185 * x282 * x3 * x50 + 2 * Bp * x1 * x11 * x192 * x224 * x41 * x50 + 2 * Bp * x1 * x14 * x185 * x27 * x340 * x36 * x50 + 2 * Bp * x1 * x14 * x192 * x38 * x50 + 4 * Bp * x1 * x185 * x235 * x27 * x340 * x453 * x50 + 8 * Bp * x1 * x185 * x235 * x27 * x36 * x50 + Bp * x1 * x192 * x211 * x38 * x50 + 4 * Bp * x1 * x192 * x235 * x340 * x50 * x522 + Bp * x11 * x18 * x203 * x209 * x3 * x35 * x50 * x9 + 6 * Bp * x11 * x18 * x203 * x224 * x41 * x50 * x6 + 4 * Bp * x11 * x18 * x203 * x224 * x453 * x50 * x9 + Bp * x11 * x185 * x209 * x27 * x3 * x41 * x50 * x9 + 3 * Bp * x11 * x192 * x209 * x3 * x35 * x50 + 4 * Bp * x11 * x192 * x224 * x453 * x50 + 4 * Bp * x11 * x203 * x282 * x3 * x50 * x6 + 3 * Bp * x11 * x209 * x269 * x275 * x3 * x41 * x50 + 2 * Bp * x11 * x224 * x275 * x35 * x50 * x6 + 4 * Bp * x11 * x269 * x275 * x3 * x38 * x50 + 4 * Bp * x14 * x18 * x203 * x35 * x50 * x9 + 2 * Bp * x14 * x18 * x203 * x38 * x50 * x6 + 8 * Bp * x14 * x18 * x27 * x281 * x50 * x9 + 4 * Bp * x14 * x185 * x27 * x41 * x50 * x9 + 2 * Bp * x14 * x185 * x310 * x50 * x9 + 4 * Bp * x14 * x192 * x35 * x50 + 2 * Bp * x14 * x203 * x269 * x310 * x50 + 4 * Bp * x14 * x269 * x275 * x41 * x50 + 2 * Bp * x14 * x275 * x340 * x36 * x50 * x6 + 3 * Bp * x18 * x203 * x211 * x38 * x50 * x6 + 8 * Bp * x18 * x203 * x262 * x3 * x453 * x50 * x9 + (3 // 2) * Bp * x185 * x210 * x27 * x38 * x50 * x9 + Bp * x210 * x269 * x275 * x38 * x50 / 2 + 4 * Bp * x235 * x275 * x340 * x36 * x50 * x6 - Bp * x266 * x37 * x390 - Bp * x268 * x410 * x522 - Bp * x275 * x340 * x422 * x453 - Bp * x282 * x359 - Bp * x335 * x36 - Bp * x451 * x478 + 4 * Bpp * x1 * x11 * x192 * x224 * x27 * x50 + 2 * Bpp * x14 * x187 * x50 + 2 * Bpp * x18 * x209 * x275 * x50 * x6 + 4 * Bpp * x18 * x262 * x275 * x3 * x50 * x6 + Bpp * x187 * x211 * x50 / 2 + 2 * Bpp * x188 * x269 * x270 * x50 + 8 * Fxc * x235 * x27 * x50 * x9 - Fxc * x236 * x237 + 8 * Fxp * Fyp * x14 * x18 * x27 * x50 * x9 + 3 * Fxp * Gp * x18 * x203 * x209 * x50 * x9 + Fxp * Gp * x192 * x209 * x50 - Fxp * Gp * x335 + 4 * Fxp * S * Sp * x14 * x18 * x203 * x9 + 4 * Fxp * S * Sp * x14 * x192 + 8 * Fxp * S * Sp * x18 * x203 * x235 * x9 + 8 * Fxp * S * x1 * x11 * x185 * x3 * x344 + 8 * Fxp * S * x1 * x18 * x209 * x27 * x284 + 8 * Fxp * S * x11 * x203 * x3 * x344 * x6 + 16 * Fxp * Sd * x1 * x18 * x23 * x3 + 4 * Fxp * x1 * x11 * x185 * x3 * x41 * x50 + 6 * Fxp * x1 * x14 * x185 * x38 * x50 + 4 * Fxp * x1 * x18 * x235 * x27 * x36 * x50 + 4 * Fxp * x11 * x203 * x3 * x41 * x50 * x6 + 2 * Fxp * x14 * x185 * x35 * x50 + 2 * Fxp * x14 * x203 * x35 * x50 * x9 + 8 * Fxp * x14 * x27 * x281 * x50 * x9 + 8 * Fxp * x18 * x235 * x310 * x50 * x9 - Fxp * x202 * x22 + 4 * Fxp * x203 * x235 * x38 * x50 * x6 - Fxp * x207 + 8 * Fxp * x27 * x282 * x50 * x9 - Fxp * x420 + 8 * Fxs * x14 * x18 * x50 * x9 - Fxs * x238 * x239 + 8 * Fyb * x14 * x27 * x50 * x9 - Fyb * x237 * x238 + 8 * Fyc * x18 * x235 * x50 * x9 - Fyc * x236 * x239 + 3 * Fyp * Gp * x185 * x209 * x27 * x50 * x9 + Fyp * Gp * x209 * x269 * x275 * x50 + 4 * Fyp * S * Sp * x14 * x185 * x27 * x9 + 4 * Fyp * S * Sp * x14 * x269 * x275 + 8 * Fyp * S * Sp * x185 * x235 * x27 * x9 + 8 * Fyp * S * x1 * x11 * x185 * x284 * x3 + 8 * Fyp * S * x11 * x203 * x284 * x3 * x6 + 8 * Fyp * S * x18 * x209 * x27 * x344 * x6 + 16 * Fyp * Sd * x23 * x27 * x3 * x6 - Fyp * Sp * x229 * x230 + 4 * Fyp * x1 * x185 * x235 * x36 * x50 + 4 * Fyp * x11 * x18 * x27 * x3 * x41 * x50 * x6 + 8 * Fyp * x14 * x18 * x27 * x35 * x50 * x9 + 8 * Fyp * x14 * x18 * x282 * x50 * x9 + 6 * Fyp * x14 * x203 * x36 * x50 * x6 + 4 * Fyp * x18 * x235 * x27 * x38 * x50 * x6 + 8 * Fyp * x18 * x281 * x50 * x9 - Fyp * x198 * x218 + 8 * Fyp * x235 * x251 * x27 * x50 * x9 - Fyp * x279 * x280 - Fyp * x302 * x40 - Fyp * x324 - Fyp * x38 * x492 - Fyp * x445 * x446 + 16 * Fypp * x185 * x235 * x27 * x50 * x9 - Gd ^ 2 * x178 + 8 * Gd * Gp * x0 * x11 * x18 * x27 * x9 + 8 * Gd * x0 * x27 * x3 * x35 * x9 + 8 * Gd * x0 * x281 * x3 * x9 + 8 * Gd * x0 * x282 * x3 * x9 + 16 * Gd * x1 * x11 * x18 * x23 * x284 - Gd * x10 * x397 * x4 + 16 * Gd * x11 * x23 * x27 * x344 * x6 - Gd * x252 * x253 - Gd * x308 * x344 - Gd * x32 * x396 + 8 * Gp * Sd * x1 * x11 * x185 * x23 + 8 * Gp * Sd * x11 * x203 * x23 * x6 + 2 * Gp * x1 * x14 * x185 * x282 * x50 + 4 * Gp * x1 * x185 * x235 * x27 * x35 * x50 + 4 * Gp * x1 * x185 * x235 * x281 * x50 - Gp * x10 * x213 * x281 + 4 * Gp * x11 * x18 * x27 * x281 * x3 * x50 * x9 + 4 * Gp * x11 * x18 * x27 * x282 * x3 * x50 * x9 + 2 * Gp * x14 * x18 * x203 * x36 * x50 * x9 + 2 * Gp * x14 * x185 * x27 * x38 * x50 * x9 + 2 * Gp * x14 * x192 * x36 * x50 + 2 * Gp * x14 * x203 * x281 * x50 * x6 + 2 * Gp * x14 * x269 * x275 * x38 * x50 + 4 * Gp * x18 * x203 * x235 * x36 * x50 * x9 + Gp * x185 * x209 * x251 * x50 + Gp * x185 * x209 * x27 * x41 * x50 * x9 + Gp * x185 * x209 * x310 * x50 * x9 + 4 * Gp * x185 * x235 * x27 * x38 * x50 * x9 + Gp * x203 * x209 * x251 * x50 * x9 + Gp * x203 * x209 * x269 * x310 * x50 + 4 * Gp * x203 * x235 * x282 * x50 * x6 + Gp * x209 * x269 * x275 * x41 * x50 - Gp * x399 * x437 - Gp * x498 * x499 + 8 * S * Sp * x1 * x11 * x185 * x27 * x3 * x35 + 4 * S * Sp * x11 * x192 * x3 * x36 + 4 * S * Sp * x11 * x269 * x275 * x3 * x38 + 4 * S * Sp * x14 * x185 * x251 + 4 * S * Sp * x14 * x185 * x27 * x41 * x9 + 4 * S * Sp * x14 * x185 * x310 * x9 + 4 * S * Sp * x14 * x203 * x251 * x9 + 4 * S * Sp * x14 * x203 * x269 * x310 + 4 * S * Sp * x14 * x269 * x275 * x41 + 6 * S * Sp * x18 * x203 * x209 * x36 * x9 + 8 * S * Sp * x18 * x235 * x27 * x281 * x9 + 8 * S * Sp * x18 * x235 * x27 * x282 * x9 + 6 * S * Sp * x185 * x209 * x27 * x38 * x9 - S * Sp * x222 * x374 + 8 * S * x1 * x11 * x18 * x27 * x3 * x542 + 8 * S * x1 * x11 * x185 * x3 * x552 + 4 * S * x1 * x14 * x185 * x284 * x38 + 4 * S * x1 * x14 * x185 * x344 * x36 + 8 * S * x1 * x18 * x235 * x27 * x284 * x36 + 8 * S * x11 * x18 * x27 * x3 * x344 * x41 * x6 + 8 * S * x11 * x18 * x27 * x3 * x540 * x6 + 8 * S * x11 * x203 * x3 * x552 * x6 + 4 * S * x14 * x18 * x203 * x283 * x9 + 16 * S * x14 * x18 * x282 * x344 * x9 + 4 * S * x14 * x185 * x27 * x362 * x9 + 4 * S * x14 * x185 * x284 * x35 - S * x14 * x189 * x437 + 4 * S * x14 * x192 * x283 + 4 * S * x14 * x203 * x284 * x35 * x9 + 4 * S * x14 * x203 * x284 * x38 * x6 + 4 * S * x14 * x203 * x344 * x36 * x6 + 4 * S * x14 * x269 * x275 * x362 + 16 * S * x14 * x27 * x281 * x284 * x9 + 8 * S * x18 * x203 * x235 * x283 * x9 + 8 * S * x18 * x235 * x27 * x344 * x38 * x6 + 16 * S * x18 * x235 * x281 * x344 * x9 + 8 * S * x185 * x235 * x27 * x362 * x9 - S * x19 * x211 * x518 - S * x22 * x235 * x508 + 16 * S * x235 * x27 * x282 * x284 * x9 - S * x319 * x35 + 32 * Sd ^ 2 * x0 * x9 + 8 * Sd * Sp * x1 * x185 * x3 * x50 + 8 * Sd * Sp * x203 * x3 * x50 * x6 + 8 * Sd * x1 * x18 * x23 * x3 * x35 + 16 * Sd * x11 * x18 * x344 * x50 * x9 + 8 * Sd * x11 * x23 * x281 * x9 + 8 * Sd * x11 * x23 * x282 * x9 + 16 * Sd * x11 * x27 * x284 * x50 * x9 - Sd * x113 * x344 * x401 + 8 * Sd * x18 * x23 * x3 * x38 * x9 - Sd * x18 * x284 * x82 - Sd * x22 * x27 * x65 * x9 + 8 * Sd * x23 * x27 * x3 * x36 * x9 - Sd * x28 * x308 + 4 * Sp * x1 * x11 * x192 * x3 * x344 + 6 * Sp * x1 * x185 * x209 * x27 * x284 + 4 * Sp * x11 * x275 * x284 * x3 * x6 - 8 * Sp * x175 * x192 * x262 + 6 * Sp * x18 * x203 * x209 * x344 * x6 - Sp * x185 * x220 * x307 - Sp * x19 * x273 * x429 - Sp * x190 * x299 - Sp * x209 * x274 * x277 + 8 * kappa * x0 * x1 * x18 * x3 * x35 + 8 * kappa * x0 * x11 * x281 * x9 + 8 * kappa * x0 * x11 * x282 * x9 + 8 * kappa * x0 * x18 * x3 * x38 * x9 + 8 * kappa * x0 * x27 * x3 * x36 * x9 - kappa * x18 * x252 * x37 - kappa * x251 * x5 - kappa * x27 * x405 + 4 * x0 * x1 * x11 * x18 * x209 * x296 + 8 * x0 * x1 * x18 * x295 * x3 + 8 * x0 * x1 * x18 * x296 * x3 + 8 * x0 * x11 * x14 * x27 * x296 * x9 + 24 * x0 * x11 * x18 * x27 * x9 + 4 * x0 * x18 * x209 * x3 * x357 * x9 + 8 * x0 * x224 * x27 * x357 * x6 + 8 * x0 * x27 * x3 * x356 * x6 + 4 * x1 * x11 * x18 * x27 * x290 * x3 * x50 + 4 * x1 * x11 * x18 * x27 * x292 * x3 * x50 + 2 * x1 * x11 * x185 * x224 * x27 * x50 * x517 + 8 * x1 * x11 * x185 * x27 * x3 * x343 * x50 + 4 * x1 * x11 * x185 * x3 * x36 * x38 * x50 + 2 * x1 * x11 * x192 * x224 * x431 * x50 + 2 * x1 * x131 * x192 * x209 * x27 * x50 + 4 * x1 * x131 * x192 * x262 * x27 * x3 * x50 + 2 * x1 * x14 * x185 * x36 * x41 * x50 + 4 * x1 * x14 * x185 * x50 * x555 + 16 * x1 * x18 * x23 * x294 * x3 + 4 * x1 * x18 * x235 * x27 * x50 * x548 + 4 * x1 * x185 * x262 * x27 * x3 * x50 * x517 - x101 * x192 * x362 * x363 - x11 * x123 * x40 * x445 + 6 * x11 * x18 * x203 * x224 * x431 * x50 * x6 + 6 * x11 * x18 * x203 * x3 * x325 * x50 * x9 + 4 * x11 * x18 * x27 * x3 * x347 * x50 * x6 + 4 * x11 * x18 * x27 * x3 * x348 * x50 * x6 + 4 * x11 * x18 * x27 * x3 * x35 * x38 * x50 * x9 + 4 * x11 * x18 * x27 * x3 * x50 * x543 * x6 + x11 * x185 * x209 * x27 * x3 * x431 * x50 * x9 + 6 * x11 * x185 * x27 * x3 * x379 * x50 * x9 + 4 * x11 * x185 * x3 * x35 * x36 * x50 + 4 * x11 * x203 * x3 * x35 * x36 * x50 * x9 + 4 * x11 * x203 * x3 * x36 * x38 * x50 * x6 + 3 * x11 * x209 * x269 * x275 * x3 * x431 * x50 + 2 * x11 * x224 * x275 * x50 * x517 * x6 - x11 * x224 * x342 * x444 - x11 * x225 * x346 * x421 - x11 * x26 * x35 * x365 - x11 * x263 * x312 - x11 * x35 * x364 * x448 - x11 * x510 * x520 - 2 * x115 * x278 * x340 * x38 - x12 * x385 * x386 - 6 * x120 * x28 * x471 * x80 - x120 * x286 * x442 * x489 - x125 * x185 * x240 - x125 * x194 - x125 * x203 * x300 - x130 * x14 * x265 * x413 - x130 * x187 * x268 - x130 * x270 * x330 * x331 - x130 * x271 * x298 + 2 * x131 * x14 * x187 * x235 * x50 + 2 * x131 * x14 * x235 * x269 * x270 * x50 - x131 * x14 * x366 + 2 * x131 * x18 * x209 * x275 * x50 * x6 + 4 * x131 * x18 * x262 * x275 * x3 * x50 * x6 + 4 * x131 * x185 * x188 * x203 * x50 * x9 + 2 * x131 * x187 * x188 * x50 + 2 * x131 * x188 * x269 * x270 * x50 - x131 * x203 * x316 * x371 - x132 * x366 - x133 * x224 * x229 * x506 - 2 * x137 * x28 * x316 + x14 * x155 * x203 * x269 * x50 - x14 * x168 * x196 + 8 * x14 * x18 * x27 * x287 * x50 * x9 + 8 * x14 * x18 * x27 * x291 * x50 * x9 + 2 * x14 * x185 * x232 * x50 + 4 * x14 * x185 * x27 * x431 * x50 * x9 + 4 * x14 * x185 * x286 + 4 * x14 * x185 * x345 * x9 + 2 * x14 * x185 * x50 * x546 + x14 * x185 * x50 * x79 - x14 * x189 * x344 * x410 - x14 * x19 * x33 * x413 * x466 - x14 * x191 * x271 + 2 * x14 * x203 * x269 * x289 * x50 + 4 * x14 * x203 * x269 * x345 + 4 * x14 * x203 * x286 * x9 + 2 * x14 * x203 * x36 * x41 * x50 * x6 + 2 * x14 * x203 * x50 * x546 * x9 + 4 * x14 * x203 * x50 * x555 * x6 + 4 * x14 * x269 * x275 * x431 * x50 - x15 * x434 * x435 + 4 * x155 * x185 * x235 * x50 * x9 + 4 * x155 * x185 * x50 * x9 - x155 * x300 * x301 - x155 * x351 * x353 - x16 * x213 * x35 - x168 * x187 * x211 - x168 * x271 * x293 - x177 * x181 * x59 - x179 * x27 * x295 + 2 * x18 * x188 * x203 * x343 * x50 * x9 - x18 * x202 * x283 - x18 * x220 * x371 * x513 - x18 * x224 * x252 * x296 + 8 * x18 * x235 * x27 * x38 * x41 * x50 * x6 + 4 * x18 * x235 * x27 * x50 * x544 * x6 - x18 * x237 * x238 * x477 - x18 * x281 * x466 * x521 - x18 * x36 * x441 * x5 - x18 * x375 * x431 - x18 * x379 * x516 - x18 * x380 * x381 - x180 * x220 * x222 * x50 - x180 * x222 * x352 - x180 * x433 * x434 - x181 * x33 * x552 - x182 * x197 * x555 - 8 * x182 * x67 - x184 * x186 - x184 * x204 - x185 * x188 * x314 * x40 * x466 + 8 * x185 * x235 * x289 * x50 * x9 - x185 * x238 * x488 - x185 * x24 * x329 * x341 - 12 * x185 * x284 * x512 * x75 - x185 * x360 * x446 - x185 * x371 * x537 - x185 * x376 * x407 - x187 * x188 * x191 + 2 * x188 * x192 * x343 * x50 - x188 * x195 * x196 - x188 * x315 * x344 * x410 - x189 * x192 * x230 * x40 - x189 * x285 - x19 * x235 * x50 * x508 - x19 * x284 * x504 - x19 * x340 * x38 * x427 - x19 * x340 * x516 * x522 - x19 * x381 * x414 - x19 * x506 * x535 - x190 * x192 * x211 * x284 + x192 * x209 * x325 * x50 + x192 * x211 * x343 * x50 - x192 * x332 * x407 - x192 * x367 * x442 * x443 - x192 * x464 * x465 * x75 - x193 * x217 * x417 * x509 - x193 * x340 * x40 * x404 - x194 * x28 * x36 - x194 * x287 - x194 * x291 - x195 * x271 * x293 - x195 * x299 / 2 - x197 * x204 * x231 - x198 * x199 - x199 * x213 * x235 - 4 * x2 * x209 * x439 * x440 - x200 * x266 * x286 - x202 * x430 - x203 * x209 * x338 * x8 + 8 * x203 * x232 * x235 * x50 * x9 + 4 * x203 * x235 * x50 * x79 * x9 - x203 * x236 * x318 - x203 * x333 * x351 * x79 - x203 * x377 * x463 + 4 * x203 * x50 * x79 * x9 - x206 * x211 * x35 - x206 * x331 * x36 - x207 * x35 - x209 * x258 * x386 + x209 * x269 * x275 * x379 * x50 - x209 * x3 * x382 * x418 * x80 - 3 * x209 * x390 * x457 - 3 // 2 * x210 * x36 * x418 * x419 - x211 * x389 * x40 * x463 * x50 - 3 * x211 * x40 * x471 * x66 - x212 * x345 * x346 - x213 * x23 * x284 * x509 * x9 - x214 * x235 - x214 - x215 * x296 * x461 - x216 * x26 - x216 * x36 - x219 * x231 - x22 * x238 * x487 - x22 * x253 * x27 * x449 - 4 * x22 * x284 * x424 * x9 - 8 * x22 * x310 * x401 * x473 - x22 * x35 * x504 - x22 * x38 * x503 * x6 - x22 * x41 * x473 * x507 - x220 * x221 - x220 * x222 * x235 * x371 - x220 * x228 - x220 * x236 * x239 * x41 - x220 * x35 * x406 * x439 - x220 * x468 - x221 * x222 - x222 * x228 - x222 * x235 * x470 - x222 * x238 * x352 - x222 * x27 * x514 - x224 * x243 * x303 - x224 * x365 * x382 - 2 * x224 * x524 * x534 - x225 * x24 * x284 * x413 - 4 * x225 * x373 * x464 * x53 - 2 * x225 * x506 * x524 - x226 * x28 - x226 * x38 - x229 * x240 * x332 - x229 * x485 * x486 + 16 * x23 * x27 * x3 * x354 * x6 - x232 * x425 * x452 - x232 * x447 - x235 * x244 + x235 * x27 * x284 * x344 * x9 * (x45 + x48) - 8 * x235 * x30 * x471 * x69 - x235 * x350 - x236 * x352 * x519 - x237 * x431 * x515 - x239 * x380 * x500 - x24 * x251 * x254 - x240 * x301 * x79 - x241 * x328 * x340 - x242 * x340 - x242 - x244 - x245 * x246 - x245 * x247 - x245 * x323 - x245 * x361 - x246 * x372 - x247 * x427 - x248 * x371 * x512 * x522 - x250 * x310 * x454 - x251 * x309 * x454 - x251 * x368 * x495 - x255 * x273 * x329 - x257 * x258 - x258 * x356 - x259 * x261 - x259 * x297 - x26 * x3 * x414 * x524 - x26 * x400 * x418 - x26 * x41 * x516 - x26 * x416 - x26 * x423 - x261 * x326 - x262 * x464 * x511 - x262 * x522 * x532 * x538 - x264 * x265 - x264 * x327 - x266 * x278 * x343 - x266 * x531 * x537 - x269 * x544 * x549 - x27 * x362 * x378 - x27 * x431 * x462 - x271 * x377 - x272 * x273 * x60 - x272 * x284 * x502 - x276 * x283 * x363 * x54 - x279 * x380 - x280 * x387 * x419 * x6 + 4 * x281 ^ 2 * x50 * x9 - x281 * x428 - x281 * x455 * x456 - x281 * x474 + 4 * x282 ^ 2 * x50 * x9 - x282 * x395 * x521 * x9 - x282 * x428 - x282 * x468 - x282 * x474 - x283 * x307 * x438 - x284 * x38 * x476 * x9 - x284 * x42 * x536 * x54 - x284 * x429 * x465 - x284 * x432 * x482 - x285 * x315 - x287 * x412 - x287 * x458 - x287 * x484 - x288 * x353 * x394 - x289 * x301 * x494 - x289 * x353 * x452 - x29 * x357 * x439 - x29 * x440 * x442 - x290 * x408 - x290 * x497 - x291 * x412 - x291 * x458 - x291 * x484 - x292 * x408 - x292 * x497 - x294 * x305 - x297 * x326 - x30 * x36 * x472 * x69 - x301 * x41 * x514 - x302 * x322 - x304 * x340 * x36 * x490 - x310 * x320 - x310 * x321 - x310 * x422 * x495 - x310 * x494 * x519 - x312 * x313 - x313 * x385 * x41 - x317 * x318 - x319 * x38 * x69 - x320 * x402 - x321 * x402 - x323 * x372 - x323 * x427 - x324 * x41 - x325 * x368 * x438 - x327 * x328 - x330 * x388 * x389 - x334 * x383 - x336 * x337 * x37 - x337 * x403 * x42 - x338 * x339 - x340 * x369 * x38 - x340 * x370 - x340 * x405 * x490 - x340 * x415 - x341 * x373 * x60 - x343 * x381 * x461 - 6 * x343 * x438 * x486 - x344 * x37 * x536 * x54 - x344 * x378 * x40 - 4 * x344 * x42 * x531 * x54 - x345 * x401 * x489 * x493 - x347 * x479 - x347 * x481 - x348 * x479 - x348 * x481 - x349 * x352 - 8 * x35 * x36 * x483 * x66 - x35 * x364 * x432 - x35 * x375 * x40 - x35 * x420 - x35 * x421 * x494 - x350 - x352 * x455 * x500 - x354 * x355 - x358 * x460 - x358 * x8 - x359 * x360 - x359 * x477 - x36 * x388 * x50 * x6 - x36 * x423 - x36 * x448 * x475 - x36 * x466 * x476 - x360 * x493 * x499 - x361 * x372 - x362 * x374 * x467 - x370 - x371 * x382 * x424 - x375 * x376 - x375 * x510 - x38 * x398 * x411 - x38 * x399 * x435 - x381 * x384 - x384 * x456 - x387 * x450 - x388 * x478 - x391 * x393 * x418 - x392 * x393 - x396 * x4 * x441 * x9 - x398 * x400 - x40 * x436 * x456 - 8 * x403 * x417 * x441 - x407 * x533 - x408 * x409 - x409 * x497 - x415 - x416 * x448 - x417 * x511 - 4 * x42 * x466 * x482 - x427 * x530 - x430 * x449 * x525 - x431 * x439 * x520 - x432 * x518 * x54 - x433 * x436 * x66 - x439 * x443 * x444 - x447 * x546 - x450 * x451 - 2 * x456 * x530 - x457 * x501 - x462 * x496 - x469 * x483 * x552 - x470 * x498 - x472 * x488 * x9 - x475 * x546 - x479 * x480 - x479 * x543 - x480 * x481 - x481 * x543 - 2 * x486 * x533 - x487 * x502 * x6 - x491 * x548 - x492 * x544 - x496 * x50 * x515 - x501 * x523 * x85 - x504 * x542 - x525 * x538 * x555 - x526 * x542 - x527 * x528 - x527 * x540 - x528 * x529 - x529 * x540 - x534 * x535 - x548 * x549 * x9
   

    return axx, ayy, axy, bx, by, cc, SS
end
