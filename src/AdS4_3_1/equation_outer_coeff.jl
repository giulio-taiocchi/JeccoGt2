
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
    #return esc( :($fs = $f_yy + (Fy + xi_y)  * ( -2*($fp_y) + (Fy + xi_y) * ($fpp) )) )
    return esc( :($fs = $f_yy + (Fy + xi_y)  * ( -2*($fp_y) - (Fy + xi_y) * ($fpp) )) )
end
macro cross_outer(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    #return esc( :($fc = $f_xy  - (Fx + xi_x) * ($fp_y) -
    #             (Fy + xi_y)  * ( $fp_x - (Fx + xi_x)  * ($fpp) ) ) )
    return esc( :($fc = $f_xy  - (Fx + xi_x) * ($fp_y) -
                 (Fy + xi_y)  * ( $fp_x - (Fx + xi_x)  * ($fpp) ) ) )
end


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (u, xi, B, Bp, G, Gp) = vars
	
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
        u, xi, xi_x, xi_y,
        B     ,        G      ,        S      ,
        Bp    ,        Gp     ,        Sp     ,
        Bpp   ,        Gpp    ,        Spp    ,
        B_x   ,        G_x    ,        S_x    ,
        B_y   ,        G_y    ,        S_y    ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,
        Bp_y  ,        Gp_y   ,        Sp_y
    ) = vars

 


x0 = exp(B)
x1 = S ^ 2
x2 = 2 * x1
x3 = u ^ 4 * x2
x4 = cosh(G)
x5 = x4 ^ 2
x6 = Bp * x5
x7 = u ^ 2
x8 = x2 * x7
x9 = x0 * x8
x10 = sinh(G)
x11 = Bp * x10 * x4
x12 = Bp ^ 2
x13 = 2 * Bp
x14 = Bpp * S + Sp * x13
x15 = S * x12 + x14
x16 = 2 * S
x17 = x16 * x5
x18 = 2 * G
x19 = sinh(x18)
x20 = x1 * x19
x21 = Gp * x13
x22 = x20 * x21
x23 = Sp ^ 2
x24 = Gp ^ 2 * x1
x25 = S * Spp
x26 = -4 * x23 + 2 * x24 + 4 * x25
x27 = S * x12 - x14
x28 = S * x21
x29 = cosh(x18)
x30 = Bp * x29
x31 = Gp * x16 * x30
x32 = 4 * Sp
x33 = Gp * x32 + Gpp * x16
x34 = x19 * x27 + x28 - x31 - x33
x35 = S * x0
x36 = 4 * x35
x37 = Gp * x2
x38 = x0 * x37
x39 = Bp * x20
x40 = S_x * x36
x41 = x2 * x5
x42 = Bp_x * x0
x43 = x0 * x2
x44 = x13 * x20
x45 = S * x5
x46 = Gp * x39
x47 = -2 * x23 + x24 + 2 * x25
x48 = x35 * (x15 * x19 - x28 + x31 - x33)
x49 = 4 * S
AA[1,1] = x0 * x3
AA[1,2] = 0
BB[1,1] = x9 * (2 * u - x6)
BB[1,2] = x8 * (Gp + x11)
CC[1,1] = x0 * (x15 * x17 + x22 + x26)
CC[1,2] = S * x34
SS[1] = -B_x * x43 * x6 - B_y * x37 - B_y * x39 + 2 * Bp * G_y * x1 * x29 + 4 * Bp * S * S_y * x10 * x4 + Bp_y * x1 * x19 - G_x * x0 * x44 - G_x * x38 + 4 * Gp * S * S_y + 2 * Gp_y * x1 + S * x34 * xi_y + 4 * S_x * Sp * x0 - Sp_x * x36 + 2 * x0 * xi_x * (x15 * x45 + x46 + x47) - x40 * x6 - x41 * x42
AA[2,1] = 0
AA[2,2] = x3
BB[2,1] = x9 * (Gp - x11)
BB[2,2] = x1 * x7 * (Bp + 4 * u + x30)
CC[2,1] = x48
CC[2,2] = x17 * x27 - x22 + x26
SS[2] = -B_x * x0 * x39 + B_x * x38 - B_y * x2 * x6 + Bp_y * x41 - G_x * x30 * x43 - G_y * x37 + G_y * x44 + Gp * x40 + Gp_x * x43 + S_y * x32 + S_y * x49 * x6 - Sp_y * x49 - x11 * x40 - x20 * x42 + x48 * xi_x + 2 * xi_y * (x27 * x45 - x46 + x47)
    
    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
          u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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
x6 = 2 * Gh
x7 = 2 * x0 * x2
x8 = 4 * Sh
x9 = Fyh + xi_yy
x10 = 4 * Sp
x11 = Fy + xi_y
x12 = 4 * x11
x13 = x11 ^ 2
x14 = 4 * Spp
x15 = Spp * x11
x16 = 4 * St
x17 = Fxh + Fyt + 2 * xi_xy
x18 = Fx + xi_x
x19 = Spp * x18 + Spt
x20 = 2 * x11
x21 = -Fxp
x22 = Gpp * x11
x23 = 2 * x18
x24 = 2 * Fy + 2 * xi_y
x25 = Gpp * x18 + Gpt
x26 = 2 * x9
x27 = 2 * Bh
x28 = 2 * x13
x29 = 2 * Bph
x30 = Fxt + xi_xx
x31 = 4 * x18
x32 = x18 ^ 2
x33 = 2 * x30
x34 = 2 * x32
x35 = 2 * Bt
x36 = 2 * Fx + 2 * xi_x
ABCS[1] = 0
ABCS[2] = -S ^ 3 * u ^ 2 * x1
ABCS[3] = Sp * x1 * x2
ABCS[4] = -12 * S ^ 4 * x0 + S * x0 * (x3 * (-Gh * x16 - Gt * x8) + x5 * (-8 * Sc + 8 * Spt * x11 + x10 * x17 + 8 * x15 * x18 - x19 * (8 * Fy + 8 * xi_y))) + S * (Gh * x5 * x8 - x3 * (Bh * x8 + Sph * x12 - 4 * Ss + x10 * x9 - x12 * (Sph + x15) + x13 * x14)) - Sh ^ 2 * x4 + Sh * St * x1 * x5 + x2 * (x3 * (2 * Bh ^ 2 + Bp * x26 + Bpp * x28 - 2 * Bs + Fyp ^ 2 + Fyp * x27 - 2 * Fyph + 2 * Gh ^ 2 + x11 * x29 - x11 * (Bpp * x20 + x29)) + x5 * (-Gp * x26 - Gph * x20 - Gpp * x28 + 2 * Gs + x24 * (Gph + x22) - x6 * (Fyp + x27))) + x3 * x7 * (-2 * Gc + Gh * (-Bt - x21) + Gp * x17 + Gpt * x20 + Gt * (Bh + Fyp) + x22 * x23 - x24 * x25) + x5 * x7 * (-Fxp * Fyp + Fxph + Fypt - Gt * x6) + (S * (Gt * x16 * x5 + x3 * (Bt * x16 + 4 * Sb - Spt * x31 - x10 * x30 - x14 * x32 + x19 * x31)) - St ^ 2 * x4 + x2 * (x3 * (2 * Bb - Bp * x33 - Bpp * x34 - Bpt * x23 + 2 * Bt ^ 2 + Fxp ^ 2 - Fxp * x35 - 2 * Fxpt + 2 * Gt ^ 2 + x36 * (Bpp * x18 + Bpt)) + x5 * (2 * Gb - Gp * x33 - Gpp * x34 - Gpt * x23 + 2 * Gt * (x21 + x35) + x25 * x36))) * exp(2 * B)

    nothing
end


# this is another coupled equation, for Bd and Gd. the notation used is
#
# ( A11 d_uu Bd + A12 d_uu Gd + B11 d_u Bd + B12 d_u Gd + C11 Bd + C12 Gd ) = -S1
# ( A21 d_uu Bd + A22 d_uu Gd + B21 d_u Bd + B22 d_u Gd + C21 Bd + C22 Gd ) = -S2

function BdGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
       u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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
x8 = sech(G)
x9 = Fyp ^ 2
x10 = S * x8
x11 = 4 * x0
x12 = exp(2 * B)
x13 = Fxp * St
x14 = Fxpt * S
x15 = Fxp ^ 2 * S
x16 = sinh(G)
x17 = 4 * Sh
x18 = cosh(G)
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x1 * x3
BB[1,2] = 0
CC[1,1] = x6 * (Gp * S * x4 + Sp)
CC[1,2] = x1 * x4 * x7
SS[1] = Bp * Sd * x6 + 8 * Fyp * Sh * x8 + x10 * x11 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp) + x10 * (-4 * Fyph + 2 * x9) - x12 * x8 * (8 * x13 - 4 * x14 + 2 * x15)
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x11 * x3
CC[2,1] = -2 * x0 * x7 * sinh(2 * G)
CC[2,2] = Sp * x11 * x5
SS[2] = -Fyp * x16 * x17 + 4 * Gp * Sd * x0 * x5 - 2 * S * x0 * x18 * (Bh * Fxp - Bp * Fxh + Bp * Fyt - Bt * Fyp - Fxp * Fyp + Fxph + Fypt) - S * x16 * (-2 * Fyph + x9) + x0 * x18 * (Fxp * x17 + 4 * Fyp * St) - x12 * x16 * (4 * x13 - 2 * x14 + x15)

    nothing
end


function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
         u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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
    @tilde_outer("Bp")
    @tilde_outer("Gp")
    @tilde_outer("Sp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")
    @hat_outer("Bp")
    @hat_outer("Gp")
    @hat_outer("Sp")
    

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
x6 = Sh * x5
x7 = 4 * Bh
x8 = Fyh + xi_yy
x9 = 4 * Sp
x10 = Fy + xi_y
x11 = 4 * x10
x12 = x10 ^ 2
x13 = 4 * Spp
x14 = Spp * x10
x15 = S ^ 2
x16 = 2 * x8
x17 = 2 * x10
x18 = 2 * x12
x19 = 2 * Fy + 2 * xi_y
x20 = Gpp * x10
x21 = 2 * Bph
x22 = Fxh + Fyt + 2 * xi_xy
x23 = Fx + xi_x
x24 = Spp * x23 + Spt
x25 = 4 * x23
x26 = Gpp * x23 + Gpt
x27 = 4 * Bt
x28 = Fxt + xi_xx
x29 = x23 ^ 2
x30 = 2 * Fx + 2 * xi_x
ABCS[1] = u ^ 4 * x0 * x2
ABCS[2] = 4 * u ^ 3 * x0 * x1
ABCS[3] = 0
ABCS[4] = S * (-4 * Gh * x6 + x3 * (Sh * x7 + Sph * x11 - 4 * Ss - x11 * (Sph + x14) + x12 * x13 + x8 * x9)) + Sh ^ 2 * x4 + x0 * (4 * S * (x3 * (Gh * St + Gt * Sh) + x5 * (2 * Sc - Sp * x22 - Spt * x17 - 2 * x14 * x23 + x19 * x24)) - 8 * St * x6 + x15 * (-8 * Sd * Sp + x3 * (-2 * Bh * Gt + 2 * Bt * Gh + 4 * Gc - 2 * Gp * x22 - Gpt * x11 + 2 * x19 * x26 - x20 * x25) + x5 * (-2 * Fxp * Fyp + 4 * Gh * Gt)) + x2 * (Bd * Bp * x3 ^ 2 + Gd * Gp)) + x15 * (x3 * (-2 * Bh ^ 2 - Bp * x16 - Bpp * x18 + 2 * Bs + Fyp ^ 2 - 2 * Gh ^ 2 - x10 * x21 + x10 * (Bpp * x17 + x21)) + x5 * (Gh * x7 + Gp * x16 + Gph * x17 + Gpp * x18 - 2 * Gs - x19 * (Gph + x20))) + (S * (-4 * Gt * St * x5 - x3 * (4 * Sb - Spt * x25 + St * x27 - x13 * x29 + x24 * x25 - x28 * x9)) + St ^ 2 * x4 + x15 * (x3 * (-2 * Bb + 2 * Bp * x28 + 2 * Bpp * x29 + 2 * Bpt * x23 - 2 * Bt ^ 2 + Fxp ^ 2 - 2 * Gt ^ 2 - x30 * (Bpp * x23 + Bpt)) + x5 * (-2 * Gb + 2 * Gp * x28 + 2 * Gpp * x29 + 2 * Gpt * x23 - Gt * x27 - x26 * x30))) * exp(2 * B)

    nothing
end

function xi_t_eq_coeff(vars::Tuple, ::Outer)
    (
        kappa, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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
    

x0 = exp(2 * B)
x1 = cosh(G)
x2 = S * x1
x3 = x0 * x2
x4 = exp(B)
x5 = sinh(G)
x6 = x4 * x5
x7 = 2 * S
x8 = 2 * Fy
x9 = Sp * x5
x10 = 2 * xi_y
x11 = x2 * x4
x12 = Bt * x11
x13 = Fxp * x11
x14 = S * x6
x15 = 2 * Fx
x16 = Gp * x5
x17 = S * x4
x18 = x16 * x17
x19 = 2 * xi_x
x20 = Fx + xi_x
x21 = Bp * x20
x22 = Bh * x2
x23 = Sp * x1
x24 = Fyp * x2
x25 = Gh * x5
x26 = S * x25
x27 = Fxp * x14
x28 = 2 * Gp
x29 = x28 * x5
x30 = S * x29
x31 = Gt * x11
x32 = Sp * x6
x33 = Fy + xi_y
x34 = Bp * x33
x35 = 2 * x2
x36 = Gp * x11
x37 = 4 * Sp
x38 = Sd * x37
x39 = S ^ 2
x40 = 2 * x39
x41 = Bp * Sd
x42 = x4 * x41
x43 = x33 ^ 2
x44 = x20 ^ 2
x45 = x0 * x1
x46 = x44 * x45
x47 = x15 + x19
x48 = x4 * x47
x49 = x33 * x5
x50 = x1 * x43 + x46 - x48 * x49
x51 = Gp ^ 2
x52 = Bp ^ 2
x53 = x1 ^ 2
x54 = x52 * x53
x55 = x51 + x54
x56 = S / 4
x57 = x55 * x56
x58 = Fyp * x33
x59 = 2 * x1
x60 = x58 * x59
x61 = x16 * x43
x62 = Fxp * x33
x63 = 2 * x6
x64 = Fxp * x20
x65 = x45 * x64
x66 = x0 * x16 * x44
x67 = x20 * x63
x68 = Bp * x1
x69 = 2 * x68
x70 = Fyp * x5
x71 = Gp * x33
x72 = x1 * x71
x73 = Fyh + xi_yy
x74 = x58 + x73
x75 = x1 * x74
x76 = Gh + x71
x77 = x49 * x76
x78 = Fxh + xi_xy
x79 = x62 + x78
x80 = x6 * x79
x81 = Fyt + xi_xy
x82 = Fyp * x20 + x81
x83 = x6 * x82
x84 = Fxt + xi_xx
x85 = x64 + x84
x86 = x45 * x85
x87 = Bh + x34
x88 = x1 * x33
x89 = x1 * x4
x90 = x76 * x89
x91 = x20 * x90
x92 = Gp * x20
x93 = Gt + x92
x94 = x4 * x88
x95 = x93 * x94
x96 = Bt + x21
x97 = x20 * x45 * x96
x98 = x0 * x5
x99 = x93 * x98
x100 = x20 * x99
x101 = x100 + x75 + x77 - x80 - x83 + x86 - x87 * x88 - x91 - x95 + x97
x102 = Fypp * x33
x103 = x1 * x87
x104 = Fxph + Fxpp * x33
x105 = Fypp * x20 + Fypt
x106 = x89 * x93
x107 = Gph + Gpp * x33
x108 = x5 * x71
x109 = x108 * x87
x110 = Gp * x1
x111 = x110 * x4
x112 = Fxpp * x20
x113 = Fxpt + x112
x114 = Gpp * x20 + Gpt
x115 = x76 * x92
x116 = x71 * x93
x117 = 1 / S
x118 = 4 * x1
x119 = Sh ^ 2 * x118
x120 = xi_y ^ 2
x121 = 4 * Spp
x122 = S * x51
x123 = x121 + x122
x124 = S * x54
x125 = x1 / 4
x126 = Sh * x1
x127 = Bh * x126 - Sh * x25 + x23 * x73
x128 = 4 * S
x129 = -Sh * x37
x130 = Spp * x33
x131 = x33 * x55
x132 = x131 * x56
x133 = Bh * Bp
x134 = Bp * Fyp
x135 = 2 * G
x136 = sinh(x135)
x137 = Bp * x136
x138 = Bpt * S
x139 = S * (2 * Gp - x137)
x140 = 4 * St
x141 = cosh(x135)
x142 = Gt * x141
x143 = Bp * x142 * x7 - Gp * x140 - Gpt * x7
x144 = -x128 * (Bp * Sh * x53 - Sph - x130 - x132) + x129 + x17 * (-Bt * x139 + Fxp * x139 + 2 * St * x137 + x136 * x138 + x143) + x40 * (-Bph * x53 + Gh * Gp - Gh * x137 + x133 * x53 + x134 * x53)
x145 = 1 / x39
x146 = x145 * x33
x147 = 2 * Fyph
x148 = 2 * x102
x149 = x147 + x148
x150 = Fyp ^ 2
x151 = 2 * Bs
x152 = 2 * Bh ^ 2
x153 = 2 * Gh ^ 2
x154 = 4 * Bh
x155 = 2 * Gh
x156 = Bh * Fyp * x59 + 2 * Gs * x5 + x1 * x150 - x1 * x151 + x1 * x152 + x1 * x153 - x154 * x25 - x155 * x70 - x29 * x73 + x69 * x73
x157 = St ^ 2 * x118
x158 = -x157
x159 = 4 * Sb
x160 = xi_x ^ 2
x161 = Fx ^ 2
x162 = St * x1
x163 = St * x5
x164 = Gt * x163
x165 = Bt * x162 + x164 - x23 * x84
x166 = x141 * x155
x167 = Bh * x28 + Fyp * x28 - 2 * Gph
x168 = Bp * x53
x169 = x168 * x40
x170 = Spp * x20
x171 = x20 * x57
x172 = x4 * (Bpt * x40 * x53 + Bt * x169 - Fxp * x169 + Gt * x137 * x40 + Gt * x28 * x39 + St * x128 * x168 - St * x37 + x128 * (Spt + x170 + x171))
x173 = -Sh * x7 * (x137 + x28) + x172 + x39 * (Bp * (x136 * (Bh + Fyp) - x166) - Bph * x136 + x167)
x174 = exp(-B)
x175 = x145 * x20
x176 = x174 * x175
x177 = 2 * Fxpt
x178 = 2 * x112 + x177
x179 = Fxp ^ 2
x180 = x1 * x179
x181 = 2 * x5
x182 = Bt ^ 2
x183 = Gt ^ 2
x184 = 4 * Bt * Gt
x185 = Bb * x59 - Bt * Fxp * x59 - Fxp * Gt * x181 + Gb * x181 + x180 + x182 * x59 + x183 * x59 + x184 * x5 - x29 * x84 - x69 * x84
x186 = x146 / 2
x187 = x175 / 2
x188 = Bh * Gt
x189 = Fxp * Gh
x190 = Fyp * Gt
x191 = -2 * Gt * x25 + x38
x192 = -Bt * Gh * x1 - Fxp * Fyp * x5 - 2 * Gc * x1 + x1 * x188 + x1 * x189 + x1 * x190 + x110 * x78 + x110 * x81 + x191
x193 = S ^ 4
x194 = 2 * x33
x195 = S * x55
x196 = x20 / 2
x197 = x195 * x196
x198 = -8 * Sh * x163 + x128 * (Gh * x162 + Gt * x126 + x5 * (-Fy * x197 + 2 * Sc - x170 * x194 - x197 * xi_y) - x78 * x9 - x81 * x9)
x199 = 12 * x193 + x198
x200 = x2 / 2
x201 = Ah / 2
x202 = Sd * x4
x203 = x2 * x33
x204 = A * Sp
x205 = Ah * Sp
x206 = x202 * x40
x207 = S * x49
x208 = S * x70
x209 = x4 / 2
x210 = At * x209
x211 = Gh * x2
x212 = x211 * x33
x213 = S ^ 3
x214 = x209 * x213
x215 = x3 / 2
x216 = Gd * x9
x217 = Ah * x200 * x34
x218 = At * Sp
x219 = x4 * x49
x220 = Gd * x207
x221 = x11 * x20
x222 = x11 * x33
x223 = x1 ^ 3
x224 = Bd * x223
x225 = S * x5
x226 = Ap / 2
x227 = x14 * x226
x228 = At * x215
x229 = S * x98
x230 = At / 2
x231 = S * x223
x232 = x231 * x33
x233 = 2 * Sd
x234 = x20 * x4
x235 = x234 * x70
x236 = A * x32
x237 = Gd * x234
x238 = x117 * x209
x239 = A / 2
x240 = Sp * x239
x241 = Bp * x2
x242 = x17 * x49
x243 = Gd * Gt
x244 = x117 * x4
x245 = x20 * x45
x246 = x20 * x3
x247 = S * x20 * x98
x248 = App * x33
x249 = x14 * x20
x250 = x92 * x98
x251 = Bd * x33
x252 = x141 * x251
x253 = Gd * x247
x254 = A * x200
x255 = x239 * x76
x256 = x0 * x44
x257 = x5 ^ 2
x258 = x21 * x228
x259 = Bd * x136
x260 = Bd * Sp
x261 = 2 * x260
x262 = 2 * Gd
x263 = Gdp * x33
x264 = Aph + x248
x265 = x200 * x33
x266 = S * x239
x267 = Bd * S
x268 = x0 * x20
x269 = 2 * St
x270 = x3 * x44
x271 = Gd * x11
x272 = x231 * x268
x273 = Bd * Bt
x274 = 2 * Bd
x275 = Bd * x5
x276 = 2 * Gd + 2 * x1 * x275
x277 = x14 * x239
x278 = Gd * x76
x279 = S * x16
x280 = Gt * x259
x281 = x239 * (Fy * Gpp + Gph + Gpp * xi_y)
x282 = App * x20 + Apt
x283 = Bp * S
x284 = x239 * x283
x285 = x239 * x36
x286 = Gd * x93
x287 = x239 * x3
x288 = x3 * x85
x289 = x268 * x275
x290 = x196 * x3
x291 = x259 - x262
x292 = x239 * (Fx * Gpp + Gpp * xi_x + Gpt)
x293 = 2 * Sh * x175
x294 = Fxh - Fyt
x295 = x257 * x294 * x4
x296 = x4 * x53
x297 = x294 * x296
x298 = x146 * x269
x299 = A * x21
x300 = Sp * x50
x301 = -x100 + x103 * x33 - x75 - x77 + x80 + x83 - x86 + x91 + x95 - x97
x302 = Fxp * x73
x303 = Fxpp * x43
x304 = x105 * x33
x305 = -Bph + x133 + x134
x306 = 2 * Bp
x307 = S * x28
x308 = -Sh * x128 * x168 + Sph * x128 + x128 * x130 + x129 + x131 * x39 + x17 * (-Bt * x307 + Fxp * x307 + x136 * (Bt * x283 - Fxp * x283 + St * x306 + x138) + x143) + x40 * (Gh * (Gp - x137) + x305 * x53)
x309 = x105 - x187 * x308
x310 = -Sh * x128 * (Gp + x5 * x68) + x172 + x39 * (-Bp * x166 + x136 * x305 + x167)
x311 = x174 * x310
x312 = x104 - x186 * x311
x313 = Fxs - Fyc + Fyp * x81 + Fypp * x20 * x33 + x104 * x33 - x302 - x303 - x304 + x309 * x33 - x312 * x33
x314 = x117 * x20
x315 = Fxp * x78
x316 = Fypp * x44
x317 = x113 * x33
x318 = x105 * x20
x319 = x33 * (x113 - x187 * x311)
x320 = Fxc - Fyb + Fyp * x84 - x112 * x33 + x20 * x309 - x315 + x316 + x317 - x318 - x319
x321 = x117 * x245
x322 = x309 * x5
x323 = x312 * x5
x324 = -Fxs * x5 + Fyc * x5 + Sh * x233 - x102 * x20 * x5 - x104 * x33 * x5 + x302 * x5 + x303 * x5 + x304 * x5 - x322 * x33 + x323 * x33 - x70 * x81
x325 = Fxc * x5 - Fyb * x5 + St * x233 - x112 * x49 + x20 * x322 - x315 * x5 + x316 * x5 + x317 * x5 - x318 * x5 - x319 * x5 + x70 * x84
x326 = 1 / x213
x327 = x174 * x33
x328 = x326 * x327
x329 = 4 * Sh
x330 = Fyt * Gp
x331 = Fxh * Gp
x332 = x39 * x5
x333 = x146 * x308
x334 = S * (-x147 - x148 + x150 + x333)
x335 = Fxp * x140 - S * x177 + S * x179 - x112 * x7 + x311 * x314
x336 = (-Fyp * x329 + x0 * x335 - x334 - x4 * x7 * (x189 - x190 + x330 - x331 + x332 * (Bd * x28 + Bp * x262) + x35 * (x260 + x41))) * sech(G) / 4
x337 = Bdh + Bdp * x33 - x328 * x336
x338 = x136 / 2
x339 = x174 * x20 * x326
x340 = Bdp * x20 + Bdt - x336 * x339
x341 = x338 * x340
x342 = x1 * x40
x343 = x161 * x39
x344 = S * x52
x345 = x119 + x128 * (-x1 * (-Fy * x132 - Spp * x43 + Ss - x132 * xi_y) + x127)
x346 = (-x0 * (-Bb * x342 - Bt * x140 * x2 + Fx * x2 * xi_x * (8 * Spp + 2 * x122 + x141 * x344 + x344) + Fxt * x2 * x37 + Fxt * x28 * x332 + Fxt * x40 * x68 - Gb * x40 * x5 + x1 * x343 * x51 + x121 * x161 * x2 - x128 * x164 + x157 - x159 * x2 + x160 * x2 * (x123 + x124) + x180 * x39 - x182 * x342 - x183 * x342 - x184 * x332 + x223 * x343 * x52 + x7 * xi_xx * (Sp * x59 + x241 + x279)) - x345 + x39 * (x1 * (Fyh * x306 - x150 - x151 + x152 + x153) + x5 * (-Fyh * x28 - Gh * x154 + 2 * Gs) + 2 * xi_yy * (-x16 + x68)) + x4 * (-2 * x193 * (Bd * x168 + Gd * Gp) - x198 + 2 * x39 * (Fxp * x70 + x1 * x28 * xi_xy + x1 * (-Bt * Gh - 2 * Gc + x188 + x330 + x331) + x191))) / x193
x347 = x264 - x327 * x346 / 2
x348 = 2 * x24
x349 = x0 * x335 * x5 / 4 + x329 * x70 / 4 + x334 * x5 / 4 + x4 * (2 * Bd * Bp * x136 * x213 + 2 * Bh * Fxp * S * x1 + 2 * Bp * S * x1 * x81 - Bt * x348 - 4 * Fxp * x126 - Fxp * x348 - Fyp * x1 * x140 - Gd * x37 * x39 - 4 * Gp * Sd * x39 + 2 * S * x1 * x309 + 2 * S * x1 * x312 - x2 * x306 * x78) / 4
x350 = Gdh + x263 - x328 * x349
x351 = Gdp * x20 + Gdt - x339 * x349
x352 = -x0 * (x128 * (x1 * (-Fx * x171 + Sb - Spp * x44 - x171 * xi_x) + x165) + x158 + x39 * (-x1 * (-x176 * x310 + x178) + x185)) / 8 + x345 / 8 - x39 * (-x1 * (x149 - x333) + x156) / 8 + x4 * (x199 + x40 * (-x192 - x322 - x323)) / 8
x353 = x10 + x8
x354 = Sdh + Sdp * x33 - x328 * x352
x355 = x174 * x326 * x352
x356 = Sdp * x20 + Sdt - x20 * x355
axx = -x3
ayy = -x2
axy = x6 * x7
bx = x4 * (2 * Fx * Sp * x1 * x4 + 2 * Fy * Gp * S * x1 + Fyp * S * x5 + Gh * S * x1 + 2 * Gp * S * x1 * xi_y - Gt * x14 + 2 * Sp * x1 * x4 * xi_x - x10 * x9 - 2 * x11 * x21 - x12 - x13 - x15 * x18 - x18 * x19 - x8 * x9)
by = -Fy * x30 + x10 * x23 - x15 * x32 + x15 * x36 - x19 * x32 + x19 * x36 + x22 + x23 * x8 - x24 - x26 + x27 - x30 * xi_y + x31 + x34 * x35
cc = -S * (2 * Bp * x0 * x1 * x20 * x96 + 2 * Bp * x0 * x1 * x85 + 2 * Bp * x0 * x20 * x5 * x93 - Bp * x80 - Bp * x83 + Fxp * x0 * x1 * x96 + Fxp * x0 * x5 * x93 - Fxp * x90 - Fyp * x103 - Fyp * x106 + Fyp * x5 * x76 + Gp * x0 * x1 * x20 * x93 + Gp * x0 * x20 * x5 * x96 + Gp * x0 * x5 * x85 + Gp * x1 * x33 * x76 + Gp * x5 * x74 + x0 * x1 * x113 + x0 * x1 * x20 * (Bpp * x20 + Bpt) + x0 * x114 * x20 * x5 + x1 * (Fyph + x102) - x104 * x6 - x105 * x6 - x106 * x34 - x107 * x20 * x89 + x107 * x33 * x5 - x109 - x111 * x79 - x111 * x82 - x114 * x94 - x115 * x6 - x116 * x6 - x21 * x90 - x88 * (Bph + Bpp * x33)) - Sp * x101 + Sp * (x0 * x44 * x69 - x34 * x67 + x48 * (-x70 - x72) + x60 + x61 - x62 * x63 + 2 * x65 + x66) + x117 * (-x0 * (x128 * (-x125 * (x121 * x160 + x122 * x160 + x123 * x15 * xi_x + x123 * x161 + x124 * x44 - x159) + x165) + x158 + x39 * (-x1 * (-x173 * x176 + x178) + x185)) + x119 + x128 * (x125 * (Fy ^ 2 * x123 - 4 * Ss + x120 * x121 + x120 * x122 + x123 * x8 * xi_y + x124 * x43) + x127) - x39 * (-x1 * (-x144 * x146 + x149) + x156) + x4 * (x199 + x40 * (-x192 - x5 * (x104 - x173 * x174 * x186) - x5 * (x105 - x144 * x187)))) / 4 + x17 * x38 + x40 * x42 - x50 * x57
SS = A * Bp * x288 - A * x117 * x352 + A * x195 * x50 / 8 - A * x209 * x24 * x93 - A * x23 * x58 - A * x39 * x42 + Ab * x215 - Ac * x14 + Ah * S * x108 - Ah * x11 * x92 - Ap * x200 * x73 + Ap * x202 * x39 - Ap * x221 * x71 - App * x200 * x43 - App * x215 * x44 + As * x200 + At * S * x250 - At * x11 * x71 - Bd ^ 2 * x214 * x53 + Bd * Bh * x232 + Bd * Sh * x53 * x67 - Bd * x141 * x211 * x234 + Bd * x206 + Bd * x221 * x76 + Bd * x222 * x93 + Bdh * x203 + Bdp * x2 * x43 - Bdp * x270 - Bdt * x246 - Bh * x220 + Bp * x226 * x270 + Bt * x228 + Bt * x253 + Fxp * x228 + Fxp * x266 * x99 + Fxp * x287 * x96 - Gd ^ 2 * x214 - Gd * x12 * x33 + Gd * x212 - Gd * x225 * x74 - Gd * x229 * x85 - Gdh * x207 + Gdh * x221 - Gdp * x225 * x43 - Gdp * x229 * x44 + Gdt * x222 - Gdt * x247 + Gt * x229 * x230 - S * x142 * x289 + S * x226 * x61 + S * x299 * x99 + Sd ^ 2 * x128 * x4 + Sd * x301 + Sd * x60 - Sh * x194 * x224 + Sh * x276 * x49 - Sp * x20 * x262 * x94 - St * x219 * x274 * x53 + St * x291 * x94 + kappa * (-S * x101 + x206 + x300) + x0 * x21 * x230 * x231 + x0 * x239 * x279 * x85 - x1 * x117 * x313 * x49 + x1 * x353 * x354 - x109 * x266 - x11 * x136 * x196 * x337 - x11 * x21 * x255 - x11 * x239 * x34 * x93 - x115 * x277 - x116 * x277 - x117 * x324 * x88 - x126 * x234 * x276 - x13 * x255 + x16 * x266 * x74 - x163 * x268 * x291 + x20 * x205 * x6 - x20 * x219 * x261 + x20 * x236 * x34 - x201 * x22 - x201 * x231 * x34 + x201 * x24 + x201 * x26 - x201 * x27 - x201 * x31 - x202 * x204 * x7 - x203 * x278 + x204 * x234 * x72 + x204 * x235 - x204 * x256 * x68 - x204 * x65 - x205 * x88 + x207 * x281 + x207 * x337 * x338 + x207 * x350 - x208 * x210 + x208 * x255 - x210 * x211 - x212 * x259 + x216 * x256 + x216 * x43 + x217 * x257 - x217 + x218 * x219 - x218 * x245 + x22 * x237 - x22 * x251 * x257 + x220 * x87 + 2 * x221 * x263 - x221 * x281 - x221 * x350 - x222 * x292 + x222 * x341 - x222 * x351 + x224 * x268 * x269 - x226 * x241 * x43 + x226 * x256 * x279 - x226 * x3 * x84 - x226 * x300 + x227 * x78 + x227 * x81 - x232 * x337 - x233 * x235 - x233 * x6 * x62 + x233 * x65 + x234 * x259 * x26 + x236 * x62 - x237 * x26 + x238 * x78 ^ 2 + x238 * x81 ^ 2 - x239 * x24 * x87 + x239 * x246 * (Bpp * Fx + Bpp * xi_x + Bpt) - x240 * x301 - x240 * x61 - x240 * x66 - x242 * x243 - x242 * x280 - x242 * x282 + x242 * x286 - x242 * x340 * x53 + x243 * x246 - x244 * x320 * x33 * x53 + x244 * x325 * x49 - x244 * x78 * x81 - x246 * x257 * x273 - x246 * x274 * x96 + x246 * x280 - x246 * x286 + x247 * x292 - x247 * x341 + x247 * x351 + x248 * x249 + x249 * x278 + x249 * x337 * x53 - x249 * x347 + x250 * x266 * x96 + x252 * x26 + x252 * x31 - x253 * x96 - x254 * x33 * (Bph + Bpp * Fy + Bpp * xi_y) + x254 * x71 * x76 + x254 * (Fy * Fypp + Fyph + Fypp * xi_y) - x257 * x258 + x258 + x261 * x46 + x264 * x265 + x265 * x347 + x267 * x80 + x267 * x83 + x271 * x79 + x271 * x82 + x272 * x273 + x272 * x340 - x274 * x288 - x277 * (Fxph + Fxpp * Fy + Fxpp * xi_y) - x277 * (Fx * Fypp + Fypp * xi_x + Fypt) + x282 * x290 - x284 * x80 - x284 * x83 - x285 * x79 - x285 * x82 + x287 * x92 * x93 + x287 * (Fx * Fxpp + Fxpp * xi_x + Fxpt) - x289 * x7 * x93 + x290 * (-x174 * x196 * x346 + x282) + x293 * x295 - x293 * x297 - x295 * x298 + x296 * x313 * x314 + x297 * x298 + x299 * x3 * x96 + x314 * x324 * x6 + x320 * x321 * x5 - x321 * x325 - x353 * x356 * x6 - x354 * x48 * x5 + x355 * x50 + x356 * x45 * x47
    return axx, ayy, axy, bx, by, cc, SS
end
