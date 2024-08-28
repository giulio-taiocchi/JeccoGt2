
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
x8 = 4 * St
x9 = 4 * Sh
x10 = Fy + xi_y
x11 = Fxh + Fyt + 2 * xi_xy
x12 = 4 * Sp
x13 = Fx + xi_x
x14 = Spp * x10
x15 = 2 * x10
x16 = -Fxp
x17 = Gpp * x10
x18 = 2 * x13
x19 = Fyh + xi_yy
x20 = 4 * x10
x21 = x10 ^ 2
x22 = 4 * Spp
x23 = 2 * x19
x24 = 2 * Bh
x25 = 2 * x21
x26 = 2 * Bph
x27 = Fxt + xi_xx
x28 = 4 * x13
x29 = x13 ^ 2
x30 = 2 * x27
x31 = 2 * x29
x32 = 2 * Bt
x33 = 2 * Fx + 2 * xi_x
ABCS[1] = 0
ABCS[2] = -S ^ 3 * u ^ 2 * x1
ABCS[3] = Sp * x1 * x2
ABCS[4] = -12 * S ^ 4 * x0 + S * x0 * (x3 * (-Gh * x8 - Gt * x9) + x5 * (-8 * Sc + 8 * Spt * x10 + x11 * x12 + 8 * x13 * x14)) + S * (Gh * x5 * x9 - x3 * (Bh * x9 + Sph * x20 - 4 * Ss + x12 * x19 - x20 * (Sph + x14) + x21 * x22)) - Sh ^ 2 * x4 + Sh * St * x1 * x5 + x2 * (x3 * (2 * Bh ^ 2 + Bp * x23 + Bpp * x25 - 2 * Bs + Fyp ^ 2 + Fyp * x24 - 2 * Fyph + 2 * Gh ^ 2 + x10 * x26 - x10 * (Bpp * x15 + x26)) + x5 * (-Gp * x23 - Gph * x15 - Gpp * x25 + 2 * Gs - x6 * (Fyp + x24) + (2 * Fy + 2 * xi_y) * (Gph + x17))) + x3 * x7 * (-2 * Gc + Gh * (-Bt - x16) + Gp * x11 + Gpt * x15 + Gt * (Bh + Fyp) + x17 * x18) + x5 * x7 * (-Fxp * Fyp + Fxph + Fypt - Gt * x6) + (S * (Gt * x5 * x8 + x3 * (Bt * x8 + 4 * Sb - Spt * x28 - x12 * x27 - x22 * x29 + x28 * (Spp * x13 + Spt))) - St ^ 2 * x4 + x2 * (x3 * (2 * Bb - Bp * x30 - Bpp * x31 - Bpt * x18 + 2 * Bt ^ 2 + Fxp ^ 2 - Fxp * x32 - 2 * Fxpt + 2 * Gt ^ 2 + x33 * (Bpp * x13 + Bpt)) + x5 * (2 * Gb - Gp * x30 - Gpp * x31 - Gpt * x18 + 2 * Gt * (x16 + x32) + x33 * (Gpp * x13 + Gpt)))) * exp(2 * B)

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
x3 = u ^ 2
x4 = x2 * x3
x5 = tanh(G)
x6 = S ^ 2
x7 = x1 * x6
x8 = Bp * x2
x9 = sech(G)
x10 = Fyp ^ 2
x11 = S * x9
x12 = 4 * x0
x13 = exp(2 * B)
x14 = Fxp * St
x15 = Fxpt * S
x16 = Fxp ^ 2 * S
x17 = sinh(G)
x18 = cosh(G)
x19 = 2 * x0
x20 = 4 * Sh
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x1 * x4
BB[1,2] = 0
CC[1,1] = x7 * (Gp * S * x5 + Sp)
CC[1,2] = x1 * x5 * x8
SS[1] = Bp * Sd * x7 + 8 * Fyp * Sh * x9 + x11 * x12 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp) + x11 * (-4 * Fyph + 2 * x10) - x13 * x9 * (8 * x14 - 4 * x15 + 2 * x16)
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = x19 * x3 * x6 * (-Bp * x17 * x18 + Gp)
BB[2,2] = -x12 * x4
CC[2,1] = -x19 * x8 * sinh(2 * G)
CC[2,2] = Sp * x12 * x6
SS[2] = -Fyp * x17 * x20 + 4 * Gp * Sd * x0 * x6 - 2 * S * x0 * x18 * (Bh * Fxp - Bp * Fxh + Bp * Fyt - Bt * Fyp - Fxp * Fyp + Fxph + Fypt) - S * x17 * (-2 * Fyph + x10) + x0 * x18 * (Fxp * x20 + 4 * Fyp * St) - x13 * x17 * (4 * x14 - 2 * x15 + x16)

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
x19 = Gpp * x10
x20 = 2 * Bph
x21 = Fxh + Fyt + 2 * xi_xy
x22 = Fx + xi_x
x23 = 4 * x22
x24 = 4 * Bt
x25 = Fxt + xi_xx
x26 = x22 ^ 2
x27 = 2 * Fx + 2 * xi_x
ABCS[1] = u ^ 4 * x0 * x2
ABCS[2] = 4 * u ^ 3 * x0 * x1
ABCS[3] = 0
ABCS[4] = S * (-4 * Gh * x6 + x3 * (Sh * x7 + Sph * x11 - 4 * Ss - x11 * (Sph + x14) + x12 * x13 + x8 * x9)) + Sh ^ 2 * x4 + x0 * (4 * S * (x3 * (Gh * St + Gt * Sh) + x5 * (2 * Sc - Sp * x21 - Spt * x17 - 2 * x14 * x22)) - 8 * St * x6 + x15 * (-8 * Sd * Sp + x3 * (-2 * Bh * Gt + 2 * Bt * Gh + 4 * Gc - 2 * Gp * x21 - Gpt * x11 - x19 * x23) + x5 * (-2 * Fxp * Fyp + 4 * Gh * Gt)) + x2 * (Bd * Bp * x3 ^ 2 + Gd * Gp)) + x15 * (x3 * (-2 * Bh ^ 2 - Bp * x16 - Bpp * x18 + 2 * Bs + Fyp ^ 2 - 2 * Gh ^ 2 - x10 * x20 + x10 * (Bpp * x17 + x20)) + x5 * (Gh * x7 + Gp * x16 + Gph * x17 + Gpp * x18 - 2 * Gs - (2 * Fy + 2 * xi_y) * (Gph + x19))) + (S * (-4 * Gt * St * x5 - x3 * (4 * Sb - Spt * x23 + St * x24 - x13 * x26 + x23 * (Spp * x22 + Spt) - x25 * x9)) + St ^ 2 * x4 + x15 * (x3 * (-2 * Bb + 2 * Bp * x25 + 2 * Bpp * x26 + 2 * Bpt * x22 - 2 * Bt ^ 2 + Fxp ^ 2 - 2 * Gt ^ 2 - x27 * (Bpp * x22 + Bpt)) + x5 * (-2 * Gb + 2 * Gp * x25 + 2 * Gpp * x26 + 2 * Gpt * x22 - Gt * x24 - x27 * (Gpp * x22 + Gpt)))) * exp(2 * B)

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
    

   x0 = cosh(G)
x1 = 8 * x0
x2 = S ^ 4
x3 = exp(3 * B)
x4 = x2 * x3
x5 = x1 * x4
x6 = exp(B)
x7 = x2 * x6
x8 = x1 * x7
x9 = exp(2 * B)
x10 = x2 * x9
x11 = sinh(G)
x12 = 16 * x11
x13 = sech(G)
x14 = x0 ^ 2
x15 = 8 * x14
x16 = Fxp * x4
x17 = Fx + xi_x
x18 = Bp * x17
x19 = x15 * x4
x20 = Fy + xi_y
x21 = S ^ 3
x22 = Gp * x17
x23 = Sp * x20
x24 = 16 * x0
x25 = x24 * x9
x26 = x11 * x21
x27 = x25 * x26
x28 = Gp * x20
x29 = Gh + x28
x30 = Bt + x18
x31 = Gt + x22
x32 = x11 * x31
x33 = x15 * x7
x34 = x10 * x11
x35 = Bp * x20
x36 = 16 * x14
x37 = x21 * x6
x38 = x23 * x37
x39 = x10 * x15
x40 = Sp * x17
x41 = Bh + x35
x42 = x11 * x29
x43 = S * x24
x44 = Sph * x43
x45 = Fy * Sp
x46 = Sh * x24
x47 = S * x12
x48 = Gh * x47
x49 = Sp * xi_y
x50 = S ^ 2
x51 = x1 * x50
x52 = x11 * x50
x53 = 8 * x52
x54 = S * x11
x55 = 32 * x6
x56 = x54 * x55
x57 = Sh * x11
x58 = Bp * Fy
x59 = Sh * x43
x60 = Bp * xi_y
x61 = Fyp * x43
x62 = Fy * Gp
x63 = Sh * x47
x64 = Gp * xi_y
x65 = Bp * x51
x66 = Bp * Sd
x67 = x21 * x66
x68 = 24 * Fy
x69 = x0 * x50
x70 = Bph * x69
x71 = 24 * xi_y
x72 = Fypp * x51
x73 = Gph * x52
x74 = Gp * x53
x75 = Gh * x53
x76 = Sd * Sp
x77 = x50 * x6
x78 = 32 * x77
x79 = Fx * Sph
x80 = x47 * x6
x81 = 16 * Fx
x82 = x6 * x81
x83 = Sp * x57
x84 = Spt * x80
x85 = St * x12 * x6
x86 = Gh * x6
x87 = St * x43
x88 = Gt * x6
x89 = Sph * xi_x
x90 = 16 * xi_x
x91 = x6 * x90
x92 = S * x0
x93 = 32 * x92
x94 = 4 * x50
x95 = Fyp ^ 2
x96 = x53 * x6
x97 = x24 * x77
x98 = S * x25
x99 = Fyp * x51
x100 = Fy ^ 2
x101 = Bp * Sp
x102 = x101 * x43
x103 = xi_y ^ 2
x104 = Bpp * Fy
x105 = 32 * xi_y
x106 = 24 * x52
x107 = Fyp * x106
x108 = 24 * Gh * x69
x109 = Fy * Gpp
x110 = Fx * Spp
x111 = x110 * x56
x112 = Fyp * Sp
x113 = x112 * x54
x114 = x6 * x92
x115 = x114 * x81
x116 = Gh * Sp
x117 = Gp * Sh
x118 = Fxp * x80
x119 = x6 * x87
x120 = x43 * x88
x121 = Spp * x56 * xi_x
x122 = x114 * x90
x123 = x24 * x50
x124 = Bpp * x103
x125 = Gpp * x100
x126 = x12 * x50
x127 = Gpp * x103
x128 = x12 * x77
x129 = Bp * x128
x130 = Bp * x52
x131 = x51 * x86
x132 = 8 * Fx
x133 = x52 * x6
x134 = Fypp * x133
x135 = Fx * Gph
x136 = x0 * x77
x137 = 24 * x136
x138 = x81 * x9
x139 = Spt * x92
x140 = Sp * x138
x141 = St * x0
x142 = x51 * x6
x143 = Gp * x142
x144 = Fxp * Fyp
x145 = Fxpp * x96
x146 = Gpt * x142
x147 = 8 * xi_x
x148 = Gph * xi_x
x149 = Gt * x9
x150 = x9 * x90
x151 = Sp * x150
x152 = 32 * x52
x153 = Fx * Sp
x154 = x153 * x56
x155 = Sp * xi_x
x156 = x155 * x56
x157 = x51 * x9
x158 = x53 * x9
x159 = Bp * Gp
x160 = x126 * x159
x161 = Bp ^ 2
x162 = x161 * x50
x163 = x0 * x162
x164 = Fy * xi_y
x165 = Gp ^ 2
x166 = x165 * x69
x167 = 40 * x166
x168 = Fyp * x130
x169 = Bp * Fx
x170 = Gh * x137
x171 = St * x138
x172 = Bp * x92
x173 = Fxp * x128
x174 = x51 * x88
x175 = Bp * xi_x
x176 = St * x150
x177 = Bt * x92
x178 = Bt * x142
x179 = Fxp * x92
x180 = Fx * Gpp
x181 = x180 * x97
x182 = Fx * Gp
x183 = Fyp * x137
x184 = x106 * x86
x185 = Gp * x54
x186 = Gt * x54
x187 = Fxp * x137
x188 = x106 * x88
x189 = x136 * x90
x190 = Gp * xi_x
x191 = 3 * x163
x192 = 20 * x166
x193 = x69 * x9
x194 = 24 * x193
x195 = Bp * x194
x196 = Bpt * x194
x197 = x126 * x9
x198 = Fxpp * x193
x199 = Fx * x9
x200 = Gpt * x106
x201 = x74 * x9
x202 = x9 * xi_x
x203 = x0 * x78
x204 = x169 * x203
x205 = x175 * x9
x206 = x175 * x203
x207 = 3 * G
x208 = x162 * cosh(x207)
x209 = Fxp ^ 2
x210 = x9 * x94
x211 = 2 * Fy
x212 = 40 * x169
x213 = Bt * x193
x214 = 40 * x175
x215 = Fxp * x193
x216 = x149 * x52
x217 = Fx ^ 2
x218 = x101 * x98
x219 = xi_x ^ 2
x220 = Fx * x133
x221 = x161 * x220
x222 = 2 * xi_y
x223 = x133 * xi_x
x224 = x161 * x223
x225 = Bpp * Fx
x226 = x106 * x9
x227 = Bt * x226
x228 = Fxp * x226
x229 = 40 * x165
x230 = x220 * x229
x231 = Gt * x194
x232 = x223 * x229
x233 = Bpp * x50
x234 = Bpp * x219
x235 = Gpp * x219
x236 = x217 * x9
x237 = 48 * x159 * x52
x238 = x219 * x9
x239 = x162 * x6 * sinh(x207)
x240 = Fx * x239
x241 = x199 * xi_x
x242 = x239 * xi_x
x243 = 35 * x163
x244 = x192 * x9
x245 = x208 * x9
x246 = 2 * x0
x247 = 3 * x54
x248 = 8 * S
x249 = S ^ 6 * x9
x250 = S ^ 5
x251 = 4 * x249
x252 = 8 * x34
x253 = 8 * x9
x254 = A * x253
x255 = A * x10
x256 = 4 * x34
x257 = 12 * x0
x258 = x20 ^ 2
x259 = x258 * x7
x260 = x20 ^ 4
x261 = x0 ^ 4
x262 = Bp * S
x263 = x261 * x262
x264 = 4 * Sp
x265 = x20 ^ 3
x266 = 2 * A
x267 = x266 * x34
x268 = 2 * x260
x269 = 4 * x0
x270 = App * x269
x271 = Fyp * x35
x272 = x269 * x7
x273 = A * x272
x274 = A * x269
x275 = x17 ^ 2
x276 = x275 * x4
x277 = 2 * x261
x278 = Bp * x265 * x50
x279 = x277 * x278
x280 = 2 * G
x281 = sinh(x280)
x282 = sinh(4 * G)
x283 = x281 ^ 2
x284 = A * x246
x285 = Bpp * x259 * x284
x286 = x0 ^ 3
x287 = x286 * x35
x288 = x10 * x269
x289 = Ah * x288
x290 = x259 * x269
x291 = Fxp * x20
x292 = Ap * x256
x293 = Fyp * x17
x294 = Ap * x0
x295 = x258 * x264 * x37
x296 = At * x288
x297 = x21 * x3
x298 = x297 * x40
x299 = 2 * S
x300 = x265 * x281
x301 = x300 * x6
x302 = x12 * x21
x303 = Sd * x302
x304 = Gd * Gp
x305 = Fyph + Fypp * x20
x306 = x11 ^ 2
x307 = x20 * x9
x308 = x15 * x50
x309 = x307 * x308
x310 = x306 * x50
x311 = x17 * x253
x312 = x310 * x311
x313 = x17 * x9
x314 = x291 * x293
x315 = x258 * x281
x316 = x11 * x159 * x259 * x266
x317 = A * x0
x318 = x246 * x255
x319 = Fxp * x28
x320 = Fyp * x22
x321 = A * x11
x322 = Gp * x321
x323 = Fyp * x20
x324 = Fyh + x323 + xi_yy
x325 = 8 * Gd
x326 = x11 * x325
x327 = x326 * x7
x328 = Sd * x1
x329 = A * x7
x330 = x246 * x329
x331 = x275 * x297
x332 = x264 * x331
x333 = x259 * x286
x334 = 4 * Bp
x335 = Bd * x334
x336 = x11 ^ 3
x337 = x11 * x94
x338 = x286 * x337
x339 = exp(4 * B)
x340 = x17 ^ 4
x341 = x339 * x340
x342 = x258 * x6
x343 = 8 * x286
x344 = x14 * x265
x345 = 2 * x77
x346 = x344 * x345
x347 = x17 ^ 3
x348 = x3 * x347
x349 = Fxh + x291 + xi_xy
x350 = Fyt + x293 + xi_xy
x351 = Sph + Spp * x20
x352 = Sh + x23
x353 = x317 * x94
x354 = x352 ^ 2
x355 = x354 * x6
x356 = x41 ^ 2
x357 = Fxph + Fxpp * x20
x358 = x29 ^ 2
x359 = Fxpp * x17 + Fxpt
x360 = Fypp * x17 + Fypt
x361 = 2 * x14
x362 = Sdh + Sdp * x20
x363 = Aph + App * x20
x364 = Bdh + Bdp * x20
x365 = A * x161
x366 = x281 * x50
x367 = x159 * x341
x368 = x283 * x341
x369 = x275 * x3
x370 = x17 * x20
x371 = x281 * x370
x372 = x18 * x286
x373 = Bp * x332
x374 = x1 * x10
x375 = x17 * x374
x376 = x302 * x313
x377 = x23 * x376
x378 = S * x264
x379 = 6 * S
x380 = x315 * x379 * x6
x381 = Sd * x17
x382 = Sc + x79 + x89
x383 = Gc + x135 + x148
x384 = Fxp * x17
x385 = Fxt + x384 + xi_xx
x386 = x269 * x4
x387 = x385 * x4
x388 = Bd * x24
x389 = S * x344
x390 = x334 * x352
x391 = x14 * x258
x392 = x248 * x391
x393 = Fyp * x352
x394 = x264 * x389
x395 = x297 * x328
x396 = kappa * x5
x397 = x276 * x286
x398 = Fxp * x29
x399 = Fyp * x31
x400 = x274 * x37
x401 = Gph + Gpp * x20
x402 = x276 * x336
x403 = 8 * Bd * Gp
x404 = x282 / 2
x405 = Bpp * x17
x406 = x348 * x361 * x50
x407 = x20 * x37
x408 = 8 * Sd
x409 = 8 * kappa
x410 = x11 * x403
x411 = cosh(x280)
x412 = x258 * x306
x413 = x37 * x412
x414 = Bpp * x20
x415 = Bph + x414
x416 = St + x40
x417 = x416 ^ 2
x418 = x3 * x417
x419 = x30 ^ 2
x420 = x31 ^ 2
x421 = x349 * x50
x422 = x253 * x291 * x421
x423 = x391 * x50
x424 = x253 * x350
x425 = x50 * x9
x426 = Sdp * x17 + Sdt
x427 = x302 * x307
x428 = Bdp * x17 + Bdt
x429 = x17 * x428
x430 = A * x256
x431 = Bp * x350
x432 = Gp * x349
x433 = Gp * x350
x434 = x3 * x370
x435 = Spp * x17 + Spt
x436 = 4 * x11
x437 = x114 * x265
x438 = x35 * x41
x439 = 2 * x329
x440 = x267 * x370
x441 = x306 * x94
x442 = x441 * x6
x443 = Gp * x265
x444 = x18 * x442 * x443
x445 = x281 * x369 * x379
x446 = x342 * x441
x447 = x345 * x391
x448 = x14 * x267
x449 = Fyp * x18
x450 = x162 * x283
x451 = x274 * x297
x452 = x20 * x415
x453 = x17 * x401
x454 = Gpp * x17 + Gpt
x455 = x20 * x454
x456 = x17 * x4
x457 = x283 * x50
x458 = Bp * x457
x459 = x339 * x347
x460 = x0 * x265
x461 = Bp * x94
x462 = x246 * x52
x463 = x258 * x462
x464 = Fyp * x29
x465 = x21 * x352
x466 = Gd * x25
x467 = x21 * x416
x468 = x41 * x462
x469 = x17 * x30
x470 = x17 * x297
x471 = x32 * x4
x472 = x17 * x471
x473 = Bd * x281
x474 = x1 * x369
x475 = x336 * x94
x476 = 2 * x423
x477 = Bp * x324
x478 = x133 * x258
x479 = x246 * x478
x480 = x18 * x416
x481 = 4 * x287
x482 = x18 * x29
x483 = x318 * x482
x484 = x10 * x317
x485 = x31 * x35
x486 = x22 * x29
x487 = x210 * x321
x488 = x277 * x425
x489 = x275 * x488
x490 = x369 * x441
x491 = x28 * x490
x492 = x14 * x275
x493 = 2 * x492
x494 = x3 * x493 * x50
x495 = x11 * x349
x496 = x114 * x258 * x264
x497 = x11 * x350
x498 = x0 * x475
x499 = x342 * x498
x500 = x352 * x41
x501 = 4 * x42
x502 = x329 * x41
x503 = Bpt + x405
x504 = x339 * x492
x505 = Fxp * x416
x506 = x29 * x348
x507 = x11 * x20
x508 = 4 * x281
x509 = Bd * x29
x510 = x265 * x6
x511 = x337 * x342
x512 = Fxp * x0 * x41
x513 = x0 * x337
x514 = x370 * x513 * x6
x515 = x11 * x281
x516 = x17 * x23
x517 = 8 * x11
x518 = x114 * x517
x519 = 6 * x492
x520 = x459 * x503
x521 = x277 * x50
x522 = Gh + x62 + x64
x523 = x308 * x9
x524 = x20 * x313
x525 = 16 * S
x526 = x306 * x525
x527 = x415 * x52
x528 = x369 * x462
x529 = x17 * x281
x530 = x11 * x4
x531 = x508 * x530
x532 = x39 * x507
x533 = A * x530
x534 = St + x153 + x155
x535 = x286 * x54
x536 = x0 * x431
x537 = x306 * x352
x538 = x525 * x537
x539 = x307 * x416
x540 = x526 * x539
x541 = x248 * x306
x542 = x29 * x40
x543 = 4 * x391
x544 = x264 * x369 * x92
x545 = x337 * x369
x546 = x0 * x507
x547 = x546 * x94
x548 = x274 * x9
x549 = x29 * x31
x550 = x30 * x459
x551 = 2 * x50 * x504
x552 = Fxp * x30
x553 = x476 * x9
x554 = x258 * x92
x555 = x306 * x370
x556 = Bpp * xi_x + Bpt + x225
x557 = 2 * x286
x558 = x265 * x557
x559 = x40 * x416
x560 = x543 * x9
x561 = 4 * x9
x562 = x492 * x561
x563 = x352 * x562
x564 = x546 * (x132 + x147)
x565 = 8 * Bd
x566 = x275 * x339
x567 = x462 * x566
x568 = Fxp * x31
x569 = x463 * x9
x570 = x275 * x9
x571 = x462 * x570
x572 = x17 * x28
x573 = x18 * x30
x574 = x3 * x385
x575 = x425 * x493
x576 = x348 * x557
x577 = x11 * x352 * x416
x578 = Sh + x45 + x49
x579 = Bt + x169 + x175
x580 = x31 * x369
x581 = x18 * x41
x582 = x30 * x35
x583 = x35 * x92
x584 = x369 * x517
x585 = x210 * x42
x586 = x0 * x32
x587 = Bph + Bpp * xi_y + x104
x588 = Bp * x31
x589 = x369 * x51
x590 = S * x36
x591 = Bh + x58 + x60
x592 = Gt + x182 + x190
x593 = 4 * x504
x594 = S * x30 * x416
x595 = x31 * x41
x596 = x42 * x92
x597 = x32 * x416
x598 = x369 * x534
x599 = x18 * x591
x600 = x20 * x369 * x556
x601 = x20 * x311
x602 = x30 * x586
x603 = x210 * x258
x604 = x369 * x579
x605 = 2 * Sph
x606 = Fy * x605 + Spp * x100 + Spp * x103 + Ss + xi_y * (Spp * x211 + x605)
x607 = 2 * Bph
x608 = Bs + Fy * (x104 + x607) + x124 + xi_y * (2 * x104 + x607)
x609 = 2 * Gph
x610 = Fy * x609 + Gs + x125 + x127 + xi_y * (Gpp * x211 + x609)
x611 = S * x606
x612 = 2 * Spt
x613 = Fx * (x110 + x612) + Sb + Spp * x219 + xi_x * (2 * x110 + x612)
x614 = 2 * Bpt
x615 = Bb + Fx * (x225 + x614) + x234 + xi_x * (2 * x225 + x614)
x616 = 2 * Gpt
x617 = Fx * (x180 + x616) + Gb + x235 + xi_x * (2 * x180 + x616)
x618 = S * x613
axx = -x5
ayy = -x8
axy = x10 * x12
bx = x13 * (8 * Fyp * x0 * x11 * x2 * x9 + 8 * Gp * x14 * x2 * x20 * x9 + 16 * Sp * x14 * x17 * x21 * x3 - x11 * x22 * x5 + 8 * x14 * x2 * x29 * x9 - x15 * x16 - x18 * x19 - x19 * x30 - x23 * x27 - x32 * x5)
by = x13 * (Fxp * x1 * x34 - Fyp * x33 - x11 * x28 * x8 + x22 * x39 - x27 * x40 + x31 * x39 + x33 * x35 + x33 * x41 + x36 * x38 - x42 * x8)
cc = -x77 * (Bb * x157 + Bh ^ 2 * x51 + Bh * x248 * (Bp * S * x0 * x20 + Fx * Gp * S * x0 * x6 - Fyp * x92 - 2 * Gh * x54 + Gp * S * x0 * x6 * xi_x + Gt * S * x0 * x6 - Sh * x246 - x246 * x45 - x246 * x49 - x247 * x62 - x247 * x64) - Bpp * x100 * x123 - Bs * x51 + Bt ^ 2 * x157 + Bt * Fxp * x157 + Bt * Gt * x197 + Bt * St * x98 - Bt * x131 + Fx * x196 + 2 * Fx * x245 * xi_x - Fxh * x129 - Fxh * x143 + Fxp * Gt * x158 - Fxp * x131 - Fxph * x96 + Fxpt * x157 + Fxt * x195 + Fxt * x201 + Fy * x111 + Fy * x121 - Fy * x145 - Fy * x146 - Fy * x181 - Fy * x230 - Fy * x232 + Fy * x44 + Fy * x72 + Fy * x84 + Fyh * x65 + Fyh * x74 + Fyp * x75 + Fyph * x51 - Fypt * x96 - Fyt * x129 - Fyt * x143 + Gb * x158 - Gc * x97 + Gh ^ 2 * x51 - Gh * Gt * x128 - Gp * x97 * xi_xy - Gpp * x189 * xi_y + Gpp * x197 * x217 + Gs * x53 + Gt ^ 2 * x157 + Sb * x98 - Sc * x56 - Sh ^ 2 * x24 + Sh * x48 + Ss * x43 - St ^ 2 * x25 + St * x149 * x47 + St * x55 * x57 + x0 * x209 * x210 + x0 * x94 * x95 - x100 * x102 - x100 * x160 + x100 * x191 + x100 * x192 + x100 * x208 - x102 * x103 - x103 * x160 + x103 * x191 + x103 * x192 + x103 * x208 - x104 * x105 * x69 + x105 * x109 * x52 + x107 * x62 + x107 * x64 + x108 * x62 + x108 * x64 - x109 * x189 + x111 * xi_y + x113 * x82 + x113 * x91 - x115 * x116 - x115 * x117 - x116 * x122 - x117 * x122 + x118 * x45 + x118 * x49 - x119 * x62 - x119 * x64 - x120 * x45 - x120 * x49 + x121 * xi_y - x123 * x124 + x125 * x126 + x126 * x127 - x130 * x55 * xi_xy - x132 * x134 + x132 * x198 - x134 * x147 - x135 * x137 - x137 * x148 + x138 * x139 + x139 * x150 - x140 * x141 + x140 * x177 - x140 * x179 + x140 * x186 - x141 * x151 - x144 * x96 - x145 * xi_y - x146 * xi_y + x147 * x198 + x151 * x177 - x151 * x179 + x151 * x186 + x152 * x180 * x202 - x152 * x60 * x62 - x153 * x205 * x93 + x154 * x58 + x154 * x60 + x156 * x58 + x156 * x60 + 6 * x163 * x164 + 70 * x163 * x241 + x164 * x167 + x167 * x241 - x168 * x82 - x168 * x91 - x169 * x170 - x170 * x175 + x171 * x172 + x171 * x185 + x172 * x176 - x173 * x58 - x173 * x60 - x174 * x58 - x174 * x60 + x176 * x185 - x178 * x62 - x178 * x64 - x181 * xi_y - x182 * x183 - x182 * x184 + 96 * x182 * x205 * x52 + x182 * x227 + x182 * x228 + x182 * x231 - x183 * x190 - x184 * x190 - x187 * x62 - x187 * x64 - x188 * x62 - x188 * x64 + x190 * x227 + x190 * x228 + x190 * x231 + 32 * x193 * x225 * xi_x + x195 * xi_xx + x196 * xi_x + x197 * x235 + x199 * x200 + x200 * x202 + x201 * xi_xx - x204 * x62 - x204 * x64 - x206 * x62 - x206 * x64 + x208 * x211 * xi_y - x211 * x221 - x211 * x224 - x211 * x240 - x211 * x242 + x212 * x213 + x212 * x215 + x212 * x216 + x213 * x214 + x214 * x215 + x214 * x216 - x217 * x218 + x217 * x233 * x25 + x217 * x244 + x217 * x245 - x218 * x219 + x219 * x244 + x219 * x245 - x221 * x222 - x222 * x224 - x222 * x240 - x222 * x242 - x230 * xi_y - x232 * xi_y + x234 * x25 * x50 + x236 * x237 + x236 * x243 + x237 * x238 + x238 * x243 + x44 * xi_y - x45 * x46 + x45 * x48 - x45 * x61 + x45 * x85 - x46 * x49 + x48 * x49 - x49 * x58 * x93 - x49 * x61 + x49 * x85 - x55 * x67 - x58 * x59 - x58 * x75 - x58 * x99 - x59 * x60 - x59 * x88 - x60 * x75 - x60 * x99 + x62 * x63 + x63 * x64 + x65 * xi_yy - x68 * x70 + x68 * x73 - 48 * x7 - x70 * x71 + x71 * x73 + x72 * xi_y + x74 * xi_yy - x76 * x78 - x79 * x80 - x80 * x89 + x82 * x83 + x83 * x91 + x84 * xi_y - x86 * x87 - x88 * x99) / 2
SS = 4 * A * Bp * Fxp * x0 * x17 * x2 * x3 + A * Bp * Fxp * x11 * x17 * x2 * x281 * x3 + A * Bp * Fyp * x0 * x17 * x2 * x281 * x9 + 2 * A * Bp * Fyp * x2 * x20 * x286 * x6 + 4 * A * Bp * Gp * x11 * x14 * x2 * x258 * x6 + 2 * A * Bp * Gp * x11 * x2 * x275 * x3 * x411 + 2 * A * Bp * Gp * x11 * x2 * x275 * x3 + 4 * A * Bp * Sp * x0 * x21 * x275 * x3 * x306 + 8 * A * Bp * Sp * x11 * x17 * x20 * x21 * x9 + 4 * A * Bp * Sp * x21 * x258 * x286 * x6 + 2 * A * Bp * x0 * x17 * x2 * x281 * x3 * x31 + A * Bp * x0 * x17 * x2 * x281 * x41 * x9 + 4 * A * Bp * x0 * x17 * x2 * x3 * x30 + A * Bp * x0 * x2 * x20 * x281 * x30 * x9 + 2 * A * Bp * x0 * x2 * x20 * x31 * x411 * x9 + 6 * A * Bp * x0 * x2 * x3 * x385 + 2 * A * Bp * x0 * x2 * x324 * x6 + 2 * A * Bp * x0 * x20 * x21 * x281 * x416 * x9 + 4 * A * Bp * x0 * x20 * x21 * x352 * x6 + 4 * A * Bp * x11 * x14 * x17 * x21 * x352 * x9 + 2 * A * Bp * x11 * x17 * x2 * x281 * x29 * x9 + 4 * A * Bp * x11 * x17 * x2 * x3 * x31 + 2 * A * Bp * x11 * x2 * x20 * x29 * x411 * x6 + 4 * A * Bp * x11 * x2 * x20 * x29 * x6 + 2 * A * Bp * x11 * x20 * x21 * x281 * x352 * x6 + 2 * A * Bp * x17 * x2 * x286 * x3 * x30 + 4 * A * Bp * x17 * x21 * x286 * x3 * x416 + 2 * A * Bp * x2 * x20 * x286 * x41 * x6 + 2 * A * Bpp * x0 * x2 * x275 * x3 + A * Bpp * x11 * x2 * x275 * x281 * x3 + 2 * A * Bpp * x2 * x258 * x286 * x6 + 2 * A * Fxp * Gp * x11 * x17 * x2 * x3 + 4 * A * Fxp * Sp * x11 * x20 * x21 * x9 + 2 * A * Fxp * x0 * x2 * x3 * x30 + 2 * A * Fxp * x11 * x2 * x3 * x31 + 2 * A * Fyp * Gp * x11 * x2 * x20 * x6 + 4 * A * Fyp * Sp * x11 * x17 * x21 * x9 + 2 * A * Fyp * x11 * x2 * x29 * x6 + 8 * A * Gp * Sp * x0 * x17 * x20 * x21 * x9 + 2 * A * Gp * x0 * x17 * x2 * x3 * x31 + 2 * A * Gp * x0 * x2 * x20 * x29 * x6 + 2 * A * Gp * x11 * x17 * x2 * x3 * x30 + 2 * A * Gp * x11 * x2 * x3 * x385 + 2 * A * Gp * x11 * x2 * x324 * x6 + 4 * A * Sp * x0 * x17 * x3 * x416 * x50 + 4 * A * Sp * x0 * x20 * x352 * x50 * x6 + 2 * A * x0 * x161 * x2 * x258 * x6 + 2 * A * x0 * x161 * x2 * x275 * x3 + A * x0 * x165 * x2 * x258 * x6 + A * x0 * x165 * x2 * x275 * x3 + A * x0 * x2 * x20 * x281 * x503 * x9 + A * x0 * x2 * x209 * x3 + 2 * A * x0 * x2 * x3 * x359 + 2 * A * x0 * x2 * x3 * x419 + 2 * A * x0 * x2 * x3 * x420 + 2 * A * x0 * x2 * x3 * x615 + 2 * A * x0 * x2 * x305 * x6 + 2 * A * x0 * x2 * x31 * x41 * x9 + 2 * A * x0 * x2 * x356 * x6 + 2 * A * x0 * x2 * x358 * x6 + A * x0 * x2 * x6 * x95 + 4 * A * x0 * x21 * x3 * x30 * x416 + 4 * A * x0 * x21 * x3 * x613 + 4 * A * x0 * x21 * x6 * x606 + 2 * A * x11 * x14 * x17 * x2 * x415 * x9 + A * x11 * x161 * x2 * x258 * x281 * x6 + A * x11 * x161 * x2 * x275 * x281 * x3 + 2 * A * x11 * x17 * x2 * x3 * x454 + 4 * A * x11 * x17 * x21 * x351 * x9 + A * x11 * x2 * x20 * x281 * x415 * x6 + 2 * A * x11 * x2 * x20 * x401 * x6 + 4 * A * x11 * x2 * x3 * x30 * x31 + 2 * A * x11 * x2 * x3 * x617 + 2 * A * x11 * x2 * x6 * x610 + 4 * A * x11 * x20 * x21 * x435 * x9 + 4 * A * x11 * x21 * x29 * x352 * x6 + 4 * A * x11 * x21 * x3 * x31 * x416 + 8 * A * x11 * x352 * x416 * x50 * x9 - A * x14 * x159 * x276 * x436 - A * x14 * x26 * x35 * x416 * x561 + 2 * A * x17 * x2 * x286 * x3 * x503 - 12 * A * x249 - A * x286 * x373 - A * x352 * x37 * x481 + 4 * Ab * x0 * x2 * x3 - Ac * x252 + 4 * Ah * Bp * x0 * x2 * x20 * x306 * x6 - Ah * Fxp * x256 + 4 * Ah * Fyp * x0 * x2 * x6 + 4 * Ah * Gp * x11 * x2 * x20 * x6 + 8 * Ah * Sp * x11 * x17 * x21 * x9 - Ah * x1 * x38 + 4 * Ah * x11 * x2 * x29 * x6 - Ah * x272 * x41 - 4 * Ah * x287 * x7 + 4 * Ap * Bp * x0 * x2 * x275 * x3 - Ap * Bp * x290 + 4 * Ap * Fxp * x0 * x17 * x2 * x3 + 4 * Ap * Fyp * x0 * x2 * x20 * x6 + 4 * Ap * Gp * x11 * x2 * x258 * x6 + 4 * Ap * Gp * x11 * x2 * x275 * x3 + 8 * Ap * Sd * x250 * x9 + 8 * Ap * Sp * x11 * x17 * x20 * x21 * x9 + 4 * Ap * x11 * x2 * x349 * x9 + 4 * Ap * x11 * x2 * x350 * x9 - Ap * x272 * x324 - Ap * x28 * x375 - Ap * x385 * x386 + 8 * App * x11 * x17 * x2 * x20 * x9 + 4 * As * x0 * x2 * x6 + 4 * At * Bp * x17 * x2 * x286 * x3 + 4 * At * Fxp * x0 * x2 * x3 - At * Fyp * x256 + 4 * At * Gp * x11 * x17 * x2 * x3 + 8 * At * Sp * x11 * x20 * x21 * x9 + 4 * At * x0 * x2 * x3 * x30 - At * x1 * x298 + 4 * At * x11 * x2 * x3 * x31 - At * x18 * x306 * x386 - Bd ^ 2 * x14 * x251 + 8 * Bd * Bp * x11 * x14 * x17 * x2 * x20 * x9 + 8 * Bd * Gp * x0 * x2 * x258 * x281 * x6 + 8 * Bd * Gp * x11 * x14 * x2 * x275 * x3 + 8 * Bd * Gp * x11 * x2 * x275 * x3 * x411 + 8 * Bd * Gp * x2 * x258 * x336 * x6 + 16 * Bd * Sd * x250 * x9 + 8 * Bd * Sp * x0 * x21 * x275 * x3 * x306 + 16 * Bd * Sp * x0 * x21 * x275 * x3 - Bd * Sp * x1 * x413 + 8 * Bd * Sp * x21 * x258 * x286 * x6 - Bd * Sp * x331 * x343 + 8 * Bd * x0 * x17 * x2 * x281 * x3 * x31 + 4 * Bd * x0 * x17 * x2 * x281 * x41 * x9 + 8 * Bd * x0 * x17 * x2 * x29 * x9 + 4 * Bd * x0 * x2 * x20 * x281 * x30 * x9 + 8 * Bd * x0 * x2 * x20 * x31 * x411 * x9 + 8 * Bd * x0 * x2 * x20 * x31 * x9 + 8 * Bd * x0 * x20 * x21 * x281 * x416 * x9 + 16 * Bd * x11 * x14 * x17 * x21 * x352 * x9 + 8 * Bd * x11 * x17 * x2 * x281 * x29 * x9 - Bd * x11 * x17 * x39 * x41 + 8 * Bd * x11 * x2 * x20 * x29 * x411 * x6 + 8 * Bd * x11 * x2 * x349 * x9 + 8 * Bd * x11 * x2 * x350 * x9 + 8 * Bd * x11 * x20 * x21 * x281 * x352 * x6 - Bd * x12 * x31 * x456 - Bd * x14 * x416 * x427 + 8 * Bd * x17 * x2 * x286 * x3 * x30 + 16 * Bd * x17 * x21 * x286 * x3 * x416 + 8 * Bd * x2 * x20 * x286 * x41 * x6 - Bd * x20 * x252 * x281 * x31 - 16 * Bd * x286 * x352 * x407 - Bd * x30 * x532 - Bd * x377 - Bd * x41 * x507 * x508 * x7 - Bd * x469 * x531 + 4 * Bp * Fxp * x0 * x20 * x275 * x3 * x336 * x50 + 2 * Bp * Fxp * x17 * x258 * x261 * x50 * x9 + 2 * Bp * Fxp * x261 * x339 * x347 * x50 + Bp * Fyp * x0 * x11 * x20 * x275 * x281 * x50 * x9 + 4 * Bp * Fyp * x11 * x17 * x258 * x286 * x50 * x6 + Bp * Fyp * x265 * x283 * x50 / 2 + 8 * Bp * Gd * x11 * x14 * x2 * x275 * x3 + 8 * Bp * Gd * x11 * x2 * x258 * x6 + 8 * Bp * Gd * x2 * x258 * x336 * x6 + 4 * Bp * Gp * x11 * x286 * x339 * x340 * x50 + 2 * Bp * Gp * x17 * x265 * x283 * x50 * x6 + 4 * Bp * Gp * x20 * x3 * x306 * x347 * x411 * x50 + 4 * Bp * Gp * x20 * x3 * x306 * x347 * x50 + Bp * Gp * x260 * x281 * x50 + Bp * Gp * x260 * x282 * x50 / 2 + 8 * Bp * S * Sp * x0 * x20 * x3 * x336 * x347 + 8 * Bp * S * Sp * x11 * x17 * x265 * x286 * x6 + 4 * Bp * S * Sp * x14 * x260 + 4 * Bp * S * Sp * x20 * x281 * x3 * x347 + Bp * S * Sp * x260 * x283 + 4 * Bp * S * Sp * x261 * x339 * x340 + 8 * Bp * S * x0 * x11 * x17 * x258 * x352 * x6 + 8 * Bp * S * x0 * x17 * x258 * x336 * x578 * x6 + 4 * Bp * S * x11 * x17 * x258 * x286 * x578 * x6 + 12 * Bp * S * x11 * x20 * x275 * x286 * x3 * x416 + 4 * Bp * S * x11 * x265 * x286 * x416 * x6 + 4 * Bp * S * x11 * x286 * x3 * x347 * x578 + 4 * Bp * S * x14 * x17 * x258 * x416 * x9 + 4 * Bp * S * x14 * x339 * x347 * x416 + Bp * S * x17 * x258 * x283 * x416 * x9 + 4 * Bp * S * x20 * x261 * x275 * x352 * x9 + 4 * Bp * S * x261 * x265 * x352 + Bp * S * x283 * x339 * x347 * x416 + 8 * Bp * Sd * x0 * x21 * x258 * x306 * x6 + 8 * Bp * Sd * x21 * x275 * x286 * x3 + 4 * Bp * x0 * x11 * x17 * x20 * x324 * x50 * x6 + Bp * x0 * x11 * x17 * x258 * x281 * x30 * x50 * x9 + Bp * x0 * x11 * x20 * x275 * x281 * x41 * x50 * x9 + 4 * Bp * x0 * x11 * x258 * x349 * x50 * x6 + 3 * Bp * x0 * x11 * x265 * x281 * x41 * x50 + 4 * Bp * x0 * x11 * x275 * x3 * x349 * x50 + 3 * Bp * x0 * x11 * x281 * x30 * x339 * x347 * x50 + 4 * Bp * x0 * x11 * x31 * x339 * x347 * x50 + 8 * Bp * x0 * x20 * x275 * x336 * x50 * x522 * x9 + 6 * Bp * x11 * x17 * x258 * x286 * x41 * x50 * x6 + 6 * Bp * x11 * x20 * x275 * x286 * x3 * x30 * x50 + 4 * Bp * x11 * x20 * x275 * x286 * x50 * x522 * x9 + 2 * Bp * x11 * x265 * x286 * x30 * x50 * x6 + 4 * Bp * x11 * x265 * x286 * x50 * x522 + 2 * Bp * x11 * x286 * x3 * x347 * x41 * x50 + 8 * Bp * x14 * x17 * x20 * x350 * x50 * x9 + 2 * Bp * x14 * x17 * x258 * x29 * x411 * x50 * x6 + 4 * Bp * x14 * x17 * x258 * x30 * x50 * x9 + 2 * Bp * x14 * x20 * x275 * x3 * x31 * x50 + 4 * Bp * x14 * x20 * x275 * x41 * x50 * x9 + 2 * Bp * x14 * x258 * x385 * x50 * x9 + 2 * Bp * x14 * x265 * x31 * x50 * x6 + 4 * Bp * x14 * x265 * x41 * x50 + 2 * Bp * x14 * x275 * x339 * x385 * x50 + 2 * Bp * x14 * x29 * x3 * x347 * x411 * x50 + 4 * Bp * x14 * x30 * x339 * x347 * x50 + (3 // 2) * Bp * x17 * x258 * x282 * x31 * x50 * x9 + 8 * Bp * x17 * x258 * x29 * x306 * x50 * x6 + 4 * Bp * x17 * x258 * x306 * x411 * x50 * x522 * x6 + 3 * Bp * x20 * x275 * x283 * x3 * x31 * x50 - Bp * x264 * x317 * x413 + Bp * x265 * x283 * x31 * x50 * x6 + 4 * Bp * x265 * x306 * x411 * x50 * x592 * x6 - Bp * x276 * x326 + Bp * x282 * x31 * x339 * x347 * x50 / 2 + 4 * Bp * x29 * x3 * x306 * x347 * x411 * x50 - Bp * x29 * x406 - Bp * x295 * x317 - Bp * x325 * x402 - Bp * x327 * x391 - Bp * x338 * x459 * x592 - Bp * x349 * x430 - Bp * x521 * x550 + 4 * Bpp * x0 * x20 * x3 * x336 * x347 * x50 + 4 * Bpp * x11 * x17 * x265 * x286 * x50 * x6 + 2 * Bpp * x14 * x260 * x50 + 2 * Bpp * x20 * x281 * x3 * x347 * x50 + Bpp * x260 * x283 * x50 / 2 + 2 * Bpp * x261 * x339 * x340 * x50 - Bpp * x266 * x397 + 8 * Fxc * x20 * x306 * x50 * x9 - Fxc * x309 + 8 * Fxp * Fyp * x14 * x17 * x20 * x50 * x9 + 3 * Fxp * Gp * x17 * x258 * x281 * x50 * x9 + Fxp * Gp * x281 * x339 * x347 * x50 - Fxp * Gp * x346 + 4 * Fxp * S * Sp * x14 * x17 * x258 * x9 + 4 * Fxp * S * Sp * x14 * x339 * x347 + 8 * Fxp * S * Sp * x17 * x258 * x306 * x9 + 8 * Fxp * S * x0 * x11 * x258 * x352 * x6 + 8 * Fxp * S * x0 * x11 * x275 * x3 * x352 + 8 * Fxp * S * x17 * x20 * x281 * x3 * x416 + 16 * Fxp * Sd * x0 * x17 * x21 * x3 - Fxp * Sp * x299 * x301 + 4 * Fxp * x0 * x11 * x17 * x20 * x3 * x30 * x50 + 8 * Fxp * x14 * x17 * x20 * x41 * x50 * x9 + 8 * Fxp * x14 * x20 * x349 * x50 * x9 + 6 * Fxp * x14 * x275 * x29 * x3 * x50 + 4 * Fxp * x17 * x20 * x3 * x306 * x31 * x50 + 8 * Fxp * x17 * x306 * x324 * x50 * x9 - Fxp * x18 * x258 * x457 * x9 / 2 + 8 * Fxp * x20 * x350 * x50 * x9 + 4 * Fxp * x258 * x29 * x306 * x50 * x6 - Fxp * x274 * x298 - Fxp * x287 * x545 - Fxp * x458 * x459 / 2 - Fxp * x491 + 8 * Fxs * x14 * x17 * x50 * x9 - Fxs * x312 + 8 * Fyb * x14 * x20 * x50 * x9 - Fyb * x20 * x253 * x310 + 8 * Fyc * x17 * x306 * x50 * x9 - Fyc * x308 * x313 + 3 * Fyp * Gp * x20 * x275 * x281 * x50 * x9 + Fyp * Gp * x265 * x281 * x50 - Fyp * Gp * x406 + 4 * Fyp * S * Sp * x14 * x20 * x275 * x9 + 4 * Fyp * S * Sp * x14 * x265 + 8 * Fyp * S * Sp * x20 * x275 * x306 * x9 + 8 * Fyp * S * x0 * x11 * x258 * x416 * x6 + 8 * Fyp * S * x0 * x11 * x275 * x3 * x416 + 8 * Fyp * S * x17 * x20 * x281 * x352 * x6 + 16 * Fyp * Sd * x0 * x20 * x21 * x6 + 4 * Fyp * x0 * x11 * x258 * x30 * x50 * x6 + 4 * Fyp * x0 * x11 * x275 * x3 * x30 * x50 + 8 * Fyp * x14 * x17 * x350 * x50 * x9 + 6 * Fyp * x14 * x258 * x31 * x50 * x6 + 2 * Fyp * x14 * x258 * x41 * x50 + 2 * Fyp * x14 * x275 * x41 * x50 * x9 + 4 * Fyp * x17 * x20 * x29 * x306 * x50 * x6 + 8 * Fyp * x17 * x349 * x50 * x9 + 8 * Fyp * x20 * x306 * x385 * x50 * x9 - Fyp * x274 * x38 + 4 * Fyp * x275 * x3 * x306 * x31 * x50 - Fyp * x279 - Fyp * x330 * x41 - Gd ^ 2 * x251 + 8 * Gd * Gp * x11 * x17 * x2 * x20 * x9 + 8 * Gd * x0 * x17 * x2 * x41 * x9 + 8 * Gd * x0 * x2 * x349 * x9 + 8 * Gd * x0 * x2 * x350 * x9 + 16 * Gd * x11 * x17 * x21 * x3 * x416 + 16 * Gd * x11 * x20 * x21 * x352 * x6 - Gd * x20 * x30 * x374 + 8 * Gp * Sd * x11 * x21 * x258 * x6 + 8 * Gp * Sd * x11 * x21 * x275 * x3 + 4 * Gp * x0 * x11 * x17 * x20 * x349 * x50 * x9 + 4 * Gp * x0 * x11 * x17 * x20 * x350 * x50 * x9 + 2 * Gp * x14 * x17 * x258 * x31 * x50 * x9 + 2 * Gp * x14 * x20 * x275 * x29 * x50 * x9 + 2 * Gp * x14 * x258 * x349 * x50 * x6 + 2 * Gp * x14 * x265 * x29 * x50 + 2 * Gp * x14 * x275 * x3 * x350 * x50 + 2 * Gp * x14 * x31 * x339 * x347 * x50 + Gp * x17 * x258 * x281 * x30 * x50 * x9 + 4 * Gp * x17 * x258 * x306 * x31 * x50 * x9 + 4 * Gp * x17 * x258 * x306 * x41 * x50 * x6 + 4 * Gp * x20 * x275 * x29 * x306 * x50 * x9 + Gp * x258 * x281 * x324 * x50 + Gp * x258 * x281 * x385 * x50 * x9 + 4 * Gp * x258 * x306 * x350 * x50 * x6 + Gp * x275 * x281 * x324 * x50 * x9 + Gp * x275 * x281 * x339 * x385 * x50 + 4 * Gp * x275 * x3 * x306 * x349 * x50 - Gp * x275 * x473 * x5 + Gp * x281 * x30 * x339 * x347 * x50 - Gp * x29 * x348 * x462 - Gp * x3 * x421 * x519 + 8 * S * Sp * x0 * x11 * x17 * x258 * x41 * x6 + 4 * S * Sp * x0 * x11 * x265 * x29 + 4 * S * Sp * x0 * x11 * x31 * x339 * x347 + 4 * S * Sp * x14 * x17 * x258 * x30 * x9 + 4 * S * Sp * x14 * x258 * x324 + 4 * S * Sp * x14 * x258 * x385 * x9 + 4 * S * Sp * x14 * x275 * x324 * x9 + 4 * S * Sp * x14 * x275 * x339 * x385 + 4 * S * Sp * x14 * x30 * x339 * x347 + 8 * S * Sp * x17 * x20 * x306 * x349 * x9 + 8 * S * Sp * x17 * x20 * x306 * x350 * x9 + 6 * S * Sp * x17 * x258 * x281 * x31 * x9 + 6 * S * Sp * x20 * x275 * x281 * x29 * x9 + 8 * S * x0 * x11 * x17 * x20 * x3 * x30 * x416 + 8 * S * x0 * x11 * x17 * x20 * x3 * x613 + 8 * S * x0 * x11 * x17 * x20 * x6 * x606 + 4 * S * x14 * x17 * x258 * x435 * x9 + 16 * S * x14 * x17 * x350 * x352 * x9 + 4 * S * x14 * x20 * x275 * x351 * x9 + 16 * S * x14 * x20 * x349 * x416 * x9 + 4 * S * x14 * x258 * x29 * x416 * x6 + 4 * S * x14 * x258 * x31 * x352 * x6 + 4 * S * x14 * x258 * x352 * x41 + 4 * S * x14 * x265 * x351 + 4 * S * x14 * x275 * x29 * x3 * x416 + 4 * S * x14 * x275 * x3 * x31 * x352 + 4 * S * x14 * x275 * x352 * x41 * x9 + 4 * S * x14 * x339 * x347 * x435 + 8 * S * x17 * x20 * x29 * x306 * x352 * x6 + 8 * S * x17 * x20 * x3 * x306 * x31 * x416 + 8 * S * x17 * x258 * x306 * x435 * x9 + 16 * S * x17 * x306 * x349 * x352 * x9 + 8 * S * x20 * x275 * x306 * x351 * x9 + 16 * S * x20 * x306 * x350 * x416 * x9 - 4 * S * x23 * x3 * x31 * x492 - S * x23 * x41 * x562 - S * x258 * x261 * x480 * x561 + 4 * S * x258 * x281 * x382 * x6 + 4 * S * x275 * x281 * x3 * x382 - S * x283 * x35 * x352 * x570 - S * x35 * x563 - S * x542 * x543 * x6 + 32 * Sd ^ 2 * x2 * x9 + 8 * Sd * Sp * x0 * x258 * x50 * x6 + 8 * Sd * Sp * x0 * x275 * x3 * x50 + 16 * Sd * kappa * x250 * x9 + 8 * Sd * x0 * x17 * x21 * x29 * x9 + 8 * Sd * x0 * x20 * x21 * x31 * x9 + 8 * Sd * x0 * x20 * x21 * x41 * x6 + 16 * Sd * x11 * x17 * x352 * x50 * x9 + 16 * Sd * x11 * x20 * x416 * x50 * x9 + 8 * Sd * x11 * x21 * x349 * x9 + 8 * Sd * x11 * x21 * x350 * x9 - Sd * x20 * x352 * x97 + 8 * Sp * kappa * x0 * x21 * x258 * x6 + 8 * Sp * kappa * x0 * x21 * x275 * x3 + 4 * Sp * x0 * x11 * x265 * x416 * x6 + 4 * Sp * x0 * x11 * x3 * x347 * x352 + 6 * Sp * x17 * x258 * x281 * x352 * x6 - 8 * Sp * x18 * x336 * x437 + 6 * Sp * x20 * x275 * x281 * x3 * x416 - Sp * x262 * x368 - Sp * x291 * x445 - Sp * x293 * x380 - Sp * x343 * x348 * x35 * x54 + 8 * kappa * x0 * x17 * x2 * x29 * x9 + 8 * kappa * x0 * x2 * x20 * x31 * x9 + 8 * kappa * x0 * x2 * x20 * x41 * x6 + 8 * kappa * x11 * x2 * x349 * x9 + 8 * kappa * x11 * x2 * x350 * x9 - kappa * x324 * x8 - kappa * x377 + 4 * x0 * x11 * x17 * x20 * x3 * x419 * x50 + 4 * x0 * x11 * x17 * x20 * x3 * x420 * x50 + 4 * x0 * x11 * x17 * x20 * x3 * x50 * x615 + 4 * x0 * x11 * x17 * x20 * x31 * x41 * x50 * x9 + 4 * x0 * x11 * x17 * x20 * x356 * x50 * x6 + 4 * x0 * x11 * x17 * x20 * x358 * x50 * x6 + x0 * x11 * x17 * x258 * x281 * x50 * x503 * x9 + 8 * x0 * x11 * x17 * x258 * x415 * x50 * x6 + 6 * x0 * x11 * x17 * x258 * x454 * x50 * x9 + 6 * x0 * x11 * x20 * x275 * x401 * x50 * x9 + 4 * x0 * x11 * x258 * x29 * x31 * x50 * x6 + 4 * x0 * x11 * x258 * x29 * x41 * x50 + 4 * x0 * x11 * x275 * x29 * x3 * x31 * x50 + 4 * x0 * x11 * x275 * x29 * x41 * x50 * x9 + 3 * x0 * x11 * x281 * x339 * x347 * x50 * x503 + 4 * x0 * x161 * x17 * x265 * x336 * x50 * x6 + 4 * x0 * x161 * x20 * x3 * x336 * x347 * x50 + 8 * x0 * x17 * x2 * x3 * (App * x17 + Apt) + 16 * x0 * x17 * x21 * x3 * x426 + 4 * x0 * x17 * x258 * x336 * x50 * x587 * x6 - x0 * x17 * x291 * x585 + 4 * x0 * x2 * x20 * x281 * x428 * x9 + 8 * x0 * x2 * x20 * x363 * x6 + 8 * x0 * x2 * x20 * x364 * x6 + 16 * x0 * x20 * x21 * x362 * x6 - x0 * x275 * x281 * x307 * x527 - 6 * x0 * x28 * x52 * x580 - x0 * x30 * x370 * x585 - 6 * x0 * x478 * x486 - x1 * x313 * x465 * x473 - x1 * x342 * x577 + 8 * x11 * x14 * x17 * x2 * x364 * x9 - x11 * x142 * x258 * x581 + 24 * x11 * x17 * x2 * x20 * x9 + 2 * x11 * x17 * x258 * x286 * x50 * x587 * x6 + 4 * x11 * x2 * x20 * x281 * x364 * x6 + 6 * x11 * x20 * x275 * x286 * x3 * x50 * x503 + 2 * x11 * x265 * x286 * x50 * x503 * x6 - x11 * x28 * x41 * x439 + 2 * x11 * x286 * x3 * x347 * x50 * x587 - x11 * x325 * x387 - 4 * x11 * x348 * x351 * x92 - x11 * x372 * x592 * x603 - x11 * x582 * x589 - x112 * x281 * x299 * x348 - x123 * x3 * x381 * x416 - x130 * x558 * x579 * x6 - x130 * x576 * x591 - x133 * x246 * x31 * x443 - x133 * x556 * x558 + 2 * x14 * x161 * x260 * x306 * x50 + 2 * x14 * x161 * x306 * x339 * x340 * x50 - x14 * x161 * x440 - x14 * x162 * x268 + 8 * x14 * x17 * x20 * x357 * x50 * x9 + 16 * x14 * x17 * x20 * x360 * x50 * x9 + 4 * x14 * x17 * x258 * x50 * x503 * x9 + x14 * x209 * x275 * x339 * x50 + 2 * x14 * x258 * x29 * x30 * x50 * x6 + 2 * x14 * x258 * x305 * x50 + 4 * x14 * x258 * x354 + 2 * x14 * x258 * x359 * x50 * x9 + 4 * x14 * x258 * x383 * x50 * x6 + 4 * x14 * x258 * x417 * x9 + 2 * x14 * x258 * x50 * x608 + x14 * x258 * x50 * x95 - x14 * x262 * x264 * x341 - x14 * x264 * x416 * x459 + 2 * x14 * x275 * x29 * x3 * x30 * x50 + 4 * x14 * x275 * x3 * x383 * x50 + 2 * x14 * x275 * x339 * x359 * x50 + 4 * x14 * x275 * x339 * x417 + 4 * x14 * x275 * x354 * x9 - 3 * x14 * x275 * x425 * x95 + 2 * x14 * x275 * x50 * x608 * x9 + 4 * x14 * x339 * x347 * x50 * x503 - x14 * x378 * x506 - x144 * x267 - x144 * x315 * x77 - x144 * x366 * x369 - x15 * x18 * x307 * x421 - x157 * x18 * x258 * x336 * x592 - x159 * x260 * x338 - x16 * x266 * x372 + 2 * x161 * x17 * x265 * x281 * x50 * x6 - x161 * x17 * x338 * x510 + 2 * x161 * x20 * x281 * x3 * x347 * x50 - x161 * x210 * x275 * x391 + 4 * x161 * x258 * x261 * x275 * x50 * x9 + 2 * x161 * x260 * x261 * x50 + 2 * x161 * x261 * x339 * x340 * x50 - x161 * x286 * x348 * x507 * x94 - x162 * x260 * x283 - x162 * x341 * x361 - x165 * x440 + 8 * x17 * x2 * x286 * x3 * x428 + 8 * x17 * x20 * x3 * x30 * x306 * x31 * x50 + 4 * x17 * x20 * x3 * x306 * x50 * x617 + 4 * x17 * x20 * x306 * x50 * x6 * x610 - x17 * x252 * x363 - 6 * x17 * x286 * x415 * x478 - x17 * x351 * x380 - x17 * x435 * x451 - x17 * x465 * x466 - x18 * x266 * x411 * x471 - x18 * x281 * x284 * x465 * x9 - x18 * x301 * x378 - x18 * x547 * x574 - x197 * x23 * x381 - x20 * x210 * x293 * x586 + 2 * x20 * x261 * x275 * x415 * x50 * x9 - x20 * x281 * x509 * x8 + x20 * x306 * x352 * x416 * x9 * (x81 + x90) - x20 * x312 * x549 - x20 * x343 * x364 * x7 - x20 * x351 * x400 - x20 * x409 * x42 * x7 - x20 * x435 * x445 - x20 * x448 * x503 - x20 * x466 * x467 + 4 * x209 * x258 * x306 * x50 * x9 + 4 * x209 * x258 * x50 * x9 - x209 * x366 * x434 - 3 * x209 * x423 * x9 - x21 * x25 * x28 * x381 - x210 * x306 * x314 - x210 * x357 * x555 - x210 * x371 * x383 - x210 * x452 * x492 - x22 * x289 - x23 * x253 * x275 * x537 - x23 * x30 * x584 * x92 - x23 * x416 * x487 - x23 * x541 * x580 - x23 * x563 - x233 * x261 * x268 - x233 * x341 * x361 - x233 * x368 / 2 - x248 * x393 * x492 * x9 - x248 * x504 * x505 - x250 * x254 * x66 - x253 * x314 * x50 - x253 * x412 * x559 - x254 * x26 * x382 - x255 * x269 * x383 - 8 * x255 * x76 - x257 * x259 - x257 * x276 - 3 * x258 * x283 * x482 * x77 - x258 * x313 * x503 * x521 - x258 * x450 * x570 - x258 * x488 * x573 - x259 * x270 - x259 * x410 * x411 - x260 * x263 * x264 + 2 * x261 * x265 * x415 * x50 - x262 * x265 * x283 * x352 - 4 * x263 * x416 * x459 - x264 * x344 * x352 + x265 * x281 * x401 * x50 + x265 * x283 * x415 * x50 - x265 * x411 * x442 * x588 - x266 * x297 * x480 * x515 - x267 * x28 * x31 - x267 * x281 * x485 - x267 * x357 - x267 * x360 - x267 * x486 - x269 * x276 * x304 - x270 * x276 - x271 * x273 - x271 * x306 * x330 - x271 * x489 - x273 * x438 - 3 // 2 * x275 * x282 * x29 * x35 * x425 + 8 * x275 * x305 * x306 * x50 * x9 + 4 * x275 * x306 * x50 * x9 * x95 - x275 * x352 * x561 * x596 + 4 * x275 * x50 * x9 * x95 - x278 * x283 * x41 - x278 * x29 * x404 - x279 * x41 - x28 * x296 - 2 * x28 * x348 * x458 - x28 * x468 * x570 - x281 * x29 * x330 * x35 + x281 * x339 * x347 * x454 * x50 - 3 * x281 * x460 * x527 - x281 * x533 * x573 - x285 * x306 - x285 - x286 * x439 * x452 - 2 * x287 * x52 * x604 - x288 * x364 * x529 - x289 * x31 - x29 * x296 - x29 * x30 * x318 - 8 * x29 * x41 * x555 * x77 - x29 * x467 * x548 - x290 * x304 - x291 * x292 - x291 * x303 * x9 - x291 * x313 * x538 - x292 * x293 - x293 * x30 * x309 - x293 * x303 * x9 - x293 * x310 * x424 - x293 * x41 * x547 * x6 - x293 * x424 * x50 - x293 * x540 - x294 * x295 - x294 * x332 - x295 * x322 - x30 * x388 * x456 - x30 * x491 - x300 * x345 * x405 - x305 * x425 * x519 - x305 * x514 - x306 * x348 * x411 * x461 * x522 - x306 * x422 - x306 * x474 * x67 - x31 * x394 * x6 - x31 * x465 * x548 - 8 * x310 * x35 * x580 - 12 * x310 * x360 * x524 - x313 * x349 * x352 * x590 - x313 * x350 * x538 - x316 * x411 - x316 - x317 * x373 - x318 * x319 - x318 * x320 - x318 * x398 - x318 * x399 - x318 * x432 - x318 * x433 - x318 * x453 - x318 * x455 - x319 * x494 - x32 * x352 * x601 * x92 - x32 * x408 * x470 - x320 * x446 - x320 * x447 - x322 * x332 - x323 * x385 * x523 - x324 * x327 - x324 * x328 * x37 - x324 * x384 * x523 - x324 * x442 * x572 - x324 * x516 * x518 - x330 * x608 - x333 * x335 - x333 * x365 - x334 * x510 * x534 * x535 - x335 * x397 - 8 * x336 * x583 * x598 - x338 * x348 * x414 - x341 * x450 - x342 * x343 * x67 - 12 * x342 * x352 * x372 * x54 - x342 * x541 * x542 - x344 * x415 * x94 - x346 * x411 * x588 - x346 * x454 - x348 * x390 * x535 + 4 * x349 ^ 2 * x50 * x9 - x349 * x540 - x35 * x411 * x490 * x592 - x35 * x498 * x604 - x35 * x502 * x515 + 4 * x350 ^ 2 * x50 * x9 - x350 * x539 * x590 - x352 * x40 * x487 - x352 * x501 * x554 - x353 * x355 - x353 * x418 - x355 * x564 - x356 * x476 - x356 * x575 - x357 * x479 - x357 * x528 - x358 * x476 - x358 * x575 - x359 * x434 * x513 - x360 * x479 - x360 * x528 - x362 * x376 - x365 * x397 - x366 * x367 - x367 * x404 * x50 - x370 * x500 * x518 - x371 * x77 * x95 - x375 * x411 * x509 - x382 * x524 * x526 - x385 * x395 - x385 * x396 - x387 * x388 - x389 * x390 - x391 * x410 * x7 - 6 * x391 * x433 * x77 - x392 * x393 - x392 * x505 * x9 - x394 * x41 - x395 * x469 - x396 * x469 - x398 * x447 - x398 * x490 - x399 * x446 - x399 * x494 - x400 * x500 - x401 * x406 - x402 * x403 - x405 * x460 * x475 * x6 - x407 * x408 * x42 - x409 * x472 - x411 * x444 - x411 * x472 * x565 - x411 * x483 - x411 * x485 * x494 - x415 * x484 * x529 - x416 * x470 * x515 * x565 - x416 * x583 * x584 - x416 * x596 * x601 - x418 * x564 - x419 * x551 - x419 * x553 - x42 * x460 * x461 - x420 * x551 - x420 * x553 - x421 * x424 - x422 - x426 * x427 - x428 * x532 - x429 * x5 - x429 * x531 - x430 * x431 - x430 * x549 - x432 * x446 - x433 * x490 - x435 * x436 * x437 - x438 * x489 - x441 * x572 * x574 - x443 * x468 - x444 - x446 * x453 - x447 * x453 - x447 * x482 - x447 * x595 - x448 * x449 - x448 * x581 - x448 * x582 - x449 * x499 - x451 * x480 - x455 * x490 - x455 * x494 - x457 * x520 - x458 * x506 - x458 * x550 - x463 * x464 - x463 * x610 - x464 * x571 - x474 * x577 - x476 * x477 - x477 * x575 - x478 * x557 * x599 - x481 * x54 * x598 - x483 - 6 * x484 * x485 - x494 * x595 - x495 * x496 - x495 * x544 - x496 * x497 - x497 * x544 - x498 * x600 - x499 * x599 - x501 * x502 - x503 * x507 * x589 - x503 * x529 * x533 - x511 * x512 - x511 * x536 - x512 * x545 - x514 * x608 - x516 * x517 * x574 * x92 - x52 * x557 * x600 - x520 * x521 - x527 * x576 - x536 * x545 - x543 * x611 - x551 * x552 - x551 * x615 - x552 * x553 - x553 * x615 - x554 * x561 * x597 - x559 * x560 - x560 * x594 - x560 * x618 - x562 * x611 - 4 * x566 * x597 * x92 - x566 * x602 * x94 - x567 * x568 - x567 * x617 - x568 * x569 - x569 * x617 - x571 * x610 - x593 * x594 - x593 * x618 - x602 * x603
   

    return axx, ayy, axy, bx, by, cc, SS
end
