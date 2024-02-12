
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
x10 = Bp * x4 * sinh(G)
x11 = Bp ^ 2
x12 = 2 * Bp
x13 = Bpp * S + Sp * x12
x14 = S * x11 + x13
x15 = 2 * S
x16 = x15 * x5
x17 = 2 * G
x18 = sinh(x17)
x19 = x1 * x18
x20 = Gp * x12
x21 = x19 * x20
x22 = Sp ^ 2
x23 = Gp ^ 2 * x1
x24 = S * Spp
x25 = -4 * x22 + 2 * x23 + 4 * x24
x26 = S * x11 - x13
x27 = S * x20
x28 = cosh(x17)
x29 = Bp * x28
x30 = Gp * x15 * x29
x31 = 4 * Sp
x32 = Gp * x31 + Gpp * x15
x33 = x18 * x26 + x27 - x30 - x32
x34 = S * x0
x35 = 4 * x34
x36 = Gp * x2
x37 = x0 * x36
x38 = Bp * x19
x39 = S_x * x35
x40 = x2 * x5
x41 = Bp_x * x0
x42 = x0 * x2
x43 = x12 * x19
x44 = S * x5
x45 = Gp * x38
x46 = -2 * x22 + x23 + 2 * x24
x47 = x34 * (x14 * x18 - x27 + x30 - x32)
x48 = 4 * S
AA[1,1] = x0 * x3
AA[1,2] = 0
BB[1,1] = x9 * (2 * u - x6)
BB[1,2] = x8 * (Gp + x10)
CC[1,1] = x0 * (x14 * x16 + x21 + x25)
CC[1,2] = S * x33
SS[1] = -B_x * x42 * x6 - B_y * x36 - B_y * x38 + 2 * Bp * G_y * x1 * x28 + 2 * Bp * S * S_y * x18 + Bp_y * x1 * x18 - G_x * x0 * x43 - G_x * x37 + 4 * Gp * S * S_y + 2 * Gp_y * x1 + S * x33 * xi_y + 4 * S_x * Sp * x0 - Sp_x * x35 + 2 * x0 * xi_x * (x14 * x44 + x45 + x46) - x39 * x6 - x40 * x41
AA[2,1] = 0
AA[2,2] = x3
BB[2,1] = x9 * (Gp - x10)
BB[2,2] = x1 * x7 * (Bp + 4 * u + x29)
CC[2,1] = x47
CC[2,2] = x16 * x26 - x21 + x25
SS[2] = -B_x * x0 * x38 + B_x * x37 - B_y * x2 * x6 + Bp_y * x40 - G_x * x29 * x42 - G_y * x36 + G_y * x43 + Gp * x39 + Gp_x * x42 + S_y * x31 + S_y * x48 * x6 - Sp_y * x48 - x10 * x39 - x19 * x41 + x47 * xi_x + 2 * xi_y * (x26 * x44 - x45 + x46)
    
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
x28 = 2 * Bph
x29 = Fxt + xi_xx
x30 = 4 * x18
x31 = x18 ^ 2
x32 = 2 * x29
x33 = 2 * x31
x34 = 2 * Bt
x35 = 2 * Fx + 2 * xi_x
ABCS[1] = 0
ABCS[2] = -S ^ 3 * u ^ 2 * x1
ABCS[3] = Sp * x1 * x2
ABCS[4] = ( -12 * S ^ 4 * x0 + S * x0 * (x3 * (-Gh * x16 - Gt * x8) + x5 * (-8 * Sc + 4 * Sp * x17 + 8 * Spt * x11 - 8 * x15 * x18 - x19 * (8 * Fy + 8 * xi_y))) + S * (Gh * x5 * x8 - x3 * (Bh * x8 + Sph * x12 - 4 * Ss + x10 * x9 - x12 * (Sph + x15) - x13 * x14)) - Sh ^ 2 * x4 + Sh * St * x1 * x5 + x2 * (x3 * (2 * Bh ^ 2 + Bp * x26 - 2 * Bpp * x13 - 2 * Bs + Fyp ^ 2 + Fyp * x27 - 2 * Fyph + 2 * Gh ^ 2 + x11 * x28 - x11 * (Bpp * x20 + x28)) + x5 * (-Gp * x26 - Gph * x20 + 2 * Gpp * x13 + 2 * Gs + x24 * (Gph + x22) - x6 * (Fyp + x27))) + x3 * x7 * (-2 * Gc + Gh * (-Bt - x21) + Gp * x17 + Gpt * x20 + Gt * (Bh + Fyp) - x22 * x23 - x24 * x25) + x5 * x7 * (-Fxp * Fyp + Fxph + Fypt - Gt * x6) + (S * (Gt * x16 * x5 + x3 * (Bt * x16 + 4 * Sb - Spt * x30 - x10 * x29 - x14 * x31 + x19 * x30)) - St ^ 2 * x4 + x2 * (x3 * (2 * Bb - Bp * x32 - Bpp * x33 - Bpt * x23 + 2 * Bt ^ 2 + Fxp ^ 2 - Fxp * x34 - 2 * Fxpt + 2 * Gt ^ 2 + x35 * (Bpp * x18 + Bpt)) + x5 * (2 * Gb - Gp * x32 - Gpp * x33 - Gpt * x23 + 2 * Gt * (x21 + x34) + x25 * x35))) * exp(2 * B))

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
SS[1] = ( Bp * Sd * x6 + 8 * Fyp * Sh * x8 + x10 * x11 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp) + x10 * (-4 * Fyph + 2 * x9) - x12 * x8 * (8 * x13 - 4 * x14 + 2 * x15))
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x11 * x3
CC[2,1] = -2 * x0 * x7 * sinh(2 * G)
CC[2,2] = Sp * x11 * x5
SS[2] = ( -Fyp * x16 * x17 + 4 * Gp * Sd * x0 * x5 - 2 * S * x0 * x18 * (Bh * Fxp - Bp * Fxh + Bp * Fyt - Bt * Fyp - Fxp * Fyp + Fxph + Fypt) - S * x16 * (-2 * Fyph + x9) + x0 * x18 * (Fxp * x17 + 4 * Fyp * St) - x12 * x16 * (4 * x13 - 2 * x14 + x15) )

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
x6 = Sh * x5
x7 = Fy * xi_y
x8 = Fy ^ 2
x9 = 2 * Spp * x3
x10 = xi_y ^ 2
x11 = 4 * S
x12 = S ^ 2
x13 = 8 * Gpp
x14 = 4 * Gpp
x15 = 2 * Bp
x16 = 4 * Bpp
x17 = Gp * S
x18 = 2 * Sp
x19 = 2 * S
x20 = Bt * S
x21 = 2 * Gt
x22 = 2 * x12
x23 = 4 * Spp
x24 = Fx * x23
x25 = xi_x * (Fy + xi_y)
x26 = 2 * Gp
x27 = Fx * x13


ABCS[1] = u ^ 4 * x0 * x2
ABCS[2] = 4 * u ^ 3 * x0 * x1
ABCS[3] = 0
ABCS[4] = (Sh ^ 2 * x4 + x0 * (-8 * St * x6 + x11 * (x3 * (Gh * St + Gt * Sh) + x5 * (-Fxh * Sp + Fy * x24 - Fyt * Sp + 2 * Sc - x18 * xi_xy + x23 * x25 + x24 * xi_y)) + x12 * (-8 * Sd * Sp + x3 * (-Bh * x21 + 2 * Bt * Gh - Fxh * x26 + Fy * x27 - Fyt * x26 + 4 * Gc - 4 * Gp * xi_xy + x13 * x25 + x27 * xi_y) - x5 * (2 * Fxp * Fyp - 4 * Gh * Gt)) + x2 * (Bd * Bp * x3 ^ 2 + Gd * Gp)) + x11 * (Bh * Sh * x3 + Fyh * Sp * x3 - Gh * x6 + Sp * x3 * xi_yy - Spp * x4 * x7 - Ss * x3 - x10 * x9 - x8 * x9) + x12 * (x3 * (-2 * Bh ^ 2 + 8 * Bpp * x7 + 2 * Bs - Fyh * x15 + Fyp ^ 2 - 2 * Gh ^ 2 + x10 * x16 - x15 * xi_yy + x16 * x8) + x5 * (4 * Bh * Gh + 2 * Fyh * Gp + 2 * Gp * xi_yy - 2 * Gs - x10 * x14 - x13 * x7 - x14 * x8)) + (-x19 * x5 * (-Fxt * x17 + Gb * S + x21 * (St + x20)) + x19 * xi_xx * (x17 * x5 + x3 * (Bp * S + x18)) + x3 * (-Bb * x22 + 2 * Bp * Fxt * x12 - Bt ^ 2 * x22 + Fxp ^ 2 * x12 + 4 * Fxt * S * Sp - Gt ^ 2 * x22 - Sb * x11 + 4 * St ^ 2 - 4 * St * x20)) * exp(2 * B))

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
x64 = S * x0
x65 = Spp * x64
x66 = 64 * xi_y
x67 = Fy * x66
x68 = Gp * xi_y
x69 = Bp * x51
x70 = Bp * Sd
x71 = x21 * x70
x72 = 24 * Fy
x73 = x0 * x50
x74 = Bph * x73
x75 = 24 * xi_y
x76 = Fypp * x51
x77 = Gph * x52
x78 = Fy ^ 2
x79 = 32 * x78
x80 = Gp * x53
x81 = Gh * x53
x82 = xi_y ^ 2
x83 = 32 * x82
x84 = Sd * Sp
x85 = x50 * x6
x86 = 32 * x85
x87 = Fx * Sph
x88 = x47 * x6
x89 = Sp * x57
x90 = 16 * Fx
x91 = x6 * x90
x92 = Fy * Spt
x93 = St * x12 * x6
x94 = Gh * x6
x95 = St * x43
x96 = Gt * x6
x97 = 16 * xi_x
x98 = x6 * x97
x99 = x54 * x98
x100 = 32 * x64
x101 = 4 * x50
x102 = Fyp ^ 2
x103 = x53 * x6
x104 = x24 * x85
x105 = S * x25
x106 = Fyp * x51
x107 = Bp * Sp
x108 = x107 * x43
x109 = Bpp * Fy
x110 = 24 * x52
x111 = Fyp * x110
x112 = 24 * Gh * x73
x113 = Gpp * x52
x114 = 64 * Fy
x115 = Fx * Spp
x116 = x54 * x6
x117 = x115 * x116
x118 = Sp * x54
x119 = x6 * x64
x120 = x119 * x90
x121 = Gh * Sp
x122 = Gp * Sh
x123 = Fxp * x88
x124 = x6 * x95
x125 = x43 * x96
x126 = x114 * xi_x
x127 = Spp * x116
x128 = Fyp * Sp
x129 = x119 * x97
x130 = x66 * xi_x
x131 = Bpp * x73
x132 = Bp * x12
x133 = x132 * x85
x134 = Bp * x52
x135 = x51 * x94
x136 = 8 * Fx
x137 = x52 * x6
x138 = Fypp * x137
x139 = Fx * Gph
x140 = x0 * x85
x141 = 24 * x140
x142 = x9 * x90
x143 = Spt * x64
x144 = Sp * x0
x145 = St * x144
x146 = x51 * x6
x147 = Gp * x146
x148 = Fxp * Fyp
x149 = Fxpp * x103
x150 = Fy * Gpt
x151 = 8 * xi_x
x152 = Gt * x12
x153 = Gt * x9
x154 = x9 * x97
x155 = 32 * x52
x156 = Fx * Sp
x157 = x156 * x56
x158 = Sp * xi_x
x159 = x158 * x56
x160 = x51 * x9
x161 = x53 * x9
x162 = Gp * x50
x163 = x132 * x162
x164 = Bp ^ 2
x165 = x164 * x50
x166 = x0 * x165
x167 = Fy * xi_y
x168 = Gp ^ 2
x169 = x168 * x73
x170 = 40 * x169
x171 = Fyp * x134
x172 = Bp * Fx
x173 = Gh * x141
x174 = x142 * x64
x175 = Bp * St
x176 = Fxp * x12 * x85
x177 = x51 * x96
x178 = Bp * xi_x
x179 = x154 * x64
x180 = Bt * Sp
x181 = Bt * x146
x182 = Fxp * Sp
x183 = Fx * Gpp
x184 = x140 * x183
x185 = Fx * Gp
x186 = Fyp * x141
x187 = x110 * x94
x188 = Gp * St * x54
x189 = Gt * x118
x190 = Fxp * x141
x191 = x110 * x96
x192 = Gpp * x140
x193 = Gp * xi_x
x194 = 3 * x166
x195 = 20 * x169
x196 = x73 * x9
x197 = 24 * x196
x198 = Bp * x197
x199 = Bpt * x197
x200 = x50 * x9
x201 = Fxpp * x196
x202 = Fx * x9
x203 = Gpt * x110
x204 = x80 * x9
x205 = x9 * xi_x
x206 = x0 * x86
x207 = x172 * x206
x208 = x178 * x9
x209 = x178 * x206
x210 = 3 * G
x211 = x165 * cosh(x210)
x212 = Fxp ^ 2
x213 = x101 * x9
x214 = 2 * Fy
x215 = 40 * x172
x216 = Bt * x196
x217 = 40 * x178
x218 = Fxp * x196
x219 = x153 * x52
x220 = Fx ^ 2
x221 = x105 * x107
x222 = xi_x ^ 2
x223 = Fx * x137
x224 = x164 * x223
x225 = 2 * xi_y
x226 = x137 * xi_x
x227 = x164 * x226
x228 = Bpp * Fx
x229 = x110 * x9
x230 = Bt * x229
x231 = Fxp * x229
x232 = 40 * x168
x233 = x223 * x232
x234 = Gt * x197
x235 = x226 * x232
x236 = Bpp * x50
x237 = Bpp * x222
x238 = x12 * x200
x239 = Gpp * x222
x240 = x220 * x9
x241 = 48 * Gp * x134
x242 = x222 * x9
x243 = x165 * x6 * sinh(x210)
x244 = Fx * x243
x245 = x202 * xi_x
x246 = x243 * xi_x
x247 = 35 * x166
x248 = x195 * x9
x249 = x211 * x9
x250 = 2 * x0
x251 = 3 * x54
x252 = 8 * S
x253 = S ^ 6 * x9
x254 = S ^ 5
x255 = 4 * x253
x256 = 8 * x34
x257 = 8 * x9
x258 = A * x257
x259 = A * x10
x260 = 4 * x34
x261 = 12 * x0
x262 = x20 ^ 2
x263 = x262 * x7
x264 = x20 ^ 4
x265 = x0 ^ 4
x266 = Bp * x264
x267 = 4 * Sp
x268 = S * x267
x269 = x20 ^ 3
x270 = 2 * A
x271 = x270 * x34
x272 = 2 * x264
x273 = Fyp * x35
x274 = 4 * x0
x275 = x274 * x7
x276 = A * x275
x277 = A * x274
x278 = x17 ^ 2
x279 = x278 * x4
x280 = 2 * x265
x281 = Bp * x269 * x50
x282 = x280 * x281
x283 = 2 * G
x284 = sinh(x283)
x285 = sinh(4 * G)
x286 = x284 ^ 2
x287 = Bpp * x262
x288 = A * x7
x289 = x250 * x288
x290 = x287 * x289
x291 = x0 ^ 3
x292 = x291 * x35
x293 = x10 * x274
x294 = Ah * x293
x295 = x263 * x274
x296 = Fxp * x20
x297 = Ap * x260
x298 = Fyp * x17
x299 = 4 * Ap * x144
x300 = x262 * x37
x301 = At * x293
x302 = x21 * x3
x303 = x302 * x40
x304 = S * x284
x305 = 2 * x304
x306 = x269 * x6
x307 = x21 * x9
x308 = x12 * x307
x309 = Sd * x308
x310 = Gd * Gp
x311 = Fypp * x20
x312 = Fyph + x311
x313 = x11 ^ 2
x314 = x15 * x200
x315 = x257 * x50
x316 = x313 * x315
x317 = x296 * x298
x318 = x284 * x85
x319 = x11 * x263
x320 = Bp * Gp * x270 * x319
x321 = 4 * Bp
x322 = x144 * x321
x323 = A * x300 * x322
x324 = x250 * x259
x325 = Fxp * x28
x326 = Fyp * x22
x327 = A * Gp * x267
x328 = Fyp * x20
x329 = Fyh + x328 + xi_yy
x330 = 8 * Gd
x331 = x11 * x330
x332 = x331 * x7
x333 = Sd * x1
x334 = x278 * x302
x335 = App * x20
x336 = x17 * x256
x337 = x263 * x291
x338 = Bd * x321
x339 = x11 ^ 3
x340 = Bp * x279
x341 = Gp * x101
x342 = x11 * x291
x343 = Bp * S
x344 = exp(4 * B)
x345 = x17 ^ 4
x346 = x344 * x345
x347 = x262 * x6
x348 = 8 * x291
x349 = x14 * x269
x350 = 2 * x85
x351 = x349 * x350
x352 = x17 ^ 3
x353 = x3 * x352
x354 = Fxh + x296 + xi_xy
x355 = Fyt + x298 + xi_xy
x356 = Sph + Spp * x20
x357 = Sh + x23
x358 = x349 * x357
x359 = A * x0
x360 = x101 * x359
x361 = x357 ^ 2
x362 = x361 * x6
x363 = x41 ^ 2
x364 = Fxph + Fxpp * x20
x365 = x29 ^ 2
x366 = Fxpp * x17 + Fxpt
x367 = Fypp * x17 + Fypt
x368 = 2 * x14
x369 = Sdh + Sdp * x20
x370 = Aph + x335
x371 = Bdh + Bdp * x20
x372 = A * x164
x373 = x284 * x50
x374 = Bp * x346
x375 = x286 * x346
x376 = x3 * x373
x377 = x17 * x20
x378 = x18 * x291
x379 = A * x334
x380 = x1 * x10
x381 = x17 * x28
x382 = x12 * x17
x383 = x307 * x382
x384 = x269 * x284
x385 = 6 * x304
x386 = x17 * x23
x387 = x313 * x35
x388 = Fxp * x17
x389 = Fxt + x388 + xi_xx
x390 = x274 * x4
x391 = x389 * x4
x392 = Bd * x24
x393 = 4 * S
x394 = Bp * x393
x395 = x14 * x262
x396 = x252 * x395
x397 = Fyp * x357
x398 = x268 * x349
x399 = x302 * x333
x400 = kappa * x5
x401 = x279 * x291
x402 = Fxp * x29
x403 = Fyp * x31
x404 = x277 * x37
x405 = Gph + Gpp * x20
x406 = x279 * x339
x407 = 8 * Bd * Gp
x408 = x285 / 2
x409 = Bpp * x17
x410 = x353 * x368 * x50
x411 = x20 * x37
x412 = 8 * Sd
x413 = 8 * kappa
x414 = x18 * x313
x415 = cosh(x283)
x416 = x1 * x313
x417 = Bpp * x20
x418 = Bph + x417
x419 = St + x40
x420 = x419 ^ 2
x421 = x3 * x420
x422 = x30 ^ 2
x423 = x31 ^ 2
x424 = x354 * x50
x425 = x257 * x296 * x424
x426 = x395 * x50
x427 = x426 * x9
x428 = x257 * x355
x429 = x298 * x428 * x50
x430 = Sdp * x17 + Sdt
x431 = x20 * x308
x432 = App * x17 + Apt
x433 = Bdp * x17 + Bdt
x434 = x17 * x433
x435 = x20 * x256
x436 = A * x260
x437 = Bp * x355
x438 = Gp * x324
x439 = Spp * x17 + Spt
x440 = 4 * x11
x441 = x35 * x41
x442 = 2 * x288
x443 = x271 * x377
x444 = x306 * x341 * x414
x445 = x278 * x3
x446 = x385 * x445
x447 = x101 * x313
x448 = x347 * x447
x449 = x350 * x395
x450 = x14 * x271
x451 = Fyp * x18
x452 = 8 * x339
x453 = x165 * x286
x454 = x277 * x302
x455 = x20 * x418
x456 = x17 * x405
x457 = Gpp * x17 + Gpt
x458 = x20 * x457
x459 = x17 * x30
x460 = Bd * x31
x461 = x286 * x50
x462 = Bp * x461
x463 = x344 * x352
x464 = x0 * x269
x465 = Bp * x101
x466 = x250 * x52
x467 = x262 * x466
x468 = Fyp * x29
x469 = x14 * x278
x470 = x200 * x469
x471 = x17 * x357
x472 = Gd * x21 * x25
x473 = x20 * x419
x474 = x20 * x30
x475 = Gp * x269
x476 = x41 * x466
x477 = x3 * x419
x478 = x32 * x4
x479 = x17 * x478
x480 = Bd * x284
x481 = 2 * x426
x482 = Bp * x329
x483 = x137 * x262
x484 = x250 * x483
x485 = x18 * x419
x486 = 4 * x292
x487 = x18 * x29
x488 = x324 * x487
x489 = x10 * x359
x490 = x31 * x35
x491 = x22 * x29
x492 = x11 * x213
x493 = x200 * x280
x494 = x278 * x493
x495 = x445 * x447
x496 = x28 * x495
x497 = x3 * x469
x498 = 2 * x497 * x50
x499 = x11 * x354
x500 = x119 * x262 * x267
x501 = x11 * x355
x502 = x0 * x101
x503 = x339 * x502
x504 = x347 * x503
x505 = x357 * x41
x506 = 4 * x42
x507 = x288 * x41
x508 = Bpt + x409
x509 = x344 * x469
x510 = Fxp * x419
x511 = x469 * x9
x512 = x29 * x353
x513 = x17 * x347
x514 = x11 * x20
x515 = 4 * x284
x516 = x101 * x353
x517 = Fxp * x41
x518 = x11 * x502
x519 = x11 * x377
x520 = x502 * x519
x521 = x520 * x6
x522 = x284 * x35
x523 = 8 * x11
x524 = x386 * x523
x525 = x463 * x508
x526 = x280 * x50
x527 = Gh + x62 + x68
x528 = x313 * x354
x529 = 6 * Gp
x530 = x418 * x52
x531 = x445 * x466
x532 = x17 * x284
x533 = x11 * x4
x534 = x515 * x533
x535 = x39 * x514
x536 = A * x533
x537 = x262 * x9
x538 = St + x156 + x158
x539 = x291 * x321
x540 = x437 * x518
x541 = x357 * x9
x542 = 16 * S
x543 = x313 * x541 * x542
x544 = x17 * x296
x545 = x473 * x9
x546 = x542 * x545
x547 = x252 * x313
x548 = x29 * x40
x549 = 4 * x395
x550 = x267 * x445 * x64
x551 = x292 * x445
x552 = x0 * x514
x553 = x101 * x552
x554 = x277 * x307
x555 = x29 * x31
x556 = x30 * x463
x557 = 2 * x50 * x509
x558 = Fxp * x30
x559 = x481 * x9
x560 = x262 * x64
x561 = x313 * x377
x562 = x213 * x561
x563 = Bpp * xi_x + Bpt + x228
x564 = 2 * x291
x565 = x563 * x564
x566 = x257 * x313
x567 = x40 * x419
x568 = x549 * x9
x569 = x278 * x357
x570 = 4 * x9
x571 = x469 * x570
x572 = x552 * (x136 + x151)
x573 = x302 * x419
x574 = x14 * x419
x575 = x278 * x344
x576 = x466 * x575
x577 = Fxp * x31
x578 = x467 * x9
x579 = x278 * x9
x580 = x466 * x579
x581 = x447 * x6
x582 = x18 * x30
x583 = 2 * x470
x584 = x353 * x564
x585 = x1 * x11
x586 = Sh + x45 + x49
x587 = Bt + x172 + x178
x588 = x31 * x445
x589 = x18 * x270
x590 = x30 * x35
x591 = x35 * x64
x592 = x0 * x213
x593 = x42 * x592
x594 = x262 * x493
x595 = Bph + Bpp * xi_y + x109
x596 = Bp * x31
x597 = x445 * x51
x598 = S * x36
x599 = Bh + x58 + x60
x600 = Gt + x185 + x193
x601 = x0 * x52
x602 = x393 * x509
x603 = x30 * x419
x604 = S * x568
x605 = x31 * x41
x606 = x32 * x419
x607 = x445 * x538
x608 = x18 * x599
x609 = x20 * x445
x610 = x257 * x377 * x64
x611 = x30 * x32
x612 = x213 * x262
x613 = x214 + x225
x614 = Spp * x262 + Ss + x356 * x613
x615 = Bs + x287 + x418 * x613
x616 = Gpp * x262 + Gs + x405 * x613
x617 = S * x614
x618 = 2 * Spt
x619 = Fx * (x115 + x618) + Sb + Spp * x222 + xi_x * (2 * x115 + x618)
x620 = 2 * Bpt
x621 = Bb + Fx * (x228 + x620) + x237 + xi_x * (2 * x228 + x620)
x622 = 2 * Gpt
x623 = Fx * (x183 + x622) + Gb + x239 + xi_x * (2 * x183 + x622)
x624 = 3 * x115
x625 = 3 * Spp
x626 = Fy * x624 + Sc + x87 + x92 + xi_x * (Fy * x625 + Sph + x625 * xi_y) + xi_y * (Spt + x624)
x627 = 3 * x183
x628 = 3 * Gpp
x629 = Fy * x627 + Gc + x139 + x150 + xi_x * (Fy * x628 + Gph + x628 * xi_y) + xi_y * (Gpt + x627)
axx = -x5
ayy = -x8
axy = x10 * x12
bx = x13 * (8 * Fyp * x0 * x11 * x2 * x9 + 8 * Gp * x14 * x2 * x20 * x9 + 16 * Sp * x14 * x17 * x21 * x3 - x11 * x22 * x5 + 8 * x14 * x2 * x29 * x9 - x15 * x16 - x18 * x19 - x19 * x30 - x23 * x27 - x32 * x5)
by = x13 * (Fxp * x1 * x34 - Fyp * x33 - x11 * x28 * x8 + x22 * x39 - x27 * x40 + x31 * x39 + x33 * x35 + x33 * x41 + x36 * x38 - x42 * x8)
cc = -x85 * (Bb * x160 + Bh ^ 2 * x51 + Bh * x252 * (Bp * S * x0 * x20 + Fx * Gp * S * x0 * x6 - Fyp * x64 - 2 * Gh * x54 + Gp * S * x0 * x6 * xi_x + Gt * S * x0 * x6 - Sh * x250 - x144 * x214 - x144 * x225 - x251 * x62 - x251 * x68) - Bs * x51 + Bt ^ 2 * x160 + Bt * Fxp * x160 + Bt * St * x105 - Bt * x135 + Bt * x152 * x200 + Fx * x199 + 2 * Fx * x249 * xi_x - Fxh * x133 - Fxh * x147 + Fxp * Gt * x161 - Fxp * x135 - Fxph * x103 + Fxpt * x160 + Fxt * x198 + Fxt * x204 - Fy * x149 - Fy * x233 - Fy * x235 + Fy * x44 + Fy * x76 + Fyh * x69 + Fyh * x80 + Fyp * x118 * x91 + Fyp * x81 + Fyph * x51 - Fypt * x103 - Fyt * x133 - Fyt * x147 + Gb * x161 - Gc * x104 + Gh ^ 2 * x51 - Gh * x152 * x85 - Gp * x104 * xi_xy - Gph * x141 * xi_x + Gpp * x220 * x238 - Gpt * x140 * x75 + Gs * x53 + Gt ^ 2 * x160 + Sb * x105 - Sc * x56 - Sh ^ 2 * x24 + Sh * x48 - Sph * x99 - Spt * x88 * xi_y + Ss * x43 - St ^ 2 * x25 + St * x153 * x47 + St * x55 * x57 + x0 * x101 * x102 + x0 * x212 * x213 - x100 * x156 * x208 - x100 * x49 * x58 - x103 * x148 - x106 * x58 - x106 * x60 - x106 * x96 - x108 * x78 - x108 * x82 - x109 * x66 * x73 + x111 * x62 + x111 * x68 + x112 * x62 + x112 * x68 + x113 * x67 + x113 * x79 + x113 * x83 - x114 * x117 - x114 * x184 - x117 * x66 - x120 * x121 - x120 * x122 - x121 * x129 - x122 * x129 + x123 * x45 + x123 * x49 - x124 * x62 - x124 * x68 - x125 * x45 - x125 * x49 - x126 * x127 - x126 * x192 - x127 * x130 + x128 * x99 - x130 * x192 - x131 * x79 - x131 * x83 - x134 * x55 * xi_xy - x136 * x138 + x136 * x201 - x138 * x151 - x139 * x141 - x141 * x150 + x142 * x143 - x142 * x145 + x142 * x188 + x142 * x189 + x143 * x154 - x145 * x154 - x149 * xi_y + x151 * x201 + x154 * x188 + x154 * x189 + x155 * x183 * x205 - x155 * x60 * x62 + x157 * x58 + x157 * x60 + x159 * x58 + x159 * x60 - x163 * x78 - x163 * x82 + 6 * x166 * x167 + 70 * x166 * x245 + x167 * x170 + x170 * x245 - x171 * x91 - x171 * x98 - x172 * x173 - x173 * x178 + x174 * x175 + x174 * x180 - x174 * x182 + x175 * x179 - x176 * x58 - x176 * x60 - x177 * x58 - x177 * x60 + x179 * x180 - x179 * x182 - x181 * x62 - x181 * x68 - x184 * x66 - x185 * x186 - x185 * x187 + 96 * x185 * x208 * x52 + x185 * x230 + x185 * x231 + x185 * x234 - x186 * x193 - x187 * x193 - x190 * x62 - x190 * x68 - x191 * x62 - x191 * x68 + x193 * x230 + x193 * x231 + x193 * x234 + x194 * x78 + x194 * x82 + x195 * x78 + x195 * x82 + 32 * x196 * x228 * xi_x + x198 * xi_xx + x199 * xi_x + x202 * x203 + x203 * x205 + x204 * xi_xx - x207 * x62 - x207 * x68 - x209 * x62 - x209 * x68 + x211 * x214 * xi_y + x211 * x78 + x211 * x82 - x214 * x224 - x214 * x227 - x214 * x244 - x214 * x246 + x215 * x216 + x215 * x218 + x215 * x219 + x216 * x217 + x217 * x218 + x217 * x219 - x220 * x221 + x220 * x236 * x25 + x220 * x248 + x220 * x249 - x221 * x222 + x222 * x248 + x222 * x249 - x224 * x225 - x225 * x227 - x225 * x244 - x225 * x246 - x233 * xi_y - x235 * xi_y + x237 * x25 * x50 + x238 * x239 + x240 * x241 + x240 * x247 + x241 * x242 + x242 * x247 + x44 * xi_y - x45 * x46 + x45 * x48 - x45 * x61 + x45 * x93 - x46 * x49 + x48 * x49 - x49 * x61 + x49 * x93 - x55 * x71 - x58 * x59 - x58 * x81 - x59 * x60 - x59 * x96 - x60 * x81 + x62 * x63 + x63 * x68 + x65 * x67 + x65 * x79 + x65 * x83 + x69 * xi_yy - 48 * x7 - x72 * x74 + x72 * x77 - x74 * x75 + x75 * x77 + x76 * xi_y + x80 * xi_yy - x84 * x86 - x87 * x88 - x88 * x92 + x89 * x91 + x89 * x98 - x94 * x95) / 2
SS = 4 * A * Bp * Fxp * x0 * x17 * x2 * x3 + A * Bp * Fxp * x11 * x17 * x2 * x284 * x3 + A * Bp * Fyp * x0 * x17 * x2 * x284 * x9 + 2 * A * Bp * Fyp * x2 * x20 * x291 * x6 + 4 * A * Bp * Gp * x11 * x14 * x2 * x262 * x6 + 2 * A * Bp * Gp * x11 * x2 * x278 * x3 * x415 + 2 * A * Bp * Gp * x11 * x2 * x278 * x3 + 4 * A * Bp * Sp * x0 * x21 * x278 * x3 * x313 + 8 * A * Bp * Sp * x11 * x17 * x20 * x21 * x9 + 4 * A * Bp * Sp * x21 * x262 * x291 * x6 + 2 * A * Bp * x0 * x17 * x2 * x284 * x3 * x31 + A * Bp * x0 * x17 * x2 * x284 * x41 * x9 + 4 * A * Bp * x0 * x17 * x2 * x3 * x30 + A * Bp * x0 * x2 * x20 * x284 * x30 * x9 + 2 * A * Bp * x0 * x2 * x20 * x31 * x415 * x9 + 6 * A * Bp * x0 * x2 * x3 * x389 + 2 * A * Bp * x0 * x2 * x329 * x6 + 2 * A * Bp * x0 * x20 * x21 * x284 * x419 * x9 + 4 * A * Bp * x0 * x20 * x21 * x357 * x6 + 4 * A * Bp * x11 * x14 * x17 * x21 * x357 * x9 + 2 * A * Bp * x11 * x17 * x2 * x284 * x29 * x9 + 4 * A * Bp * x11 * x17 * x2 * x3 * x31 + 2 * A * Bp * x11 * x2 * x20 * x29 * x415 * x6 + 4 * A * Bp * x11 * x2 * x20 * x29 * x6 + 2 * A * Bp * x11 * x20 * x21 * x284 * x357 * x6 + 2 * A * Bp * x17 * x2 * x291 * x3 * x30 + 4 * A * Bp * x17 * x21 * x291 * x3 * x419 + 2 * A * Bp * x2 * x20 * x291 * x41 * x6 + 2 * A * Bpp * x0 * x2 * x278 * x3 + A * Bpp * x11 * x2 * x278 * x284 * x3 + 2 * A * Bpp * x2 * x262 * x291 * x6 + 2 * A * Fxp * Gp * x11 * x17 * x2 * x3 + 4 * A * Fxp * Sp * x11 * x20 * x21 * x9 + 2 * A * Fxp * x0 * x2 * x3 * x30 + 2 * A * Fxp * x11 * x2 * x3 * x31 + 2 * A * Fyp * Gp * x11 * x2 * x20 * x6 + 4 * A * Fyp * Sp * x11 * x17 * x21 * x9 + 2 * A * Fyp * x11 * x2 * x29 * x6 + 8 * A * Gp * Sp * x0 * x17 * x20 * x21 * x9 + 2 * A * Gp * x0 * x17 * x2 * x3 * x31 + 2 * A * Gp * x0 * x2 * x20 * x29 * x6 + 2 * A * Gp * x11 * x17 * x2 * x3 * x30 + 2 * A * Gp * x11 * x2 * x3 * x389 + 2 * A * Gp * x11 * x2 * x329 * x6 - A * Gp * x14 * x340 * x440 + 4 * A * Sp * x0 * x17 * x3 * x419 * x50 + 4 * A * Sp * x0 * x20 * x357 * x50 * x6 + A * x0 * x102 * x2 * x6 + 2 * A * x0 * x164 * x2 * x262 * x6 + 2 * A * x0 * x164 * x2 * x278 * x3 + A * x0 * x168 * x2 * x262 * x6 + A * x0 * x168 * x2 * x278 * x3 + A * x0 * x2 * x20 * x284 * x508 * x9 + A * x0 * x2 * x212 * x3 + 2 * A * x0 * x2 * x3 * x366 + 2 * A * x0 * x2 * x3 * x422 + 2 * A * x0 * x2 * x3 * x423 + 2 * A * x0 * x2 * x3 * x621 + 2 * A * x0 * x2 * x31 * x41 * x9 + 2 * A * x0 * x2 * x312 * x6 + 2 * A * x0 * x2 * x363 * x6 + 2 * A * x0 * x2 * x365 * x6 + 4 * A * x0 * x21 * x3 * x30 * x419 + 4 * A * x0 * x21 * x3 * x619 + 4 * A * x0 * x21 * x6 * x614 + 2 * A * x11 * x14 * x17 * x2 * x418 * x9 + A * x11 * x164 * x2 * x262 * x284 * x6 + A * x11 * x164 * x2 * x278 * x284 * x3 + 2 * A * x11 * x17 * x2 * x3 * x457 + 4 * A * x11 * x17 * x21 * x356 * x9 + A * x11 * x2 * x20 * x284 * x418 * x6 + 2 * A * x11 * x2 * x20 * x405 * x6 + 4 * A * x11 * x2 * x3 * x30 * x31 + 2 * A * x11 * x2 * x3 * x623 + 2 * A * x11 * x2 * x6 * x616 + 4 * A * x11 * x20 * x21 * x439 * x9 + 4 * A * x11 * x21 * x29 * x357 * x6 + 4 * A * x11 * x21 * x3 * x31 * x419 + 8 * A * x11 * x357 * x419 * x50 * x9 + 2 * A * x17 * x2 * x291 * x3 * x508 - A * x18 * x250 * x284 * x307 * x357 - A * x23 * x419 * x492 - 12 * A * x253 - A * x26 * x35 * x570 * x574 - A * x357 * x37 * x486 - A * x357 * x40 * x492 + 4 * Ab * x0 * x2 * x3 - Ac * x256 + 4 * Ah * Bp * x0 * x2 * x20 * x313 * x6 - Ah * Fxp * x260 + 4 * Ah * Fyp * x0 * x2 * x6 + 4 * Ah * Gp * x11 * x2 * x20 * x6 + 8 * Ah * Sp * x11 * x17 * x21 * x9 - Ah * x1 * x38 + 4 * Ah * x11 * x2 * x29 * x6 - Ah * x275 * x41 - 4 * Ah * x292 * x7 + 4 * Ap * Bp * x0 * x2 * x278 * x3 - Ap * Bp * x295 + 4 * Ap * Fxp * x0 * x17 * x2 * x3 + 4 * Ap * Fyp * x0 * x2 * x20 * x6 + 4 * Ap * Gp * x11 * x2 * x262 * x6 + 4 * Ap * Gp * x11 * x2 * x278 * x3 + 8 * Ap * Sd * x254 * x9 + 8 * Ap * Sp * x11 * x17 * x20 * x21 * x9 + 4 * Ap * x11 * x2 * x354 * x9 + 4 * Ap * x11 * x2 * x355 * x9 - Ap * x275 * x329 - Ap * x380 * x381 - Ap * x389 * x390 + 4 * App * x0 * x2 * x262 * x6 - App * x274 * x279 + 4 * As * x0 * x2 * x6 + 4 * At * Bp * x17 * x2 * x291 * x3 + 4 * At * Fxp * x0 * x2 * x3 - At * Fyp * x260 + 4 * At * Gp * x11 * x17 * x2 * x3 + 8 * At * Sp * x11 * x20 * x21 * x9 + 4 * At * x0 * x2 * x3 * x30 - At * x1 * x303 + 4 * At * x11 * x2 * x3 * x31 - At * x390 * x414 - Bd ^ 2 * x14 * x255 + 8 * Bd * Bp * x11 * x14 * x17 * x2 * x20 * x9 + 8 * Bd * Gp * x0 * x2 * x262 * x284 * x6 + 8 * Bd * Gp * x11 * x14 * x2 * x278 * x3 + 8 * Bd * Gp * x11 * x2 * x278 * x3 * x415 + 8 * Bd * Gp * x2 * x262 * x339 * x6 + 16 * Bd * Sd * x254 * x9 + 8 * Bd * Sp * x0 * x21 * x278 * x3 * x313 + 16 * Bd * Sp * x0 * x21 * x278 * x3 + 8 * Bd * Sp * x21 * x262 * x291 * x6 - Bd * Sp * x300 * x416 - Bd * Sp * x334 * x348 + 8 * Bd * x0 * x17 * x2 * x284 * x3 * x31 + 4 * Bd * x0 * x17 * x2 * x284 * x41 * x9 + 8 * Bd * x0 * x17 * x2 * x29 * x9 + 4 * Bd * x0 * x2 * x20 * x284 * x30 * x9 + 8 * Bd * x0 * x2 * x20 * x31 * x415 * x9 + 8 * Bd * x0 * x2 * x20 * x31 * x9 + 8 * Bd * x0 * x20 * x21 * x284 * x419 * x9 + 16 * Bd * x11 * x14 * x17 * x21 * x357 * x9 + 8 * Bd * x11 * x17 * x2 * x284 * x29 * x9 - Bd * x11 * x17 * x39 * x41 + 8 * Bd * x11 * x2 * x20 * x29 * x415 * x6 + 8 * Bd * x11 * x2 * x354 * x9 + 8 * Bd * x11 * x2 * x355 * x9 + 8 * Bd * x11 * x20 * x21 * x284 * x357 * x6 + 8 * Bd * x17 * x2 * x291 * x3 * x30 + 16 * Bd * x17 * x21 * x291 * x3 * x419 - Bd * x17 * x29 * x380 * x415 + 8 * Bd * x2 * x20 * x291 * x41 * x6 - Bd * x23 * x383 - 16 * Bd * x291 * x357 * x411 - Bd * x30 * x535 - Bd * x41 * x514 * x515 * x7 - 8 * Bd * x415 * x479 - Bd * x431 * x574 - Bd * x459 * x534 - Bd * x523 * x532 * x573 + 4 * Bp * Fxp * x0 * x20 * x278 * x3 * x339 * x50 + 2 * Bp * Fxp * x17 * x262 * x265 * x50 * x9 + 2 * Bp * Fxp * x265 * x344 * x352 * x50 + Bp * Fyp * x0 * x11 * x20 * x278 * x284 * x50 * x9 + 4 * Bp * Fyp * x11 * x17 * x262 * x291 * x50 * x6 + Bp * Fyp * x269 * x286 * x50 / 2 + 8 * Bp * Gd * x11 * x14 * x2 * x278 * x3 + 8 * Bp * Gd * x11 * x2 * x262 * x6 + 8 * Bp * Gd * x2 * x262 * x339 * x6 + 4 * Bp * Gp * x11 * x291 * x344 * x345 * x50 + 2 * Bp * Gp * x17 * x269 * x286 * x50 * x6 + 4 * Bp * Gp * x20 * x3 * x313 * x352 * x415 * x50 + 4 * Bp * Gp * x20 * x3 * x313 * x352 * x50 + Bp * Gp * x264 * x284 * x50 + Bp * Gp * x264 * x285 * x50 / 2 + 8 * Bp * S * Sp * x0 * x20 * x3 * x339 * x352 + 8 * Bp * S * Sp * x11 * x17 * x269 * x291 * x6 + 4 * Bp * S * Sp * x14 * x264 + 4 * Bp * S * Sp * x20 * x284 * x3 * x352 + Bp * S * Sp * x264 * x286 + 4 * Bp * S * Sp * x265 * x344 * x345 + 8 * Bp * S * x0 * x11 * x17 * x262 * x357 * x6 + 8 * Bp * S * x0 * x17 * x262 * x339 * x586 * x6 + 4 * Bp * S * x11 * x17 * x262 * x291 * x586 * x6 + 12 * Bp * S * x11 * x20 * x278 * x291 * x3 * x419 + 4 * Bp * S * x11 * x269 * x291 * x419 * x6 + 4 * Bp * S * x11 * x291 * x3 * x352 * x586 + 4 * Bp * S * x14 * x17 * x262 * x419 * x9 + 4 * Bp * S * x14 * x344 * x352 * x419 + Bp * S * x17 * x262 * x286 * x419 * x9 + 4 * Bp * S * x20 * x265 * x278 * x357 * x9 + 4 * Bp * S * x265 * x269 * x357 + Bp * S * x286 * x344 * x352 * x419 + 8 * Bp * Sd * x0 * x21 * x262 * x313 * x6 + 8 * Bp * Sd * x21 * x278 * x291 * x3 + 4 * Bp * x0 * x11 * x17 * x20 * x329 * x50 * x6 + Bp * x0 * x11 * x17 * x262 * x284 * x30 * x50 * x9 + Bp * x0 * x11 * x20 * x278 * x284 * x41 * x50 * x9 + 4 * Bp * x0 * x11 * x262 * x354 * x50 * x6 + 3 * Bp * x0 * x11 * x269 * x284 * x41 * x50 + 4 * Bp * x0 * x11 * x278 * x3 * x354 * x50 + 3 * Bp * x0 * x11 * x284 * x30 * x344 * x352 * x50 + 4 * Bp * x0 * x11 * x31 * x344 * x352 * x50 + 8 * Bp * x0 * x20 * x278 * x339 * x50 * x527 * x9 + 6 * Bp * x11 * x17 * x262 * x291 * x41 * x50 * x6 + 6 * Bp * x11 * x20 * x278 * x291 * x3 * x30 * x50 + 4 * Bp * x11 * x20 * x278 * x291 * x50 * x527 * x9 + 2 * Bp * x11 * x269 * x291 * x30 * x50 * x6 + 4 * Bp * x11 * x269 * x291 * x50 * x527 + 2 * Bp * x11 * x291 * x3 * x352 * x41 * x50 + 8 * Bp * x14 * x17 * x20 * x355 * x50 * x9 + 2 * Bp * x14 * x17 * x262 * x29 * x415 * x50 * x6 + 4 * Bp * x14 * x17 * x262 * x30 * x50 * x9 + 2 * Bp * x14 * x20 * x278 * x3 * x31 * x50 + 4 * Bp * x14 * x20 * x278 * x41 * x50 * x9 + 2 * Bp * x14 * x262 * x389 * x50 * x9 + 2 * Bp * x14 * x269 * x31 * x50 * x6 + 4 * Bp * x14 * x269 * x41 * x50 + 2 * Bp * x14 * x278 * x344 * x389 * x50 + 2 * Bp * x14 * x29 * x3 * x352 * x415 * x50 + 4 * Bp * x14 * x30 * x344 * x352 * x50 + (3 // 2) * Bp * x17 * x262 * x285 * x31 * x50 * x9 + 8 * Bp * x17 * x262 * x29 * x313 * x50 * x6 + 4 * Bp * x17 * x262 * x313 * x415 * x50 * x527 * x6 + 3 * Bp * x20 * x278 * x286 * x3 * x31 * x50 - Bp * x267 * x291 * x379 + Bp * x269 * x286 * x31 * x50 * x6 + 4 * Bp * x269 * x313 * x415 * x50 * x6 * x600 + Bp * x285 * x31 * x344 * x352 * x50 / 2 + 4 * Bp * x29 * x3 * x313 * x352 * x415 * x50 - Bp * x29 * x410 - Bp * x330 * x406 - Bp * x332 * x395 - Bp * x354 * x436 - Bp * x526 * x556 + 4 * Bpp * x0 * x20 * x3 * x339 * x352 * x50 + 4 * Bpp * x11 * x17 * x269 * x291 * x50 * x6 + 2 * Bpp * x14 * x264 * x50 + 2 * Bpp * x20 * x284 * x3 * x352 * x50 + Bpp * x264 * x286 * x50 / 2 + 2 * Bpp * x265 * x344 * x345 * x50 - Bpp * x270 * x401 + 8 * Fxc * x20 * x313 * x50 * x9 - Fxc * x20 * x314 + 8 * Fxp * Fyp * x14 * x17 * x20 * x50 * x9 + 3 * Fxp * Gp * x17 * x262 * x284 * x50 * x9 + Fxp * Gp * x284 * x344 * x352 * x50 - Fxp * Gp * x351 + 4 * Fxp * S * Sp * x14 * x17 * x262 * x9 + 4 * Fxp * S * Sp * x14 * x344 * x352 + 8 * Fxp * S * Sp * x17 * x262 * x313 * x9 + 8 * Fxp * S * x0 * x11 * x262 * x357 * x6 + 8 * Fxp * S * x0 * x11 * x278 * x3 * x357 + 8 * Fxp * S * x17 * x20 * x284 * x3 * x419 + 16 * Fxp * Sd * x0 * x17 * x21 * x3 + 4 * Fxp * x0 * x11 * x17 * x20 * x3 * x30 * x50 - Fxp * x101 * x11 * x551 + 8 * Fxp * x14 * x17 * x20 * x41 * x50 * x9 + 8 * Fxp * x14 * x20 * x354 * x50 * x9 + 6 * Fxp * x14 * x278 * x29 * x3 * x50 + 4 * Fxp * x17 * x20 * x3 * x31 * x313 * x50 + 8 * Fxp * x17 * x313 * x329 * x50 * x9 - Fxp * x18 * x461 * x537 / 2 + 8 * Fxp * x20 * x355 * x50 * x9 + 4 * Fxp * x262 * x29 * x313 * x50 * x6 - Fxp * x277 * x303 - Fxp * x462 * x463 / 2 - Fxp * x496 + 8 * Fxs * x14 * x17 * x50 * x9 - Fxs * x17 * x316 + 8 * Fyb * x14 * x20 * x50 * x9 - Fyb * x20 * x316 + 8 * Fyc * x17 * x313 * x50 * x9 - Fyc * x17 * x314 + 3 * Fyp * Gp * x20 * x278 * x284 * x50 * x9 + Fyp * Gp * x269 * x284 * x50 - Fyp * Gp * x410 + 4 * Fyp * S * Sp * x14 * x20 * x278 * x9 + 4 * Fyp * S * Sp * x14 * x269 + 8 * Fyp * S * Sp * x20 * x278 * x313 * x9 + 8 * Fyp * S * x0 * x11 * x262 * x419 * x6 + 8 * Fyp * S * x0 * x11 * x278 * x3 * x419 + 8 * Fyp * S * x17 * x20 * x284 * x357 * x6 + 16 * Fyp * Sd * x0 * x20 * x21 * x6 + 4 * Fyp * x0 * x11 * x262 * x30 * x50 * x6 + 4 * Fyp * x0 * x11 * x278 * x3 * x30 * x50 + 8 * Fyp * x14 * x17 * x355 * x50 * x9 + 6 * Fyp * x14 * x262 * x31 * x50 * x6 + 2 * Fyp * x14 * x262 * x41 * x50 + 2 * Fyp * x14 * x278 * x41 * x50 * x9 + 4 * Fyp * x17 * x20 * x29 * x313 * x50 * x6 + 8 * Fyp * x17 * x354 * x50 * x9 + 8 * Fyp * x20 * x313 * x389 * x50 * x9 - Fyp * x277 * x38 + 4 * Fyp * x278 * x3 * x31 * x313 * x50 - Fyp * x282 - Fyp * x289 * x387 - Fyp * x289 * x41 + 16 * Fypp * x20 * x278 * x313 * x50 * x9 - Gd ^ 2 * x255 + 8 * Gd * Gp * x11 * x17 * x2 * x20 * x9 + 8 * Gd * x0 * x17 * x2 * x41 * x9 + 8 * Gd * x0 * x2 * x354 * x9 + 8 * Gd * x0 * x2 * x355 * x9 + 16 * Gd * x11 * x17 * x21 * x3 * x419 + 16 * Gd * x11 * x20 * x21 * x357 * x6 - Gd * x380 * x474 + 8 * Gp * Sd * x11 * x21 * x262 * x6 + 8 * Gp * Sd * x11 * x21 * x278 * x3 + 4 * Gp * x0 * x11 * x17 * x20 * x354 * x50 * x9 + 4 * Gp * x0 * x11 * x17 * x20 * x355 * x50 * x9 + 2 * Gp * x14 * x17 * x262 * x31 * x50 * x9 + 2 * Gp * x14 * x20 * x278 * x29 * x50 * x9 + 2 * Gp * x14 * x262 * x354 * x50 * x6 + 2 * Gp * x14 * x269 * x29 * x50 + 2 * Gp * x14 * x278 * x3 * x355 * x50 + 2 * Gp * x14 * x31 * x344 * x352 * x50 + Gp * x17 * x262 * x284 * x30 * x50 * x9 + 4 * Gp * x17 * x262 * x31 * x313 * x50 * x9 + 4 * Gp * x17 * x262 * x313 * x41 * x50 * x6 + 4 * Gp * x20 * x278 * x29 * x313 * x50 * x9 + Gp * x262 * x284 * x329 * x50 + Gp * x262 * x284 * x389 * x50 * x9 + 4 * Gp * x262 * x313 * x355 * x50 * x6 + Gp * x278 * x284 * x329 * x50 * x9 + Gp * x278 * x284 * x344 * x389 * x50 + 4 * Gp * x278 * x3 * x313 * x354 * x50 - Gp * x278 * x480 * x5 + Gp * x284 * x30 * x344 * x352 * x50 - Gp * x29 * x353 * x466 - Gp * x373 * x374 + 8 * S * Sp * x0 * x11 * x17 * x262 * x41 * x6 + 4 * S * Sp * x0 * x11 * x269 * x29 + 4 * S * Sp * x0 * x11 * x31 * x344 * x352 + 4 * S * Sp * x14 * x17 * x262 * x30 * x9 + 4 * S * Sp * x14 * x262 * x329 + 4 * S * Sp * x14 * x262 * x389 * x9 + 4 * S * Sp * x14 * x278 * x329 * x9 + 4 * S * Sp * x14 * x278 * x344 * x389 + 4 * S * Sp * x14 * x30 * x344 * x352 + 8 * S * Sp * x17 * x20 * x313 * x354 * x9 + 8 * S * Sp * x17 * x20 * x313 * x355 * x9 + 6 * S * Sp * x17 * x262 * x284 * x31 * x9 + 6 * S * Sp * x20 * x278 * x284 * x29 * x9 + 8 * S * x0 * x11 * x17 * x20 * x3 * x30 * x419 + 8 * S * x0 * x11 * x17 * x20 * x3 * x619 + 8 * S * x0 * x11 * x17 * x20 * x6 * x614 + 8 * S * x0 * x11 * x262 * x6 * x626 + 8 * S * x0 * x11 * x278 * x3 * x626 + 4 * S * x14 * x17 * x262 * x439 * x9 + 16 * S * x14 * x17 * x355 * x357 * x9 + 4 * S * x14 * x20 * x278 * x356 * x9 + 16 * S * x14 * x20 * x354 * x419 * x9 + 4 * S * x14 * x262 * x29 * x419 * x6 + 4 * S * x14 * x262 * x31 * x357 * x6 + 4 * S * x14 * x262 * x357 * x41 + 4 * S * x14 * x269 * x356 + 4 * S * x14 * x278 * x29 * x3 * x419 + 4 * S * x14 * x278 * x3 * x31 * x357 + 4 * S * x14 * x278 * x357 * x41 * x9 + 4 * S * x14 * x344 * x352 * x439 + 8 * S * x17 * x20 * x29 * x313 * x357 * x6 + 8 * S * x17 * x20 * x3 * x31 * x313 * x419 + 8 * S * x17 * x262 * x313 * x439 * x9 + 16 * S * x17 * x313 * x354 * x357 * x9 + 8 * S * x20 * x278 * x313 * x356 * x9 + 16 * S * x20 * x313 * x355 * x419 * x9 - S * x286 * x35 * x357 * x579 - S * x548 * x549 * x6 + 32 * Sd ^ 2 * x2 * x9 + 8 * Sd * Sp * x0 * x262 * x50 * x6 + 8 * Sd * Sp * x0 * x278 * x3 * x50 + 8 * Sd * x0 * x17 * x21 * x29 * x9 + 8 * Sd * x0 * x20 * x21 * x31 * x9 + 8 * Sd * x0 * x20 * x21 * x41 * x6 - Sd * x104 * x20 * x357 + 16 * Sd * x11 * x17 * x357 * x50 * x9 + 16 * Sd * x11 * x20 * x419 * x50 * x9 + 8 * Sd * x11 * x21 * x354 * x9 + 8 * Sd * x11 * x21 * x355 * x9 - Sd * x17 * x24 * x477 * x50 - Sd * x21 * x25 * x381 - Sd * x238 * x386 + 4 * Sp * x0 * x11 * x269 * x419 * x6 + 4 * Sp * x0 * x11 * x3 * x352 * x357 - Sp * x119 * x18 * x269 * x452 + 6 * Sp * x17 * x262 * x284 * x357 * x6 + 6 * Sp * x20 * x278 * x284 * x3 * x419 - Sp * x296 * x446 - Sp * x298 * x347 * x385 - Sp * x343 * x375 + 8 * kappa * x0 * x17 * x2 * x29 * x9 + 8 * kappa * x0 * x2 * x20 * x31 * x9 + 8 * kappa * x0 * x2 * x20 * x41 * x6 + 8 * kappa * x11 * x2 * x354 * x9 + 8 * kappa * x11 * x2 * x355 * x9 - kappa * x329 * x8 + 4 * x0 * x11 * x17 * x20 * x3 * x422 * x50 + 4 * x0 * x11 * x17 * x20 * x3 * x423 * x50 + 4 * x0 * x11 * x17 * x20 * x3 * x50 * x621 + 4 * x0 * x11 * x17 * x20 * x31 * x41 * x50 * x9 + 4 * x0 * x11 * x17 * x20 * x363 * x50 * x6 + 4 * x0 * x11 * x17 * x20 * x365 * x50 * x6 + x0 * x11 * x17 * x262 * x284 * x50 * x508 * x9 + 8 * x0 * x11 * x17 * x262 * x418 * x50 * x6 + 6 * x0 * x11 * x17 * x262 * x457 * x50 * x9 + 6 * x0 * x11 * x20 * x278 * x405 * x50 * x9 + 4 * x0 * x11 * x262 * x29 * x31 * x50 * x6 + 4 * x0 * x11 * x262 * x29 * x41 * x50 + 4 * x0 * x11 * x278 * x29 * x3 * x31 * x50 + 4 * x0 * x11 * x278 * x29 * x41 * x50 * x9 + 3 * x0 * x11 * x284 * x344 * x352 * x50 * x508 + 4 * x0 * x164 * x17 * x269 * x339 * x50 * x6 + 4 * x0 * x164 * x20 * x3 * x339 * x352 * x50 + 8 * x0 * x17 * x2 * x3 * x432 + 16 * x0 * x17 * x21 * x3 * x430 + 4 * x0 * x17 * x262 * x339 * x50 * x595 * x6 + 4 * x0 * x2 * x20 * x284 * x433 * x9 + 8 * x0 * x2 * x20 * x370 * x6 + 8 * x0 * x2 * x20 * x371 * x6 + 16 * x0 * x20 * x21 * x369 * x6 - 6 * x0 * x483 * x491 - x0 * x611 * x612 - x1 * x307 * x471 * x480 - x101 * x164 * x17 * x306 * x342 - x101 * x339 * x409 * x464 * x6 - x101 * x349 * x418 - 3 * x102 * x14 * x200 * x278 + x102 * x14 * x262 * x50 + 4 * x102 * x278 * x313 * x50 * x9 + 4 * x102 * x278 * x50 * x9 - x102 * x318 * x377 + 8 * x11 * x14 * x17 * x2 * x371 * x9 - x11 * x146 * x18 * x262 * x41 + 24 * x11 * x17 * x2 * x20 * x9 + 2 * x11 * x17 * x262 * x291 * x50 * x595 * x6 + 4 * x11 * x2 * x20 * x284 * x371 * x6 + 6 * x11 * x20 * x278 * x291 * x3 * x50 * x508 + 2 * x11 * x269 * x291 * x50 * x508 * x6 - x11 * x28 * x41 * x442 - x11 * x284 * x573 * x589 + 2 * x11 * x291 * x3 * x352 * x50 * x595 - x11 * x300 * x327 - x11 * x327 * x334 - x11 * x330 * x391 - x11 * x395 * x407 * x7 - x11 * x507 * x522 - x11 * x590 * x597 - 12 * x116 * x262 * x357 * x378 - x116 * x269 * x538 * x539 - x118 * x348 * x35 * x353 - x119 * x269 * x439 * x440 - x119 * x329 * x524 - x119 * x377 * x505 * x523 - x128 * x305 * x353 - x134 * x306 * x564 * x587 - x134 * x584 * x599 - x137 * x250 * x31 * x475 - x137 * x269 * x565 + 2 * x14 * x164 * x264 * x313 * x50 + 2 * x14 * x164 * x313 * x344 * x345 * x50 - x14 * x164 * x443 - x14 * x165 * x272 + 8 * x14 * x17 * x20 * x364 * x50 * x9 + 8 * x14 * x17 * x20 * x367 * x50 * x9 + 4 * x14 * x17 * x262 * x50 * x508 * x9 + x14 * x212 * x278 * x344 * x50 + 2 * x14 * x262 * x29 * x30 * x50 * x6 + 2 * x14 * x262 * x312 * x50 + 4 * x14 * x262 * x361 + 4 * x14 * x262 * x420 * x9 + 4 * x14 * x262 * x50 * x6 * x629 + 2 * x14 * x262 * x50 * x615 - x14 * x267 * x343 * x346 - x14 * x267 * x419 * x463 - x14 * x268 * x512 + 2 * x14 * x278 * x29 * x3 * x30 * x50 + 4 * x14 * x278 * x3 * x50 * x629 + 2 * x14 * x278 * x344 * x366 * x50 + 4 * x14 * x278 * x344 * x420 + 4 * x14 * x278 * x361 * x9 + 2 * x14 * x278 * x50 * x615 * x9 + 4 * x14 * x344 * x352 * x50 * x508 - x148 * x262 * x318 - x148 * x271 - x148 * x278 * x376 - x15 * x18 * x20 * x424 * x9 - x16 * x270 * x378 - x160 * x18 * x262 * x339 * x600 - x160 * x519 * x629 - x162 * x374 * x408 + 2 * x164 * x17 * x269 * x284 * x50 * x6 + 2 * x164 * x20 * x284 * x3 * x352 * x50 - x164 * x213 * x278 * x395 + 4 * x164 * x262 * x265 * x278 * x50 * x9 + 2 * x164 * x264 * x265 * x50 + 2 * x164 * x265 * x344 * x345 * x50 - x164 * x291 * x514 * x516 - x165 * x264 * x286 - x165 * x346 * x368 - x168 * x443 + 8 * x17 * x2 * x291 * x3 * x433 + 8 * x17 * x20 * x3 * x30 * x31 * x313 * x50 + 4 * x17 * x20 * x3 * x313 * x50 * x623 + 4 * x17 * x20 * x313 * x50 * x6 * x616 - 6 * x17 * x291 * x418 * x483 - x17 * x302 * x32 * x412 - x17 * x354 * x541 * x598 - x17 * x355 * x543 - x17 * x439 * x454 - x17 * x508 * x594 - x18 * x268 * x384 * x6 - x18 * x3 * x389 * x553 - x18 * x342 * x600 * x612 - x18 * x41 * x450 - x182 * x305 * x306 + 2 * x20 * x265 * x278 * x418 * x50 * x9 - x20 * x29 * x480 * x8 - x20 * x298 * x32 * x592 + x20 * x313 * x357 * x419 * x9 * (x90 + x97) - x20 * x348 * x371 * x7 - x20 * x356 * x404 - x20 * x413 * x42 * x7 - x20 * x439 * x446 - x20 * x450 * x508 - 3 // 2 * x200 * x278 * x285 * x29 * x35 + 4 * x212 * x262 * x313 * x50 * x9 + 4 * x212 * x262 * x50 * x9 - x212 * x376 * x377 - 3 * x212 * x427 - x213 * x313 * x317 - x213 * x455 * x469 - x22 * x294 - x23 * x30 * x445 * x523 * x64 - x23 * x31 * x393 * x497 - x23 * x357 * x571 - x23 * x393 * x41 * x511 - x23 * x547 * x588 - x23 * x566 * x569 - x236 * x265 * x272 - x236 * x346 * x368 - x236 * x375 / 2 - x252 * x397 * x511 - x252 * x509 * x510 - x254 * x258 * x70 - x258 * x26 * x626 - x259 * x274 * x629 - 8 * x259 * x84 - x261 * x263 - x261 * x279 - 3 * x262 * x286 * x487 * x85 + 8 * x262 * x313 * x366 * x50 * x9 - x262 * x453 * x579 - x262 * x566 * x567 - x265 * x266 * x268 + 2 * x265 * x269 * x418 * x50 - x265 * x393 * x485 * x537 - x265 * x394 * x419 * x463 - x266 * x341 * x342 - x267 * x358 + x269 * x284 * x405 * x50 - x269 * x286 * x343 * x357 + x269 * x286 * x418 * x50 - x269 * x415 * x581 * x596 - x271 * x28 * x31 - x271 * x284 * x490 - x271 * x364 - x271 * x367 - x271 * x491 - x273 * x276 - x273 * x494 - x274 * x279 * x310 - x276 * x441 + 8 * x278 * x312 * x313 * x50 * x9 - x278 * x477 * x523 * x591 - x28 * x301 - 2 * x28 * x353 * x462 - x28 * x476 * x579 - 6 * x28 * x588 * x601 - x281 * x286 * x41 - x281 * x29 * x408 - x282 * x41 + x284 * x344 * x352 * x457 * x50 - x284 * x435 * x460 - x284 * x455 * x579 * x601 - 3 * x284 * x464 * x530 - x284 * x536 * x582 - x289 * x29 * x522 - x289 * x615 - x29 * x30 * x324 - x29 * x301 - 8 * x29 * x41 * x561 * x85 - x29 * x419 * x554 - x290 * x313 - x290 - x291 * x442 * x455 - x293 * x371 * x532 - x294 * x31 - x295 * x310 - x296 * x297 - x296 * x309 - x297 * x298 - x298 * x309 - x298 * x313 * x546 - x298 * x314 * x474 - x298 * x41 * x553 * x6 - x299 * x300 - x299 * x334 - x3 * x366 * x520 - x3 * x381 * x389 * x447 - x3 * x389 * x524 * x64 - x30 * x377 * x593 - x30 * x496 - x31 * x357 * x554 - x31 * x398 * x6 - 16 * x311 * x470 - 6 * x312 * x470 - x312 * x521 - x313 * x323 - x313 * x341 * x355 * x445 - x313 * x353 * x415 * x465 * x527 - x313 * x425 - x313 * x429 - x314 * x328 * x389 - x314 * x329 * x388 - x315 * x317 - x316 * x377 * x555 - x319 * x407 * x415 - x32 * x357 * x610 - x320 * x415 - x320 - x322 * x379 - x323 - x324 * x325 - x324 * x326 - x324 * x402 - x324 * x403 - x324 * x456 - x324 * x458 - x325 * x498 - x326 * x448 - x326 * x449 - x329 * x332 - x329 * x333 * x37 - x329 * x381 * x581 - x331 * x340 - x335 * x336 - x336 * x370 - x337 * x338 - x337 * x372 - x338 * x401 - x341 * x347 * x528 - x342 * x417 * x516 - x342 * x463 * x465 * x600 - x346 * x453 - x347 * x348 * x71 - x347 * x357 * x419 * x585 - x347 * x517 * x518 - x347 * x540 - x347 * x547 * x548 - x35 * x393 * x469 * x541 - x35 * x415 * x495 * x600 - x35 * x445 * x503 * x587 - x350 * x384 * x409 - x351 * x415 * x596 - x351 * x457 - x353 * x356 * x440 * x64 - x353 * x357 * x539 * x54 + 4 * x354 ^ 2 * x50 * x9 - x354 * x438 + 4 * x355 ^ 2 * x50 * x9 - x355 * x395 * x529 * x85 - x355 * x438 - x355 * x545 * x598 - x356 * x385 * x513 - x357 * x506 * x560 - x358 * x394 - x360 * x362 - x360 * x421 - x362 * x572 - x363 * x481 - x363 * x583 - x364 * x484 - x364 * x531 - x364 * x562 - x365 * x481 - x365 * x583 - 6 * x366 * x427 - x367 * x484 - x367 * x531 - x367 * x562 - x369 * x383 - x372 * x401 - x382 * x4 * x460 - 8 * x387 * x50 * x588 - x389 * x399 - x389 * x400 - x391 * x392 - x392 * x4 * x459 - x396 * x397 - x396 * x510 * x9 - x398 * x41 - x399 * x459 - x400 * x459 - x402 * x449 - x402 * x495 - x403 * x448 - x403 * x498 - x404 * x505 - x405 * x410 - x405 * x447 * x513 - x406 * x407 - x411 * x412 * x42 - x413 * x479 - x415 * x444 - x415 * x478 * x589 - x415 * x488 - x415 * x490 * x498 - x416 * x445 * x71 - x418 * x489 * x532 - x419 * x42 * x610 - x42 * x464 * x465 - x42 * x569 * x570 * x64 - x421 * x572 - x422 * x557 - x422 * x559 - x423 * x557 - x423 * x559 - x424 * x428 - x424 * x497 * x529 - x425 - x429 - x430 * x431 - x432 * x435 - x433 * x535 - x434 * x5 - x434 * x534 - x436 * x437 - x436 * x555 - x441 * x494 - x444 - x445 * x517 * x518 - x445 * x540 - x449 * x456 - x449 * x487 - x449 * x605 - x450 * x451 - x450 * x590 - x451 * x504 - x452 * x591 * x607 - x454 * x485 - x458 * x495 - x458 * x498 - x461 * x525 - x462 * x512 - x462 * x556 - x467 * x468 - x467 * x616 - x468 * x580 - x471 * x472 - x472 * x473 - x475 * x476 - x477 * x569 * x585 - x481 * x482 - x482 * x583 - x483 * x564 * x608 - x486 * x54 * x607 - x488 - 6 * x489 * x490 - x498 * x605 - x499 * x500 - x499 * x550 - x500 * x501 - x501 * x550 - x502 * x575 * x611 - x503 * x563 * x609 - x504 * x608 - x506 * x507 - x508 * x514 * x597 - x508 * x532 * x536 - 2 * x52 * x551 * x587 - x52 * x565 * x609 - x521 * x615 - x525 * x526 - x528 * x546 - x530 * x584 - x542 * x561 * x626 * x9 - x543 * x544 - x544 * x593 - x549 * x617 - x557 * x558 - x557 * x621 - x558 * x559 - x559 * x621 - x560 * x570 * x606 - x567 * x568 - x571 * x617 - 4 * x575 * x606 * x64 - x576 * x577 - x576 * x623 - x577 * x578 - x578 * x623 - x580 * x616 - x582 * x594 - x602 * x603 - x602 * x619 - x603 * x604 - x604 * x619

   

    return axx, ayy, axy, bx, by, cc, SS
end
