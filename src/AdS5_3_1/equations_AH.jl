
function AH_eq_coeff(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,
        B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
        B1p_x, B2p_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        B1p_y, B2p_y, Gp_y,  Sp_y , Fxp_y , Fyp_y ,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

    @tilde_outer("B1p")
    @tilde_outer("B2p")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("B1p")
    @hat_outer("B2p")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")

    x0 = exp(B1 + B2)
    x1 = cosh(G)
    x2 = S*x1
    x3 = sinh(G)
    x4 = exp(B2)
    x5 = x3*x4
    x6 = 2*x5
    x7 = -B1
    x8 = exp(B2 + x7)
    x9 = B2h*S
    x10 = Fyp*S
    x11 = Sp*sigma0_y
    x12 = Sp*xi_y
    x13 = 2*B2p
    x14 = S*x13
    x15 = sigma0_y*x14
    x16 = x14*xi_y
    x17 = exp(B1)
    x18 = S*x17
    x19 = Gt*x18
    x20 = Fx*Gp
    x21 = 2*x18*x20
    x22 = 2*sigma0_x
    x23 = Gp*x18
    x24 = x22*x23
    x25 = 2*xi_x
    x26 = x23*x25
    x27 = Fy + xi_y
    x28 = sigma0_y + x27
    x29 = Gp*x28
    x30 = S*(Gh + 2*x29)
    x31 = 2*B1p
    x32 = x13 + x31
    x33 = Fx*S
    x34 = S*x31
    x35 = Fx*Sp
    x36 = B2t*S + Fxp*S - Sp*sigma0_x - Sp*xi_x + St + sigma0_x*x14 + x14*xi_x - x35
    x37 = -Sh
    x38 = Sp*x28
    x39 = x1*x4
    x40 = x28^2
    x41 = x39*x40
    x42 = 3*Spp
    x43 = 3*B2p
    x44 = 3*Sp
    x45 = sigma0_x + xi_x
    x46 = Fx + x45
    x47 = 6*Sp
    x48 = exp(2*B1 + B2)
    x49 = x1*x48
    x50 = x46^2
    x51 = x49*x50
    x52 = Gp*x3
    x53 = Sp*x45
    x54 = St - 2*x35 - 2*x53
    x55 = x17*x3
    x56 = x54*x55
    x57 = B1h + B1p*x28
    x58 = B2p*x28
    x59 = B2h + x58
    x60 = x1*x59
    x61 = Gh + x29
    x62 = x3*x61
    x63 = B2p*x46 + B2t
    x64 = x55*x63
    x65 = Gp*x45 + Gt + x20
    x66 = x1*x17
    x67 = x65*x66
    x68 = x1*x57 - x60 - x62 + x64 + x67
    x69 = S*x68 + x1*(x37 - x38) + x56
    x70 = x4*x69
    x71 = Sh + x38
    x72 = x3*x71
    x73 = St + x35 + x53
    x74 = x66*x73
    x75 = B1p*x46 + B1t
    x76 = x1*x63
    x77 = x3*x65
    x78 = x17*(x1*x75 + x76 + x77)
    x79 = x1*(-Gh - x29) + x3*(-B2h - x58) + x78
    x80 = S*x79 - x72 + x74
    x81 = x0*x80
    x82 = Fxp*x28
    x83 = 2*sigma0_xy + 2*xi_xy
    x84 = x0*x3
    x85 = Fyp*sigma0_y
    x86 = Fyp*x27
    x87 = 6*S
    x88 = Fyp*sigma0_x
    x89 = Fx + xi_x
    x90 = Fyp*x89
    x91 = Fxp*sigma0_x
    x92 = Fxp*x89
    x93 = x0*x1
    x94 = x17*(Sd*x87 + x5*(2*Fyt + x83 + 2*x88 + 2*x90) - x93*(2*Fxt + 2*sigma0_xx + 2*x91 + 2*x92 + 2*xi_xx))
    x95 = -x39*(2*Fyh + 2*sigma0_yy + 2*x85 + 2*x86 + 2*xi_yy) + x84*(2*Fxh + 2*x82 + x83) + x94
    x96 = 2*x28
    x97 = Fyh + sigma0_yy + x85 + x86 + xi_yy
    x98 = sigma0_xy + xi_xy
    x99 = Fxh + x82 + x98
    x100 = 2*x93
    x101 = Fyt + x88 + x90 + x98
    x102 = 2*Gp
    x103 = Fxt + sigma0_xx + x91 + x92 + xi_xx
    x104 = Sph + Spp*x28
    x105 = Gp*x1
    x106 = B2ph + B2pp*x28
    x107 = Gph + Gpp*x28
    x108 = B2pp*x46 + B2pt
    x109 = Fx*Gpp + Gpp*x45 + Gpt
    x110 = 2*Fx
    x111 = x0*(x110 + x22 + x25)
    x112 = Gp*x17
    x113 = x4*(2*Fy + 2*sigma0_y + 2*xi_y)

    axx = -x0*x2

    axy = S*x6

    ayy = -x2*x8

    bx = x4*(x1*(-x17*(B1t*S + sigma0_x*x34 + x32*x33 + x34*xi_x + x36) + x30) + x3*(-Fy*Sp + Fy*x14 + Sh + x10 - x11 - x12 + x15 + x16 - x19 - x21 - x24 - x26 + x9))

    by = x8*(x1*(B1h*S + Fy*(Sp - x14 + x34) + sigma0_y*x34 - x10 + x11 + x12 - x15 - x16 + x19 + x21 + x24 + x26 + x34*xi_y + x37 - x9) + x3*(x17*(x13*x33 + x36) - x30))

    cc = (-B1p*(S*x95 - x111*x80 + x113*x69 + x41*x44 + x44*x51) + Fxp*x46*x47*x49 - 2*Fxp*x81 + 6*Fyp*x38*x39 + 2*Fyp*x70 + Gp*x40*x44*x5 + S*(B1p*x94 + Gp*x100*x99 - Gp*x6*x97 - x13*x39*x97 + x17*(Sd*x47 + Sdp*x87 - x100*x103*(B1p + B2p) - x100*(Fx*Fxpp + Fxpp*x45 + Fxpt) + x101*x102*x39 + x101*x13*x5 - x102*x103*x84 + x6*(Fx*Fypp + Fypp*x45 + Fypt)) + x32*x84*x99 - x39*(2*Fyph + Fypp*x96) + x84*(2*Fxph + Fxpp*x96)) + Sp*x41*x43 + Sp*x51*(6*B1p + x43) + Sp*x95 - x111*(B1p*x74 + Gp*x55*x73 + S*(B1p*x78 - Gp*x60 - Gp*x62 - x1*x107 - x106*x3 + x17*(x1*x108 + x1*(B1pp*x46 + B1pt) + x105*x65 + x109*x3 + x52*x63 + x52*x75)) + Sp*x79 - x104*x3 - x105*x71 + x66*(Fx*Spp + Spp*x45 + Spt)) + x113*(B1p*x56 + Gp*x54*x66 - Gp*x72 + S*(B1p*x64 + B1p*x67 - x1*x106 + x1*(B1ph + B1pp*x28) - x105*x61 - x107*x3 + x108*x55 + x109*x66 + x112*x76 + x112*x77 + x52*x57 - x52*x59) + Sp*x68 - x1*x104 - x55*(Fxp*x44 + Spp*x110 + Spp*x22 + Spp*x25 - Spt)) + x13*x28*x70 - x32*x46*x81 + x41*x42 + x42*x51 + x44*x48*x50*x52)*exp(x7)/2



    return axx, ayy, axy, bx, by, cc
end

function AH_eq_res(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   , B2   , G   ,  S    , Fx    , Fy    , Sd ,
        B1p  , B2p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp , B2pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x , B2_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
        B1_y , B2_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("B2")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("B2")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

    x0 = cosh(G)
x1 = S * x0
x2 = Bh * x1
x3 = sigma0_y ^ 2
x4 = Sp * x0
x5 = Fyp * x1
x6 = sinh(G)
x7 = S * x6
x8 = Gh * x7
x9 = Fy + xi_y
x10 = Bp * x1
x11 = x9 ^ 2
x12 = Gp * x7
x13 = 2 * sigma0_y * x9
x14 = Bt * x1
x15 = Fxp * x1
x16 = Gt * x7
x17 = sigma0_x ^ 2
x18 = Fx + xi_x
x19 = x18 ^ 2
x20 = 2 * sigma0_x
x21 = x18 * x20
x22 = Fxp * x7
x23 = Fyp * x7
x24 = Gh * x1
x25 = Gt * x1
x26 = Sp * x6
x27 = sigma0_y * x20
x28 = Gp * x1
x29 = x26 * x9
x30 = 2 * x18
x31 = sigma0_y * x30
x32 = x28 * x9

sigma0_y * x2 - sigma0_y * x5 - sigma0_y * x8 - sigma0_yy * x1 - x1 * (Fyh + xi_yy) + x10 * x11 + x10 * x13 + x10 * x3 - x11 * x12 + x11 * x4 - x12 * x13 - x12 * x3 + x13 * x4 + x2 * x9 + x3 * x4 - x5 * x9 - x8 * x9 + (2 * Sp * sigma0_x * x0 * x18 + Sp * x0 * x17 + Sp * x0 * x19 - sigma0_x * x14 - sigma0_x * x15 - sigma0_x * x16 - sigma0_xx * x1 - x1 * (Fxt + xi_xx) - x10 * x17 - x10 * x19 - x10 * x21 - x12 * x17 - x12 * x19 - x12 * x21 - x14 * x18 - x15 * x18 - x16 * x18) * exp(2 * B) + (2 * S ^ 2 * Sd + sigma0_x * x23 + sigma0_x * x24 + 2 * sigma0_xy * x7 + sigma0_y * x22 + sigma0_y * x25 + x18 * x23 + x18 * x24 - x20 * x29 + x20 * x32 + x22 * x9 + x25 * x9 - x26 * x27 - x26 * x31 + x27 * x28 + x28 * x31 - x29 * x30 + x30 * x32 + x7 * (Fxh + xi_xy) + x7 * (Fyt + xi_xy)) * exp(B)

end
