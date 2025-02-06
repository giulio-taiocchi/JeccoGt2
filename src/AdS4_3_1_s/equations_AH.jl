
function AH_eq_coeff(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd ,
        Bp  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        Bpp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	B_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
        Bp_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        Bp_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

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


x0 = exp(2 * B)
x1 = cosh(G)
x2 = S * x1
x3 = exp(B)
x4 = sinh(G)
x5 = S * x4
x6 = 2 * x5
x7 = Fyp * x5
x8 = Gh * x2
x9 = Sp * x4
x10 = 2 * x9
x11 = sigma0_y * x10
x12 = Gp * x2
x13 = 2 * x12
x14 = Fy + xi_y
x15 = x10 * x14
x16 = Bt * x2
x17 = Fxp * x2
x18 = Gt * x5
x19 = 2 * sigma0_x
x20 = Bp * x2
x21 = Gp * x6
x22 = Fx + xi_x
x23 = 2 * x22
x24 = Bh * x2
x25 = Fyp * x2
x26 = Gh * x5
x27 = Sp * x1
x28 = 2 * x27
x29 = sigma0_y * x28
x30 = 2 * x20
x31 = sigma0_y * x21
x32 = Fxp * x5
x33 = Gt * x2
x34 = x12 * x19
x35 = x12 * x23
x36 = Bh * x27
x37 = Bph * x2
x38 = Fyp * x27
x39 = sigma0_y ^ 2
x40 = Spp * x1
x41 = Fypp * x2
x42 = Gh * x9
x43 = Gp * x5
x44 = Gph * x5
x45 = Bh * x43
x46 = Fyh + xi_yy
x47 = Bp * x27
x48 = Bpp * x2
x49 = Bp * sigma0_y
x50 = 2 * x25
x51 = 3 * Gp
x52 = x51 * x7
x53 = Gp * x8
x54 = x14 ^ 2
x55 = Gpp * x5
x56 = sigma0_y * x14
x57 = 2 * x56
x58 = Bp * x43
x59 = Bp * x14
x60 = Gp ^ 2
x61 = x2 * x60
x62 = Bpt * x2
x63 = Bt * x27
x64 = Fxpp * x2
x65 = Gpt * x5
x66 = Gt * x9
x67 = Fxt + xi_xx
x68 = Bt * x43
x69 = Gp * x33
x70 = sigma0_x ^ 2
x71 = Bp * x16
x72 = 4 * sigma0_x
x73 = Bp * x17
x74 = Bp * x18
x75 = x32 * x51
x76 = x22 ^ 2
x77 = 3 * x58
x78 = x19 * x22
x79 = sigma0_x * x6
x80 = Bp ^ 2 * x2
x81 = 2 * x80
x82 = 2 * S ^ 2
x83 = Fxpp * x5
x84 = Fypp * x5
x85 = Gh * x27
x86 = Gph * x2
x87 = Gpt * x2
x88 = Gt * x27
x89 = Fxh + xi_xy
x90 = Fyt + xi_xy
x91 = Fxp * x9
x92 = Fyp * x9
x93 = Spp * x4
x94 = sigma0_y * x19
x95 = Bp * x32
x96 = Bp * x7
x97 = Bp * x8
x98 = Bp * x33
x99 = Gp * x26
x100 = Gp * x18
x101 = Bp * x5
x102 = Bp * sigma0_x
x103 = x17 * x51
x104 = x25 * x51
x105 = Gpp * x2
x106 = x14 * x93
x107 = sigma0_y * x23
x108 = Bp * x22
x109 = x60 * x79
x110 = x105 * x14
x111 = x22 * x6 * x60
axx = -x0 * x2
ayy = -x2
axy = x3 * x6
bx = x3 * (sigma0_y * x13 - x11 + x13 * x14 - x15 + x3 * (2 * Sp * sigma0_x * x1 + 2 * Sp * x1 * x22 - sigma0_x * x21 - x16 - x17 - x18 - x19 * x20 - x20 * x23 - x21 * x22) + x7 + x8)
by = sigma0_y * x30 - x14 * x21 + x14 * x28 + x14 * x30 + x24 - x25 - x26 + x29 + x3 * (-sigma0_x * x10 - x10 * x22 + x32 + x33 + x34 + x35) - x31
cc = -Fyp ^ 2 * x2 + Fyp * x24 - Fyph * x2 - Gh * x7 - Gpp * x56 * x6 + sigma0_y * x36 + sigma0_y * x37 + sigma0_y * x38 - sigma0_y * x41 - sigma0_y * x42 - sigma0_y * x44 + sigma0_y * x45 - sigma0_y * x52 - sigma0_y * x53 - sigma0_yy * x27 - sigma0_yy * x43 + x0 * (2 * Bp * Sp * sigma0_x * x1 * x22 + Bp * Sp * x1 * x70 + Bp * Sp * x1 * x76 - Fxp ^ 2 * x2 + Fxp * Sp * sigma0_x * x1 + Fxp * Sp * x1 * x22 - Fxp * x16 - Fxp * x18 - Fxpt * x2 - Gpp * x22 * x79 + 2 * Spp * sigma0_x * x1 * x22 + Spp * x1 * x70 + Spp * x1 * x76 - 6 * sigma0_x * x22 * x58 - sigma0_x * x62 - sigma0_x * x63 - sigma0_x * x64 - sigma0_x * x65 - sigma0_x * x66 - sigma0_x * x68 - sigma0_x * x69 - sigma0_x * x75 - sigma0_xx * x27 - sigma0_xx * x30 - sigma0_xx * x43 - x19 * x71 - x19 * x74 - x22 * x62 - x22 * x63 - x22 * x64 - x22 * x65 - x22 * x66 - x22 * x68 - x22 * x69 - x22 * x72 * x80 - 4 * x22 * x73 - x22 * x75 - x23 * x71 - x23 * x74 - x27 * x67 - x30 * x67 - x43 * x67 - x48 * x70 - x48 * x76 - x48 * x78 - x55 * x70 - x55 * x76 - x61 * x70 - x61 * x76 - x61 * x78 - x70 * x77 - x70 * x81 - x72 * x73 - x76 * x77 - x76 * x81) + x14 * x36 + x14 * x37 + x14 * x38 - x14 * x41 - x14 * x42 - x14 * x44 + x14 * x45 - x14 * x52 - x14 * x53 - x27 * x46 + x29 * x59 + x3 * (Bp * Sd * x82 + Bp * sigma0_xy * x6 + 2 * Fxp * x7 + Fxph * x5 + Fyp * x33 + Fypt * x5 + Gh * x17 + 4 * S * Sd * Sp + Sdp * x82 + sigma0_x * x104 + sigma0_x * x84 + sigma0_x * x85 + sigma0_x * x86 - sigma0_x * x92 + sigma0_x * x96 + sigma0_x * x97 + sigma0_x * x99 + sigma0_xy * x10 + sigma0_xy * x13 + sigma0_y * x100 + sigma0_y * x103 + sigma0_y * x109 + sigma0_y * x111 + sigma0_y * x83 + sigma0_y * x87 + sigma0_y * x88 - sigma0_y * x91 + sigma0_y * x95 + sigma0_y * x98 + x100 * x14 + x101 * x89 + x101 * x90 - x102 * x11 - x102 * x15 + x103 * x14 + x104 * x22 + x105 * x107 + x105 * x94 - x106 * x19 - x106 * x23 - x107 * x93 - x108 * x11 - x108 * x15 + x109 * x14 + x110 * x19 + x110 * x23 + x111 * x14 + x12 * x89 + x12 * x90 + x14 * x83 + x14 * x87 + x14 * x88 - x14 * x91 + x14 * x95 + x14 * x98 + x22 * x84 + x22 * x85 + x22 * x86 - x22 * x92 + x22 * x96 + x22 * x97 + x22 * x99 + x34 * x49 + x34 * x59 + x35 * x49 + x35 * x59 + x89 * x9 + x9 * x90 - x93 * x94) + x31 * x59 + x39 * x40 + x39 * x47 + x39 * x48 - x39 * x55 + x39 * x58 - x39 * x61 + x40 * x54 + x40 * x57 - x43 * x46 + x47 * x54 + x48 * x54 + x48 * x57 + x49 * x50 + x50 * x59 - x54 * x55 + x54 * x58 - x54 * x61 - x57 * x61



    return axx, ayy, axy, bx, by, cc
end

function AH_eq_res(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd ,
        Bp  ,  Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        Bpp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,
        B_x ,  G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
        B_y ,  G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B")
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
res = sigma0_y * x2 - sigma0_y * x5 - sigma0_y * x8 - sigma0_yy * x1 - x1 * (Fyh + xi_yy) + x10 * x11 + x10 * x13 + x10 * x3 - x11 * x12 + x11 * x4 - x12 * x13 - x12 * x3 + x13 * x4 + x2 * x9 + x3 * x4 - x5 * x9 - x8 * x9 + (2 * Sp * sigma0_x * x0 * x18 + Sp * x0 * x17 + Sp * x0 * x19 - sigma0_x * x14 - sigma0_x * x15 - sigma0_x * x16 - sigma0_xx * x1 - x1 * (Fxt + xi_xx) - x10 * x17 - x10 * x19 - x10 * x21 - x12 * x17 - x12 * x19 - x12 * x21 - x14 * x18 - x15 * x18 - x16 * x18) * exp(2 * B) + (2 * S ^ 2 * Sd + sigma0_x * x23 + sigma0_x * x24 + 2 * sigma0_xy * x7 + sigma0_y * x22 + sigma0_y * x25 + x18 * x23 + x18 * x24 - x20 * x29 + x20 * x32 + x22 * x9 + x25 * x9 - x26 * x27 - x26 * x31 + x27 * x28 + x28 * x31 - x29 * x30 + x30 * x32 + x7 * (Fxh + xi_xy) + x7 * (Fyt + xi_xy)) * exp(B)

end
