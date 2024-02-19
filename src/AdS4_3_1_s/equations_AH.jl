
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
x0 = cosh(G)
x1 = S * x0
x2 = exp(B) / 2
x3 = exp(-B) / 2
x4 = sinh(G)
x5 = S * x4
x6 = x5 / 2
x7 = x1 / 2
x8 = Sp * x4
x9 = sigma0_y * x8
x10 = Gp * x1
x11 = Fy + xi_y
x12 = x11 * x8
x13 = Bt * x1
x14 = Fxp * x1
x15 = Gt * x5
x16 = Sp * x0
x17 = 2 * Fx
x18 = Bp * S
x19 = x0 * x18
x20 = Gp * x5
x21 = -Sp + x18
x22 = 2 * sigma0_x
x23 = sigma0_x * x8
x24 = sigma0_x * x10
x25 = Fx + xi_x
x26 = x25 * x8
x27 = x10 * x25
x28 = 2 * Fy
x29 = Fyp * x1
x30 = Gh * x5
x31 = 2 * x0 * (Sp + x18) - 2 * x20
x32 = Gh / 2
x33 = Fxp / 2
x34 = Fxpp * x6
x35 = Gt / 2
x36 = sigma0_x / 2
x37 = Fyp * Sp
x38 = x37 * x4
x39 = Fypp * x6
x40 = x16 * x32
x41 = Gph * S
x42 = x0 * x41
x43 = Gpt * x7
x44 = x16 * x35
x45 = Spp * x4
x46 = sigma0_x * sigma0_y
x47 = Gpp * x1
x48 = Fxh + xi_xy
x49 = x8 / 2
x50 = Fyt + xi_xy
x51 = (3 // 2) * Gp
x52 = x14 * x51
x53 = x29 * x51
x54 = x20 * x32
x55 = x20 * x35
x56 = x25 / 2
x57 = Gp * x7
x58 = x11 * x45
x59 = sigma0_y * x25
x60 = Gp ^ 2
x61 = x5 * x60
x62 = x11 * x47
x63 = x11 * x61
x64 = Bpt * x1
x65 = Bt * x16
x66 = Fxpp * x1
x67 = Gpt * x5
x68 = Fxt + xi_xx
x69 = Bt * x19
x70 = Bp * x15
x71 = Bt * x20
x72 = sigma0_x ^ 2
x73 = Bpp * S
x74 = x0 * x73
x75 = Gpp * x5
x76 = 3 * Fxp
x77 = x19 * x76
x78 = x20 * x76
x79 = x25 ^ 2
x80 = Bp ^ 2 * S
x81 = x0 * x80
x82 = S * x60
x83 = x0 * x82
x84 = 2 * Bp
x85 = x20 * x84
x86 = x22 * x25
x87 = Fypp * S
x88 = x0 * x87
x89 = Gh * Sp
x90 = x4 * x89
x91 = Gph * x5
x92 = Bh * x19
x93 = Fy * Gp
x94 = Fy ^ 2
x95 = xi_y ^ 2
x96 = 3 * Fyp
x97 = x20 * x96
x98 = x28 * xi_y
x99 = Gp * x84 - Gpp
x100 = -Spp - x73 + x80 + x82
x101 = Gp * S
axx = -x1 * x2
ayy = -x1 * x3
axy = x5
bx = Fyp * x6 + Gh * x7 + sigma0_y * x10 + x10 * x11 - x12 - x2 * (x13 + x14 + x15 - x16 * x17 + x17 * x19 + x17 * x20 + x22 * (x0 * x21 + x20) + 2 * xi_x * (-x16 + x19 + x20)) - x9
by = Fxp * x6 + Gt * x7 - x23 + x24 - x26 + x27 + x3 * (Bh * x1 + sigma0_y * x31 + x16 * x28 + x19 * x28 - x20 * x28 - x29 - x30 + x31 * xi_y)
cc = Fxp * Fyp * x5 + Fxph * x6 + Fypt * x6 + S ^ 2 * Sdp + 2 * S * Sd * Sp + sigma0_x * x39 + sigma0_x * x40 + sigma0_x * x53 + sigma0_x * x54 - sigma0_x * x58 + sigma0_x * x62 + sigma0_x * x63 + sigma0_xy * x10 + sigma0_xy * x8 + sigma0_y * x34 + sigma0_y * x43 + sigma0_y * x44 + sigma0_y * x52 + sigma0_y * x55 + x11 * x34 + x11 * x43 + x11 * x44 + x11 * x52 + x11 * x55 - x12 * x33 + x14 * x32 + x2 * (-4 * Bp * sigma0_x * x20 * x25 - Fxp ^ 2 * x1 + Fxp * Sp * sigma0_x * x0 + Fxp * Sp * x0 * x25 - Fxp * x13 - Fxp * x15 - Fxpt * x1 - Gt * x23 - Gt * x24 - Gt * x26 - Gt * x27 + 2 * Spp * sigma0_x * x0 * x25 + Spp * x0 * x72 + Spp * x0 * x79 - sigma0_x * x64 - sigma0_x * x65 - sigma0_x * x66 - sigma0_x * x67 - sigma0_x * x69 - sigma0_x * x70 - sigma0_x * x71 - sigma0_x * x77 - sigma0_x * x78 - sigma0_xx * x16 - sigma0_xx * x19 - sigma0_xx * x20 - x16 * x68 - x19 * x68 - x20 * x68 - x25 * x64 - x25 * x65 - x25 * x66 - x25 * x67 - x25 * x69 - x25 * x70 - x25 * x71 - x25 * x77 - x25 * x78 - x72 * x74 - x72 * x75 - x72 * x81 - x72 * x83 - x72 * x85 - x74 * x79 - x74 * x86 - x75 * x79 - x75 * x86 - x79 * x81 - x79 * x83 - x79 * x85 - x81 * x86 - x83 * x86) + x25 * x39 + x25 * x40 + x25 * x53 + x25 * x54 - x25 * x58 + x25 * x62 + x25 * x63 + x29 * x35 + x3 * (Bh * Fy * Gp * S * x4 + Bh * Fy * Sp * x0 + Bh * Fyp * S * x0 + Bh * Gp * S * x4 * xi_y + Bh * Sp * x0 * xi_y + 3 * Bp * Fy * Fyp * S * x0 + Bp * Fy * Gh * S * x4 + 4 * Bp * Fy * Gp * S * x4 * xi_y + Bp * Fyh * S * x0 + 3 * Bp * Fyp * S * x0 * xi_y + Bp * Gh * S * x4 * xi_y + 2 * Bp * Gp * S * x4 * x94 + 2 * Bp * Gp * S * x4 * x95 + Bp * S * x0 * xi_yy + Bph * Fy * S * x0 + Bph * S * x0 * xi_y + 2 * Bpp * Fy * S * x0 * xi_y + Bpp * S * x0 * x94 + Bpp * S * x0 * x95 + Fy * Fyp * Sp * x0 + 2 * Fy * Spp * x0 * xi_y - Fy * x88 - Fy * x90 - Fy * x91 - Fy * x92 - Fy * x97 - Fyh * x16 - Fyh * x20 - Fyp ^ 2 * x1 + Fyp * Sp * x0 * xi_y - Fyp * x30 - Fyph * x1 - Gh * x1 * x93 - Gh * x10 * xi_y + Spp * x0 * x94 + Spp * x0 * x95 + sigma0_y ^ 2 * (-x0 * x100 + x5 * x99) + sigma0_y * (x0 * (-Bh * x21 + Bph * S - Gh * x101 + Spp * x28 + x18 * x96 + x28 * x73 - x28 * x80 - x28 * x82 + x37 - x87) + x4 * (Bh * Gp * S + Bp * S * (Gh + 4 * x93) - Gpp * S * x28 - x101 * x96 - x41 - x89) - 2 * xi_y * (x0 * x100 - x5 * x99)) + sigma0_yy * (x0 * x21 - x20) - x16 * xi_yy - x20 * xi_yy - x75 * x94 - x75 * x95 - x75 * x98 - x81 * x94 - x81 * x95 - x81 * x98 - x83 * x94 - x83 * x95 - x83 * x98 - x88 * xi_y - x90 * xi_y - x91 * xi_y - x92 * xi_y - x97 * xi_y) - x33 * x9 - x36 * x38 + x36 * x42 - x38 * x56 + x42 * x56 - x45 * x46 - x45 * x59 + x46 * x47 + x46 * x61 + x47 * x59 + x48 * x49 + x48 * x57 + x49 * x50 + x50 * x57 + x59 * x61



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

 x0 = sinh(G)
x1 = S * x0
x2 = x1 / 2
x3 = Fxp * x2
x4 = Fyp * x2
x5 = cosh(G)
x6 = S * x5
x7 = x6 / 2
x8 = Gh * x7
x9 = Gt * x7
x10 = Sp * x0
x11 = sigma0_x * sigma0_y
x12 = Gp * x6
x13 = Fy + xi_y
x14 = Fx + xi_x
x15 = x10 * x13
x16 = sigma0_y * x14
x17 = x12 * x13
x18 = Fyp * x6
x19 = Gh * x1
x20 = sigma0_y ^ 2
x21 = Gp * x1
x22 = x13 ^ 2
x23 = Bt * x6
x24 = Fxp * x6
x25 = Gt * x1
x26 = sigma0_x ^ 2
x27 = Bp * x6
x28 = x14 ^ 2
x29 = 2 * sigma0_x * x14


(S ^ 2 * Sd - sigma0_x * x15 + sigma0_x * x17 + sigma0_x * x4 + sigma0_x * x8 + sigma0_xy * x1 + sigma0_y * x3 + sigma0_y * x9 - x10 * x11 - x10 * x16 + x11 * x12 + x12 * x16 + x13 * x3 + x13 * x9 - x14 * x15 + x14 * x17 + x14 * x4 + x14 * x8 + x2 * (Fxh + xi_xy) + x2 * (Fyt + xi_xy) + (Bh * S * sigma0_y * x5 + Bh * S * x13 * x5 + 2 * Bp * S * sigma0_y * x13 * x5 + Bp * S * x20 * x5 + Bp * S * x22 * x5 + 2 * Sp * sigma0_y * x13 * x5 + Sp * x20 * x5 + Sp * x22 * x5 - 2 * sigma0_y * x13 * x21 - sigma0_y * x18 - sigma0_y * x19 - sigma0_yy * x6 - x13 * x18 - x13 * x19 - x20 * x21 - x21 * x22 - x6 * (Fyh + xi_yy)) * exp(-B) / 2 + (2 * Sp * sigma0_x * x14 * x5 + Sp * x26 * x5 + Sp * x28 * x5 - sigma0_x * x23 - sigma0_x * x24 - sigma0_x * x25 - sigma0_xx * x6 - x14 * x23 - x14 * x24 - x14 * x25 - x21 * x26 - x21 * x28 - x21 * x29 - x26 * x27 - x27 * x28 - x27 * x29 - x6 * (Fxt + xi_xx)) * exp(B) / 2)

end
