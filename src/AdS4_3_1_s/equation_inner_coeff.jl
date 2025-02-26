
#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_inner("B")

should expand to

  Bt = B_x -  (Fx * u + xi_x) * Bp

etc.

=#
macro tilde_inner(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx * u + xi_x) * $fp) )
end
macro hat_inner(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy * u + xi_y) * $fp) )
end

macro bar_inner(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx * u + xi_x) * ( -2*($fp_x) + (Fx * u + xi_x) * ($fpp) )) )
end

macro star_inner(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy * u + xi_y) * ( -2*($fp_y) + (Fy * u + xi_y)* ($fpp) )) )
end

macro cross_inner(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx * u + xi_x) * ($fp_y) -
                  (Fy * u + xi_y) * ( $fp_x -(Fx * u + xi_x) * ($fpp) ) ) )
end



# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (S0, S0_t, u, xi, B, Bp, G, Gp) = vars

  
x0 = u ^ 4
x1 = 3 * u
x2 = (-B * x1 + Bp) ^ 2 * cosh(G * u ^ 3) ^ 2
ABCS[1] = 4 * u ^ 2
ABCS[2] = 24 * u
ABCS[3] = 9 * G ^ 2 * u ^ 6 - 6 * G * Gp * u ^ 5 + Gp ^ 2 * x0 + x0 * x2 + 24
ABCS[4] = u * (x2 + (-G * x1 + Gp) ^ 2) * (u * xi + 1)



    nothing
end


# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        S0, S0_x, S0_y, S0_t, S0_tx, S0_ty, u, xi, xi_x, xi_y,
        B     ,        G      ,       S      ,
        Bp    ,        Gp     ,       Sp     ,
        Bpp   ,        Gpp    ,       Spp    ,
        B_x   ,        G_x    ,       S_x    ,
        B_y   ,        G_y    ,       S_y    ,
        Bp_x  ,        Gp_x   ,       Sp_x   ,
        Bp_y  ,        Gp_y   ,       Sp_y   ,
       
    ) = vars
     
x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = u ^ 2
x4 = u * xi
x5 = S * x0
x6 = x4 + x5 + 1
x7 = x6 ^ 2
x8 = 2 * x3 * x7
x9 = Bp * x3
x10 = 3 * x1
x11 = G * x0
x12 = 2 * x11
x13 = cosh(x12)
x14 = -Bp
x15 = 3 * u
x16 = B * x15
x17 = x14 + x16
x18 = x10 + x13 * x17 * x3 - x9
x19 = u * x2
x20 = 2 * Gp
x21 = 6 * u
x22 = G * x21
x23 = -x17
x24 = sinh(x12)
x25 = x23 * x24
x26 = -x20 + x22 - x25
x27 = x0 * x7
x28 = 4 * Spp
x29 = 8 * xi
x30 = Sp * u
x31 = 2 * Bpp
x32 = B * x3
x33 = S * x3
x34 = xi ^ 2
x35 = u * x34
x36 = Bp ^ 2
x37 = u ^ 4
x38 = Bp * x37
x39 = 6 * B
x40 = u ^ 5
x41 = S * x40
x42 = Sp * x39
x43 = 24 * xi
x44 = Bp * x0
x45 = 2 * Sp
x46 = 12 * xi
x47 = 2 * x5
x48 = G * Gp
x49 = x37 * x48
x50 = Sp * x3
x51 = 16 * xi
x52 = B ^ 2
x53 = 9 * x40
x54 = G ^ 2
x55 = 18 * x54
x56 = Gp ^ 2
x57 = x0 * x56
x58 = S ^ 2
x59 = x40 * x58
x60 = Sp ^ 2 * x0
x61 = u ^ 6
x62 = u ^ 7
x63 = Bp * x62
x64 = B * Bp * x46
x65 = S * x62
x66 = S * x61
x67 = x37 * xi
x68 = x48 * x65
x69 = x40 * x48
x70 = u ^ 8
x71 = 18 * x52
x72 = S * x70
x73 = x61 * xi
x74 = 2 * x36
x75 = 36 * x54
x76 = x56 * x66
x77 = 4 * xi
x78 = x56 * x77
x79 = x48 * x72
x80 = u ^ 9
x81 = x58 * x80
x82 = x34 * x40
x83 = u ^ 10 * x58
x84 = Bp * x39
x85 = x34 * x61
x86 = S * x80 * xi
x87 = x65 * xi
x88 = 12 * x48
x89 = 9 * x52
x90 = u ^ 11 * x58
x91 = x34 * x62
x92 = x56 * x81
x93 = x56 * x82
x94 = -Gp
x95 = G * x15
x96 = x94 + x95
x97 = -x96
x98 = x25 * x27 * x97
x99 = Bpp * x6
x100 = x3 * x36 * x6
x101 = 6 * x1 * x6
x102 = 2 * x50
x103 = -x102
x104 = 7 * x4
x105 = 11 * x5
x106 = x103 + x104 + x105 + 5
x107 = 5 * x4
x108 = 9 * x5
x109 = x103 + x107 + x108 + 3
x110 = u * (Bp * (-x101 - x106) + x100 + x16 * (x10 * x6 + x109)) + x99
x111 = x13 * x6
x112 = 2 * Gpp
x113 = x112 * x6
x114 = x0 * x23
x115 = 2 * x114
x116 = -u * (Bp * (-x101 + x106) + x100 + x16 * (3 * B * x0 * x6 - x109)) + x99
x117 = 3 * B
x118 = x10 + x117 * x66 + x117 * x67 - x6 * x9
x119 = x102 + x118
x120 = 2 * u
x121 = -S * x120 + Sp
x122 = S_x - x121 * xi_x
x123 = -x47 + x50 + 1
x124 = S_y - x121 * xi_y
x125 = 2 * x0
x126 = Spp - 4 * x30 + 6 * x33
x127 = cosh(x11) ^ 2
x128 = G * x120 + x94
x129 = B_y - x23 * xi_y
x130 = G_y - x97 * xi_y
x131 = x115 * x13
x132 = Bpp + x21 * (B * x120 + x14)
x133 = -Bp * x21 + Bpp + 12 * x32
x134 = B_x * x15 - Bp_x
x135 = 2 * x127
x136 = B_x - x23 * xi_x
x137 = x0 * x135 * x17
x138 = G_x - x97 * xi_x
x139 = -x20 + x22 + x25
x140 = 2 * Spp
x141 = 2 * x56
x142 = 6 * x48
x143 = 9 * x54
AA[1,1] = x2 * x8
AA[1,2] = 0
BB[1,1] = x19 * x7 * (x18 + 8)
BB[1,2] = -x26 * x27
CC[1,1] = x19 * (-12 * B * S * x63 + 15 * B * x34 * x37 + 36 * B * x41 + 27 * B * x58 * x70 + 42 * B * x66 * xi - 5 * Bp * u - 18 * Bp * x41 * xi + Bp * x45 * x66 + Bpp * x3 * x34 + Bpp * x47 + Bpp * x58 * x61 + Bpp + S * x31 * x67 - 16 * S * x38 + x0 * x36 + x1 * x43 + x110 * x111 + x28 * x4 + x28 * x5 + x28 + x29 - 24 * x30 + x31 * x4 + 9 * x32 + 48 * x33 - 7 * x34 * x44 + 4 * x35 + x36 * x81 + x36 * x82 - x37 * x42 + x37 * x78 - x38 * x39 + x38 * x45 * xi - x40 * x42 * xi + x40 * x55 - x40 * x64 - x42 * x65 - x43 * x69 - x43 * x79 + x44 * x45 - x46 * x9 - 12 * x49 + 32 * x5 * xi - x50 * x51 + x52 * x53 + x55 * x90 + x55 * x91 + 2 * x57 - 11 * x58 * x63 + 12 * x59 - 4 * x60 - x64 * x72 + x65 * x78 + x66 * x74 + x67 * x74 - 24 * x68 + x71 * x72 + x71 * x73 + x71 * x86 + x72 * x75 + x73 * x75 + x74 * x87 + x75 * x86 + 4 * x76 - x83 * x84 - x83 * x88 - x84 * x85 - x85 * x88 + x89 * x90 + x89 * x91 + 2 * x92 + 2 * x93 + 2 * x98)
CC[1,2] = -u * x6 * (x111 * x115 * x97 + x113 + x116 * x24 - x120 * (Gp * (x104 + x105 - x119 + 5) + x95 * (-x107 - x108 + x119 - 3)))
SS[1] = 4 * u * x122 * x123 * x2 - x124 * x125 * x26 * x6 - x2 * (4 * x4 + 4 * x5 + 4) * (-S_x * x120 + Sp_x + x0 * x122 * x127 * x23 - x126 * xi_x) + x7 * (-G_y * x21 + 2 * Gp_y - x125 * x129 * x97 + x130 * x131 + x2 * (-x125 * x138 * (Gp + x25 - x95) + x135 * (x133 * xi_x + x134) + x136 * x137) + x24 * (-B_y * x15 + Bp_y - x114 * x129 - x132 * xi_y) - xi_y * (12 * u * x128 + x112))
AA[2,1] = 0
AA[2,2] = x8
BB[2,1] = -x139 * x2 * x27
BB[2,2] = -u * x7 * (x18 - 8)
CC[2,1] = x19 * x6 * (2 * x0 * x13 * x23 * x6 * x97 + x110 * x24 - x113 - x120 * (Gp * (-x106 - x118) + x95 * (x109 + x118)))
CC[2,2] = x120 * (-x116 * x127 * x6 + x140 * x4 + x140 * x5 + x140 + x141 * x67 + x141 * x87 - x142 * x83 - x142 * x85 + x143 * x90 + x143 * x91 - x29 * x50 - 12 * x30 + 24 * x33 + 2 * x35 - x46 * x69 - x46 * x79 - 6 * x49 + x5 * x51 + x53 * x54 + x55 * x72 + x55 * x73 + x55 * x86 + x57 + 6 * x59 - 2 * x60 - 12 * x68 + 2 * x76 + x77 + x92 + x93 - x98)
SS[2] = 4 * u * x123 * x124 + x6 * (8 * S_y * u - 4 * Sp_y + 4 * x114 * x124 * x127 - x122 * x125 * x139 * x2 + 4 * x126 * xi_y) + x7 * (x125 * x130 * (x25 + x96) - x127 * (B_y * x21 - 2 * Bp_y + 2 * x133 * xi_y) + x129 * x137 + x2 * (-G_x * x21 + 2 * Gp_x + 2 * x0 * x136 * x97 - x131 * x138 - x24 * (x0 * x136 * x23 - x132 * xi_x - x134) - 2 * xi_x * (Gpp + x128 * x21)))


    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
       S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Sp")
    @tilde_inner("Fxp")
    @tilde_inner("Fyp")
    @tilde_inner("Bp")
    @tilde_inner("Sp")
    @tilde_inner("Gp")

    @hat_inner("Sp")
    @hat_inner("Fxp")
    @hat_inner("Fyp")
    @hat_inner("Gp")
    @hat_inner("Bp")
    @hat_inner("Gp")

    @cross_inner("G")
    @cross_inner("S")

x0 = u ^ 3
x1 = S * x0
x2 = u * xi
x3 = x2 + 1
x4 = x1 + x3
x5 = u ^ 2
x6 = B * x0
x7 = exp(x6)
x8 = 8 * x7
x9 = x5 * x8
x10 = x4 ^ 2
x11 = Sp * u
x12 = S * x5
x13 = 4 * Sp
x14 = x13 * x7
x15 = G * x0
x16 = cosh(x15)
x17 = 2 * x16
x18 = u * x17
x19 = S * u
x20 = 24 * x7
x21 = 16 * x2
x22 = Bs * x17
x23 = x16 * x5
x24 = 2 * x23
x25 = 2 * x0
x26 = sinh(x15)
x27 = Gs * x26
x28 = 4 * x16
x29 = x0 * x28
x30 = Fyp ^ 2
x31 = 6 * B
x32 = u ^ 5
x33 = x16 * x32
x34 = Fyh * x33
x35 = u ^ 4
x36 = x16 * x35
x37 = x36 * xi_yy
x38 = Bh * Fy
x39 = 10 * x33
x40 = Bh * x36
x41 = 2 * Fyp
x42 = u ^ 6
x43 = 4 * x42
x44 = Gh * x26
x45 = Bh * x44
x46 = x28 * x42
x47 = Bh * Sh
x48 = 12 * xi_y
x49 = Fyh * x36
x50 = 2 * Bp
x51 = x0 * xi_yy
x52 = Bp * x17
x53 = S * x46
x54 = x35 * xi
x55 = x28 * x54
x56 = Fy * Fyp
x57 = Fy * x32
x58 = Fy * Sh
x59 = 24 * x33
x60 = 6 * G
x61 = x26 * x60
x62 = Fyh * x32
x63 = 2 * x35
x64 = Gp * x26
x65 = Fyh * x64
x66 = 8 * S
x67 = Fyp * x44
x68 = Fyp * xi_y
x69 = S * x35
x70 = x28 * x69
x71 = 4 * xi
x72 = x23 * x71
x73 = x61 * xi_yy
x74 = x43 * x44
x75 = x44 * x48
x76 = x64 * xi_yy
x77 = S * x43
x78 = 4 * x54
x79 = Sp * x8
x80 = x7 * xi
x81 = 16 * xi_y
x82 = Sh * x81
x83 = x13 * x16
x84 = Bh ^ 2
x85 = 2 * x42
x86 = x16 * x85
x87 = Fy ^ 2
x88 = 3 * x36
x89 = Gh ^ 2
x90 = S ^ 2
x91 = x7 * x90
x92 = S ^ 3
x93 = x7 * x92
x94 = u ^ 7
x95 = 48 * x94
x96 = u ^ 10
x97 = 12 * x7
x98 = x96 * x97
x99 = u ^ 8
x100 = B * x99
x101 = x16 * x38
x102 = 12 * x101
x103 = B * x94
x104 = x16 * x48
x105 = x103 * x104
x106 = x16 * x31
x107 = x42 * x56
x108 = 12 * x99
x109 = B * Fy
x110 = x109 * x44
x111 = 12 * x16
x112 = x111 * x58
x113 = 48 * x42
x114 = x16 * xi_y
x115 = x109 * x114
x116 = S * x99
x117 = Fyh * x116
x118 = B * x111
x119 = x42 * xi
x120 = Fyh * x119
x121 = x31 * x33
x122 = x16 * xi_yy
x123 = B * S
x124 = 12 * x122 * x123
x125 = x33 * xi
x126 = x125 * xi_yy
x127 = 12 * B
x128 = 12 * G
x129 = x26 * x99
x130 = x128 * x129
x131 = 32 * x116
x132 = 24 * x119
x133 = Bh * x94
x134 = Fyp * x28
x135 = S * x134
x136 = x32 * xi
x137 = G * x26
x138 = x137 * x48
x139 = u ^ 9
x140 = S * x139
x141 = 8 * x45
x142 = x94 * xi
x143 = x28 * x47
x144 = 32 * x94
x145 = S * x16
x146 = Bh * xi_y
x147 = x145 * x146
x148 = x33 * xi_y
x149 = 2 * Fy
x150 = Bp * x149
x151 = S * x94
x152 = Fyh * x151
x153 = Bp * x28
x154 = x153 * xi
x155 = Bp * xi_yy
x156 = S * x142
x157 = x156 * x28
x158 = 8 * x56
x159 = S * x42
x160 = x159 * x16
x161 = x36 * xi
x162 = G * x16
x163 = Gh * x162
x164 = Fy * x163
x165 = Fy * xi_y
x166 = x137 * x165
x167 = Fy * x44
x168 = x32 * xi_y
x169 = x149 * x64
x170 = Fy * x114
x171 = 76 * x159
x172 = x119 * x16
x173 = x128 * x26
x174 = 4 * x151
x175 = x32 * x71
x176 = x159 * xi
x177 = x111 * x176
x178 = x13 * xi
x179 = x32 * x61
x180 = S * x32
x181 = x180 * x28
x182 = x68 * xi
x183 = x181 * xi
x184 = x163 * x48
x185 = x173 * xi_yy
x186 = Sh * x94
x187 = 4 * x140
x188 = Sh * x44
x189 = x44 * xi_y
x190 = S * x189
x191 = x71 * x94
x192 = 24 * xi
x193 = x151 * x71
x194 = Sp * x20
x195 = x42 * x83
x196 = S * xi_yy
x197 = 24 * B
x198 = x16 * x197
x199 = x87 * x94
x200 = xi_y ^ 2
x201 = x200 * x32
x202 = x28 * x84
x203 = Bp * x86
x204 = x139 * x90
x205 = xi ^ 2
x206 = x205 * x32
x207 = 24 * x137
x208 = x64 * x87
x209 = 54 * x145
x210 = 6 * x125
x211 = x90 * x99
x212 = x16 * x211
x213 = 10 * x212
x214 = x205 * x36
x215 = 2 * x214
x216 = x17 * x180
x217 = x0 * x17
x218 = x30 * xi
x219 = x90 * x94
x220 = Fyph * x17
x221 = x0 * x205
x222 = x28 * x89
x223 = 2 * x27
x224 = x205 * x7
x225 = xi ^ 3
x226 = x225 * x7
x227 = xi ^ 4
x228 = x180 * x8
x229 = x14 * x42
x230 = x32 * x90
x231 = x99 * xi
x232 = x194 * x205
x233 = x14 * x35
x234 = u ^ 11
x235 = S * x234
x236 = x198 * x38
x237 = x139 * xi
x238 = S * x96
x239 = x146 * x198
x240 = x142 * x56
x241 = 72 * x139
x242 = x109 * xi_y
x243 = x137 * x242
x244 = x167 * x197
x245 = x145 * x234
x246 = 156 * x140
x247 = B * xi
x248 = 108 * x115
x249 = x140 * xi
x250 = Fyh * x249
x251 = B * x104
x252 = Fyp * x251
x253 = x189 * x197
x254 = B * x96
x255 = x145 * x254
x256 = Sh * x231
x257 = x207 * x38
x258 = 36 * x249
x259 = x146 * x207
x260 = x96 * xi
x261 = x260 * x66
x262 = 32 * x231
x263 = x116 * x165
x264 = x165 * xi
x265 = S * x231
x266 = x173 * x56
x267 = x16 * x66
x268 = 24 * x235
x269 = x164 * x192
x270 = x173 * x58
x271 = 108 * x166
x272 = 4 * Fy
x273 = x142 * x145
x274 = 4 * x265
x275 = Fyp * x138
x276 = x163 * xi_y
x277 = 24 * x276
x278 = x26 * x7
x279 = 2 * x278
x280 = u * x279
x281 = x137 * x87
x282 = 36 * x281
x283 = x16 * x87
x284 = x283 * x96
x285 = B * x283
x286 = 60 * x231
x287 = x234 * x90
x288 = Fyh * x287
x289 = x205 * x94
x290 = Fyh * x289
x291 = x137 * x200
x292 = 36 * x291
x293 = x16 * x200
x294 = B * x293
x295 = 72 * x116
x296 = x90 * x96
x297 = x122 * x31
x298 = x113 * xi
x299 = x205 * x42
x300 = B ^ 2
x301 = 36 * x300
x302 = x139 * x170
x303 = 22 * x287
x304 = 14 * x289
x305 = Bh * x296
x306 = Fyp * x17
x307 = Bh * x86
x308 = u ^ 12
x309 = x308 * x90
x310 = 4 * x45
x311 = x205 * x99
x312 = 20 * x114
x313 = Bh * x104
x314 = x238 * xi
x315 = x140 * x87
x316 = x199 * xi
x317 = x296 * x52
x318 = x203 * x205
x319 = x52 * xi_yy
x320 = x28 * x56
x321 = G ^ 2
x322 = 36 * x321
x323 = x206 * x28
x324 = 84 * x238
x325 = x64 * x71
x326 = 2 * x296
x327 = x205 * x85
x328 = S * x86
x329 = 2 * x76
x330 = x159 * x79
x331 = x308 * xi
x332 = S * x331
x333 = x235 * xi
x334 = 144 * x123 * x166
x335 = x145 * x260
x336 = x140 * x165
x337 = S * x308
x338 = 18 * x300
x339 = x293 * x99
x340 = x17 * x84
x341 = exp(2 * x6)
x342 = Fxpt * x341
x343 = 18 * x321
x344 = 21 * x90
x345 = 5 * x299
x346 = x17 * x89
x347 = x91 * x94
x348 = 16 * x347
x349 = x139 * x205
x350 = u ^ 14
x351 = x350 * x90
x352 = x205 * x254
x353 = u ^ 13
x354 = x353 * x90
x355 = B * x313
x356 = x106 * x56
x357 = 12 * x110
x358 = x205 * x96
x359 = 60 * x311
x360 = 72 * x123
x361 = x281 * x360
x362 = x234 * x281
x363 = 72 * x247
x364 = x234 * x283
x365 = 96 * xi
x366 = x106 * x68
x367 = x291 * x360
x368 = B * x75
x369 = 72 * x249
x370 = 72 * x300
x371 = x145 * x308
x372 = x165 * x371
x373 = x170 * x260
x374 = x173 * x38
x375 = Bh * x138
x376 = x114 * x150
x377 = x154 * x238
x378 = x56 * x61
x379 = 12 * x164
x380 = 72 * x321
x381 = x169 * xi_y
x382 = x61 * x68
x383 = x279 * x5
x384 = x29 * x7
x385 = x0 * x26
x386 = Sc * x8
x387 = 60 * x354
x388 = 36 * x349
x389 = 48 * x287
x390 = x200 * x289
x391 = S * x301
x392 = x283 * x353
x393 = x301 * xi
x394 = x234 * x293
x395 = x139 * x293
x396 = x52 * x87
x397 = S * x322
x398 = x322 * xi
x399 = 2 * x208
x400 = u ^ 15
x401 = x400 * x90
x402 = 72 * x243
x403 = x205 * x234
x404 = x350 * xi
x405 = x145 * x264 * x353
x406 = Bb * x341
x407 = Fxt * x341
x408 = x26 * x341
x409 = Gb * x408
x410 = x29 * x341
x411 = Gt * x7
x412 = Gh * x7
x413 = Bt * x412
x414 = x278 * x35
x415 = Fx * Fy
x416 = Fx * x7
x417 = 4 * Fyp
x418 = x385 * x417
x419 = x39 * x7
x420 = Fx * Gh
x421 = x20 * x26
x422 = Fx * x32
x423 = x60 * x7
x424 = x33 * x423
x425 = Fxh * Gp
x426 = 2 * x7
x427 = x36 * x426
x428 = x228 * x26
x429 = x233 * x26
x430 = Fxp * x272
x431 = Fxp * x278
x432 = Fxp * x36
x433 = Fxp * xi_y
x434 = x278 * x69
x435 = 4 * x434
x436 = Fxph * x26
x437 = 4 * x80
x438 = x437 * x5
x439 = Fy * Gt
x440 = Fypt * x26
x441 = Fyt * Gp
x442 = x36 * x97
x443 = Gc * x8
x444 = St * x412
x445 = Gp * xi_xy
x446 = Sh * x411
x447 = Gt * xi_y
x448 = x26 * x386
x449 = 16 * xi_xy
x450 = 8 * Sh
x451 = St * x42
x452 = St * x81
x453 = u ^ 16
x454 = x453 * x90
x455 = B * x282
x456 = x205 * x308
x457 = x170 * x301
x458 = x391 * xi
x459 = x283 * x350
x460 = x293 * x308
x461 = x170 * x322
x462 = x397 * xi
x463 = Fxp ^ 2 * x341
x464 = x283 * x338
x465 = x293 * x351
x466 = x338 * x96
x467 = x205 * x293
x468 = x36 * x407
x469 = Fx * x341
x470 = Bt * x469
x471 = 2 * x341
x472 = Gt * x408
x473 = x43 * x472
x474 = x341 * x46
x475 = Bt * St
x476 = Fx * Fxp
x477 = Fxp * x472
x478 = x407 * x64
x479 = x283 * x343
x480 = x343 * x96
x481 = Sb * x341
x482 = x106 * x420
x483 = x7 * x99
x484 = x106 * x439
x485 = x114 * x94
x486 = x16 * x99
x487 = Fx * x423
x488 = x140 * x28
x489 = x28 * x80
x490 = Gt * x489
x491 = Bt * x423
x492 = Fy * x486
x493 = Bt * Gh
x494 = x489 * x493
x495 = x162 * x416
x496 = Fy * x495
x497 = Gp * x46
x498 = x278 * x415
x499 = x26 * x415
x500 = x499 * x79
x501 = Fx * xi
x502 = x26 * x97
x503 = x16 * x42
x504 = Fx * Fyp
x505 = x26 * x8
x506 = x159 * x505
x507 = x505 * x54
x508 = x97 * x99
x509 = G * x44
x510 = Fx * x509
x511 = Fx * Sh
x512 = x162 * x508
x513 = x495 * xi_y
x514 = Fx * x116
x515 = 32 * x514
x516 = x16 * x515
x517 = Gh * x16
x518 = x42 * x501
x519 = 2 * Fx
x520 = Gp * x519
x521 = Fx * xi_y
x522 = x278 * x521
x523 = x26 * x80
x524 = 16 * x523
x525 = x14 * x26
x526 = x525 * xi_y
x527 = x26 * x521
x528 = x162 * x97
x529 = Fxh * x528
x530 = x151 * x28
x531 = x530 * x7
x532 = x32 * x489
x533 = Fxh * x525
x534 = x176 * x502
x535 = Fxp * Fy
x536 = x180 * x431
x537 = Gh * x489
x538 = Fxp * x537
x539 = 4 * xi_y
x540 = x180 * x437
x541 = x137 * x439
x542 = Fy * St
x543 = x16 * x439
x544 = x20 * x439
x545 = Gt * x134 * x80
x546 = Fyt * x528
x547 = Fyt * x525
x548 = x137 * x447 * x97
x549 = x162 * xi_xy
x550 = x20 * x549
x551 = St * x94
x552 = x528 * xi_y
x553 = Gt * x44
x554 = x140 * x8
x555 = x553 * x8
x556 = x445 * x8
x557 = x145 * x447
x558 = x20 * x447
x559 = x26 * xi_xy
x560 = Bt ^ 2
x561 = x341 * x86
x562 = Fx ^ 2 * x341
x563 = Gt ^ 2
x564 = Fxh * x91
x565 = 10 * x129
x566 = x224 * x26
x567 = x566 * x63
x568 = 2 * x347
x569 = x25 * x566
x570 = Fyt * x91
x571 = Gc * x28
x572 = x139 * x91
x573 = x224 * x32
x574 = B * Fx
x575 = x108 * x574
x576 = x16 * x341
x577 = Bt * x576
x578 = x42 * x476
x579 = x106 * x341
x580 = x118 * x407
x581 = x153 * x407
x582 = x128 * x408
x583 = Fx * x582
x584 = x583 * x99
x585 = Bt * x341
x586 = 24 * x518
x587 = x28 * x341
x588 = Bt * Fxp
x589 = x587 * x588
x590 = Bt * x472
x591 = 8 * x590
x592 = x140 * x587
x593 = x142 * x587
x594 = x408 * x60
x595 = 8 * x341
x596 = Gt * x162
x597 = x469 * x596
x598 = x173 * x407
x599 = St * x472
x600 = x574 * x97
x601 = Gh * x600
x602 = Gt * x109 * x97
x603 = Gt * x114
x604 = x235 * x528
x605 = x139 * x501
x606 = x528 * x605
x607 = Bt * Fy
x608 = x237 * x528
x609 = Bt * S
x610 = x162 * x98
x611 = x610 * xi_y
x612 = x231 * x552
x613 = 168 * x238
x614 = x162 * x80
x615 = 120 * x415
x616 = Gp * x415
x617 = x501 * x94
x618 = Fy * x8
x619 = Fy * x26
x620 = x140 * x528
x621 = Fyp * x501
x622 = x20 * x509
x623 = 108 * x521
x624 = 36 * x140 * x80
x625 = Gp * x7
x626 = x521 * x80
x627 = x527 * x80
x628 = x116 * x489
x629 = Fxp * Fyp
x630 = x523 * x77
x631 = x433 * x528
x632 = x137 * x544
x633 = x137 * x558
x634 = x562 * x94
x635 = x17 * x406
x636 = x562 * x64
x637 = x463 * xi
x638 = x17 * x342
x639 = 2 * x409
x640 = x308 * x91
x641 = Bh * Gt * x17
x642 = x224 * x99
x643 = x17 * x493
x644 = x91 * x96
x645 = x415 * x566
x646 = x26 * x572
x647 = x422 * x566
x648 = x16 * x234
x649 = x420 * x91
x650 = 14 * Fx
x651 = x224 * x94
x652 = x60 * x648
x653 = x16 * x60
x654 = x651 * x653
x655 = Gp * x96
x656 = x17 * x655
x657 = x224 * x86
x658 = x129 * x91
x659 = Fxp * Gh
x660 = x234 * x91
x661 = Gt * x644
x662 = x299 * x97
x663 = 4 * x553
x664 = x28 * x445
x665 = x198 * x585
x666 = Fxp * x341
x667 = x111 * x574
x668 = x666 * x667
x669 = Fx * x472
x670 = x197 * x472
x671 = St * x341
x672 = 12 * x574
x673 = G * x408
x674 = Bt * x673
x675 = 24 * x605
x676 = 36 * x501
x677 = x140 * x676
x678 = Fxp * x582
x679 = x341 * x596
x680 = x371 * xi
x681 = x337 * x501
x682 = S * x610
x683 = x137 * x562
x684 = 36 * x683
x685 = x16 * x562
x686 = B * x685
x687 = x106 * x407
x688 = x153 * x562
x689 = x289 * x650
x690 = x17 * x341
x691 = 4 * x590
x692 = x314 * x587
x693 = x562 * x83
x694 = x145 * x562
x695 = x407 * x61
x696 = x224 * x96
x697 = x350 * x91
x698 = x353 * x91
x699 = x31 * x603
x700 = x139 * x224
x701 = Fx * x697
x702 = Bh * x653
x703 = x607 * x653
x704 = x114 * x60
x705 = Bt * x704
x706 = S * x380
x707 = x353 * x706
x708 = x380 * x80
x709 = x28 * x616
x710 = x640 * x653
x711 = x224 * x60
x712 = x205 * x98
x713 = x205 * x519
x714 = Fxp * x704
x715 = Gt * x128
x716 = x560 * x690
x717 = x563 * x690
x718 = x205 * x574
x719 = x718 * x96
x720 = x476 * x579
x721 = 12 * x472
x722 = x360 * x683
x723 = Bt * x583
x724 = x476 * x594
x725 = 12 * x597
x726 = x353 * x694
x727 = x234 * x685
x728 = x52 * x562
x729 = 2 * x636
x730 = x322 * x91
x731 = B * x684
x732 = x350 * x694
x733 = x338 * x685
x734 = x343 * x685
x735 = 3 * G
x736 = x26 * x4
x737 = 3 * x6
x738 = x4 * x737
x739 = 36 * x1 * x3 + 12 * x3 ^ 2 + 12 * x85 * x90
x740 = 9 * x0
x741 = x321 * x42
x742 = 9 * x741
x743 = 6 * S * x3 * (3 * x741 + 2) + x0 * x90 * (x742 + 4) + x740 * (G * x2 + G) ^ 2
x744 = 6 * x15
x745 = x32 * x321
x746 = 10 * x1 + 6 * x2 + 2 * x738 + 6
x747 = Gp * x4
x748 = x4 * x735
x749 = 6 * x4 * x6
x750 = 9 * x1 + 5 * x2 + 4
x751 = u * x60
x752 = x519 * xi
x753 = Bp * Fx
x754 = Fx * x5
x755 = Bp * x752
x756 = 2 * Sp
x757 = 18 * Fx
x758 = x300 * x757
x759 = x301 * x501
x760 = Gt * x60
x761 = x321 * x757
ABCS[1] = 0
ABCS[2] = -x4 ^ 3 * x9
ABCS[3] = -x10 * x9 * (-x11 + 3 * x12 + xi)
ABCS[4] = B * x231 * x603 * x97 + B * x245 * x365 * x562 - B * x292 * x351 + B * x557 * x98 + Bh * Fx * x604 + Bh * x105 + Bh * x134 * x136 + Bh * x135 * x231 + Bh * x238 * x490 + Bh * x411 * x488 + Bh * x486 * x487 + Bh * x528 * x681 + Bh * x606 + Bp * x264 * x46 - Bs * x157 - Bs * x53 - Bs * x55 - Bt * x333 * x552 - Bt * x432 * x471 + Bt * x473 + Bt * x584 - Bt * x612 + Fx * x268 * x674 + Fx * x417 * x646 + Fx * x696 * x702 + Fxh * x383 - Fxh * x424 - Fxh * x428 + Fxh * x429 - Fxh * x534 - Fxh * x567 - Fxh * x654 - Fxp * x41 * x658 + Fxp * x412 * x530 - Fxp * x418 * x80 + Fxp * x492 * x711 + 4 * Fxp * x566 * x57 + Fxph * x280 + Fxph * x435 + Fxph * x569 - Fy * x13 * x148 + Fy * x416 * x497 - Fy * x451 * x524 - 84 * Fy * x514 * x523 + Fyh * x153 * x265 + Fyh * x177 + Fyh * x213 + Fyh * x215 - Fyh * x24 + Fyh * x317 + Fyh * x318 + Fyp * Gt * x657 + Fyp * x104 * x140 * x247 + Fyp * x205 * x307 + Fyp * x487 * x503 + Fyp * x528 * x617 - Fyph * x18 - Fyph * x183 - Fyph * x70 - Fyph * x72 + Fypt * x280 + Fypt * x435 + Fypt * x569 + Fyt * x383 - Fyt * x424 - Fyt * x428 + Fyt * x429 - Fyt * x534 - Fyt * x567 - Fyt * x654 - G * x442 * xi_xy - Gc * x384 + Gh * x426 * x432 + Gp * x16 * x617 * x618 + Gp * x488 * x626 - Gt * x138 * x698 - S ^ 4 * x98 - S * St * x611 + S * x362 * x365 + Sb * x410 - Sh ^ 2 * x46 - Sh * x105 + Sh * x138 * x238 - Sh * x255 * x48 - Sh * x421 * x422 - Sh * x606 + Sh * x74 + 16 * Sp * x0 * x226 + Sp * x21 * x7 + Ss * x29 + Ss * x53 + Ss * x55 - St ^ 2 * x474 + St * x235 * x583 - St * x421 * x57 + St * x469 * x59 + St * x473 + St * x575 * x576 + St * x582 * x605 + St * x584 - St * x612 - x0 * x22 - 96 * x1 * x224 + x1 * x79 + x100 * x102 - x100 * x112 - x100 * x292 + x101 * x127 * x351 - x101 * x131 - x101 * x132 - x101 * x258 - x101 * x303 - x101 * x304 + x102 * x352 - x103 * x75 + x106 * x107 - x106 * x288 - x106 * x290 - x106 * x350 * x649 - x107 * x61 - x108 * x110 + x108 * x164 + x108 * x597 + x111 * x585 * x719 - x112 * x139 * x247 - x113 * x115 + x113 * x166 - x113 * x513 + x114 * x520 * x660 - x115 * x246 - x115 * x359 + x116 * x252 - x116 * x272 * x64 * xi_y - x116 * x275 - x116 * x529 + x116 * x538 + x116 * x545 - x116 * x546 + x116 * x580 + x116 * x598 + x116 * x631 - x117 * x118 + x117 * x173 - x118 * x120 + x118 * x140 * x56 + x118 * x240 - x118 * x250 + x119 * x252 - x119 * x275 - x119 * x529 - x119 * x546 + x119 * x580 + x119 * x598 + x119 * x631 - 80 * x12 * x80 + x120 * x173 + x121 * x407 + x121 * x68 + 8 * x122 * x219 - 84 * x123 * x284 - x123 * x364 * x365 - x124 * x231 - x124 * x94 - x125 * x13 * x407 - 24 * x125 * x146 - x125 * x558 + x125 * x82 - x126 * x127 + x126 * x66 - x127 * x245 * x58 - x128 * x44 * x701 - x13 * x468 - x13 * x49 - x130 * x38 + x130 * x58 + x131 * x167 - x131 * x543 * x7 - x131 * x603 * x80 + x132 * x167 + x133 * x135 - x133 * x138 + x133 * x490 + x134 * x151 * x411 + x136 * x185 + x136 * x533 + x136 * x547 - x136 * x550 - x136 * x581 - x136 * x589 + x138 * x186 + x138 * x256 + x139 * x269 - x139 * x322 * x522 + x14 * x205 * x211 + x14 - x140 * x141 - x140 * x143 + x140 * x202 + x140 * x222 - x140 * x266 - x140 * x476 * x582 + x140 * x500 + x140 * x591 - x140 * x668 - x140 * x688 - x140 * x693 - x141 * x142 - x142 * x143 + x142 * x202 + x142 * x222 - x142 * x248 - x142 * x266 + x142 * x271 + x142 * x528 * x535 - x142 * x555 + x142 * x591 - x142 * x668 - x142 * x688 - x142 * x693 - x144 * x147 + x144 * x190 - x144 * x557 * x7 + 24 * x145 * x201 + x145 * x501 * x618 * x655 - x147 * x262 + x148 * x150 + x148 * x520 * x7 + x151 * x185 - x151 * x407 * x83 - 108 * x151 * x498 + x151 * x505 * x621 + x151 * x533 + x151 * x547 - x151 * x550 - x151 * x581 - x151 * x589 - 64 * x151 * x627 + x152 * x153 - x152 * x83 + x153 * x263 + x153 * x315 + x153 * x316 + x154 * x336 + x154 * x62 + x155 * x157 + x155 * x53 + x155 * x55 + x156 * x505 * x535 + x157 * x406 - x158 * x160 - x158 * x161 - x159 * x448 + x16 * x237 * x602 + x16 * x303 * x470 - x16 * x420 * x624 + x16 * x554 * x616 - x160 * x443 - x160 * x476 * x595 + x160 * x556 - x161 * x443 + x161 * x556 - x162 * x224 * x241 * x415 - 60 * x162 * x521 * x642 - x162 * x615 * x698 - x162 * x623 * x640 + x164 * x268 + 64 * x165 * x273 + x165 * x323 + x165 * x55 + x166 * x246 + 168 * x166 * x314 + x166 * x359 + x167 * x258 + x167 * x303 + x167 * x304 - x168 * x169 + x168 * x192 * x44 - x17 * x211 * x68 + x17 * x644 * x659 + x170 * x171 + 28 * x170 * x204 - x171 * x522 - x172 * x544 + 16 * x172 * x58 + x173 * x196 * x231 + x173 * x250 - x174 * x477 - x174 * x478 - x174 * x65 - x174 * x67 - x175 * x477 - x175 * x478 - x175 * x65 - x175 * x67 + 24 * x176 * x293 + x177 * x407 - x178 * x34 - x178 * x37 + x179 * x407 - x179 * x68 - x18 * x342 + x180 * x232 - x180 * x449 * x523 - x181 * x68 - x182 * x29 - x182 * x53 - x183 * x342 + x184 * x349 + x184 * x354 + x184 * x94 - x186 * x490 + x187 * x188 - x187 * x208 + x187 * x599 - x187 * x636 + x188 * x191 + 20 * x189 * x296 - x19 * x20 - x190 * x197 * x234 * xi + x190 * x262 + x191 * x599 + x192 * x235 * x276 + x193 * x27 + x193 * x409 - x193 * x76 + x194 * x69 * xi - x195 * x196 - x195 * x264 - x195 * x562 - x195 * x87 + x197 * x235 * x669 - x198 * x199 - x198 * x201 + x198 * x235 * x470 - x198 * x390 + x198 * x634 + x199 * x207 + x199 * x209 - x199 * x325 - x20 * x235 * x510 - x20 * x517 * x518 + 8 * x200 * x212 + x201 * x207 + x202 * x314 - x203 * x562 + x203 * x87 - x204 * x22 + x204 * x223 + x204 * x319 - x204 * x320 - x204 * x329 - x204 * x476 * x587 + x204 * x635 + x204 * x639 - x205 * x561 * x588 - x206 * x22 + x206 * x223 + x206 * x319 - x206 * x320 - x206 * x329 + x206 * x635 + x206 * x639 + x207 * x390 + x207 * x634 - x208 * x238 * x71 - x208 * x85 + x209 * x634 + x210 * x562 + x210 * x87 + x212 * x30 + x212 * x463 + x213 * x407 + x214 * x30 + x214 * x463 + x215 * x407 - x215 * x68 + x216 * x30 + x216 * x463 + x217 * x218 + x217 * x406 + x217 * x637 + x218 * x328 - x219 * x220 - x219 * x638 + x219 * x79 * xi - x220 * x221 - x221 * x638 + x222 * x314 + x225 * x330 - x225 * x348 - 48 * x226 * x69 - x227 * x228 + x227 * x233 + x229 * x26 * x501 * xi_y + x229 * x90 + x23 * x30 + x23 * x463 - 136 * x230 * x80 + x231 * x239 - x231 * x253 - x231 * x259 + x231 * x277 - x231 * x633 + 42 * x231 * x694 - 56 * x231 * x93 + x232 * x5 + x234 * x247 * x557 * x97 - x234 * x322 * x521 * x566 + x234 * x363 * x683 - x234 * x367 - x234 * x499 * x708 + x235 * x236 - x235 * x244 - x235 * x257 + x235 * x270 + x235 * x365 * x683 - 192 * x235 * x415 * x614 - x235 * x632 + x236 * x237 + x236 * x332 - x237 * x244 - x237 * x257 + x237 * x270 - x237 * x517 * x600 - x237 * x632 + x237 * x667 * x671 + x238 * x239 - x238 * x253 - x238 * x259 + x238 * x277 - x238 * x325 * x562 - x238 * x494 - x238 * x501 * x678 - x238 * x633 + x239 * x333 - x24 * x407 - x24 * x68 - x240 * x267 - x241 * x243 - x241 * x247 * x291 - 168 * x242 * x335 - 144 * x243 * x260 - x244 * x332 - x245 * x601 + x245 * x602 + x245 * x671 * x672 - x246 * x513 - x248 * x309 - x249 * x275 - x249 * x529 - x249 * x546 + x249 * x580 + x249 * x598 + x249 * x631 + x25 * x27 - x25 * x341 * x4 * xi_xx * (x16 * (-3 * B * u * x4 + Bp * x4 + 2 * Sp - 4 * x19) + x736 * (Gp - u * x735)) + x25 * x409 - x25 * x76 - x251 * x256 - x254 * x282 + x254 * x684 + 12 * x255 * x56 * xi + 84 * x255 * x562 - x257 * x332 - x259 * x333 - x261 * x45 + x261 * x590 - x263 * x83 - x264 * x43 * x64 + 42 * x265 * x283 - x265 * x550 - x265 * x581 - x265 * x589 - x266 * x314 - x267 * x617 * x666 + x268 * x597 + x269 * x337 + x27 * x77 + x27 * x78 + x271 * x309 - x273 * x443 + x273 * x556 - x274 * x477 - x274 * x478 - x274 * x65 - x274 * x67 + x278 * x450 * x451 + x28 * x514 * x625 * xi_y + x281 * x286 + x281 * x324 + x281 * x387 + x282 * x349 + x283 * x345 + x284 * x338 + x284 * x343 + x284 * x344 - x285 * x286 - x285 * x387 - x285 * x388 + x286 * x683 + x286 * x686 + x287 * x366 + x287 * x376 - x287 * x381 - x287 * x382 + x287 * x687 + x287 * x695 + x288 * x61 + x289 * x366 + x289 * x376 - x289 * x381 - x289 * x382 + x289 * x687 + x289 * x695 - x29 * x56 + x290 * x61 + x291 * x295 + x291 * x298 + x291 * x369 + x291 * x389 - x292 * x352 - x294 * x295 - x294 * x298 - x294 * x369 - x294 * x389 - x296 * x297 - x296 * x588 * x690 + x296 * x73 - x297 * x299 - x299 * x313 + x299 * x73 + x299 * x75 - 92 * x299 * x91 + x301 * x302 + x301 * x726 + x302 * x322 + x303 * x669 + x305 * x306 - x305 * x312 + x306 * x661 + x307 * x411 - x308 * x322 * x645 - x308 * x334 + 24 * x308 * x501 * x609 * x673 - x308 * x522 * x706 - x309 * x310 + x309 * x340 + x309 * x346 + x309 * x356 - x309 * x378 + x309 * x396 - x309 * x399 + x309 * x691 + x309 * x716 + x309 * x717 - x309 * x720 - x309 * x724 - x309 * x728 - x309 * x729 - x31 * x34 - x31 * x37 + x31 * x411 * x485 - x310 * x311 + x311 * x340 + x311 * x346 + x311 * x356 - x311 * x378 + x311 * x396 - x311 * x399 + x311 * x691 + x311 * x716 + x311 * x717 - x311 * x720 - x311 * x724 - x311 * x728 - x311 * x729 - x312 * x661 - x314 * x555 - x315 * x83 - x316 * x83 - x317 * x407 - x318 * x407 - x32 * x452 * x523 + x32 * x471 * xi_x ^ 2 * (x16 * (B * x739 + x10 * x300 * x740 + x743) + x60 * x736 * (4 * x1 + 2 * x2 + x738 + 2)) + x32 * x538 + x32 * x545 - x322 * x498 * x96 + x322 * x726 - x323 * x341 * x476 + x324 * x683 - x325 * x336 - x325 * x634 - x326 * x477 - x326 * x478 - x326 * x65 - x326 * x67 - x327 * x477 - x327 * x478 - x327 * x65 - x327 * x67 + x328 * x637 + x33 * x407 * x66 + x330 * x559 - x331 * x367 - x332 * x528 * x607 - x332 * x632 - x333 * x633 - x334 * x353 * xi - x335 * x666 * x672 + x338 * x339 + x338 * x465 + x339 * x343 + x34 * x66 - x342 * x70 - x342 * x72 + x343 * x465 + x344 * x685 * x96 + x345 * x685 - x348 * x559 + x349 * x355 - x349 * x368 - x349 * x375 - x349 * x548 - x349 * x8 * x92 - x35 * x437 * x527 + x35 * x73 + x35 * x75 - 60 * x35 * x91 - x350 * x499 * x706 * x80 - x351 * x357 - x351 * x374 + x351 * x379 + x351 * x574 * x721 + x351 * x585 * x667 + x351 * x723 + x351 * x725 - x353 * x361 + x353 * x722 + x354 * x355 - x354 * x368 - x354 * x375 - x357 * x358 - x358 * x374 + x358 * x379 + x358 * x723 + x358 * x725 + x36 * x41 * x411 + x36 * x82 - x361 * x404 - x362 * x363 + x364 * x393 + x364 * x398 + x37 * x66 + x370 * x372 + x370 * x373 + x370 * x405 + x372 * x380 + x373 * x380 - x377 * x562 + x377 * x87 - x38 * x39 + x380 * x405 + x383 * x433 + x384 * x445 - x385 * x386 + x385 * x430 * x7 + x385 * x433 * x437 + x385 * x79 * xi_xy + x387 * x683 + x387 * x686 + x388 * x683 + x388 * x686 + x39 * x470 + x391 * x392 + x391 * x394 + x392 * x397 + x393 * x395 + x393 * x727 + x393 * x732 + x394 * x397 + x395 * x398 + x398 * x727 + x398 * x732 + x40 * x41 - x40 * x48 - x400 * x527 * x730 - x401 * x402 + x401 * x457 + x401 * x461 - x402 * x403 + x403 * x457 + x403 * x461 + x404 * x722 + x406 * x53 + x406 * x55 + x409 * x77 + x409 * x78 - x41 * x431 * x5 - x410 * x476 - x411 * x74 - x412 * x516 - x413 * x488 - x413 * x86 + 6 * x414 * x415 - x414 * x452 + x416 * x418 - x417 * x536 + x417 * x647 - x419 * x420 - x419 * x439 + x42 * x500 - x42 * x511 * x524 - 10 * x42 * x645 + 10 * x422 * x472 + x422 * x526 + x423 * x503 * x535 + x424 * x433 + x425 * x427 + x425 * x531 + x425 * x532 + x425 * x628 + x425 * x657 + x426 * x5 * xi_x * (x16 * (-x4 * x5 * (Gh * (10 * x1 + 6 * x2 + x738 + 6) - u * (Fy * x747 - Fy * x750 * x751 + Fyp * x748 + x5 * x735 * (Bh * x1 + Bh * x2 + Bh - 2 * Sh))) + x7 * (-Fxp * x10 * (x737 + 1) + x5 * (-Bp * x519 * x69 + Bt * x4 * x746 + 38 * Fx * x12 - Fx * x69 * x756 + Gt * x744 + St * (8 * x2 + x749 + 8) + u * x713 - u * x753 + 32 * x1 * x501 - x11 * x519 + x140 * x322 * x501 + x140 * x759 + x156 * x715 + x159 * x715 + 84 * x176 * x574 + 78 * x180 * x574 - x180 * x755 + x197 * x754 + x204 * x760 + x206 * x760 + 54 * x211 * x574 - x219 * x753 - x221 * x753 + x230 * x650 + x287 * x758 + x287 * x761 + x289 * x758 + x289 * x761 + x301 * x514 + x32 * x758 + x322 * x514 + 30 * x35 * x718 + x42 * x759 - x5 * x755 + 54 * x501 * x6 + x54 * x715 + x676 * x741 + x745 * x757 + x752 - x754 * x756 * xi))) - x25 * xi_y * (x162 * x739 + x26 * x743) + x26 * (Fyp * x10 + x5 * (-Gh * x10 * x744 - x149 * (-x11 * x4 + x12 * (x21 + 18 * x3 * x741 + 19) + x230 * (x742 + 7) + x3 * (9 * x3 * x745 + xi)) - x3 * x450 + x4 * x7 * (Gt * x746 + u * (-Fx * x747 + Fx * x751 * (x749 + x750) - Fxp * x748 + x5 * x60 * (Bt * x1 + Bt * x2 + Bt + St)))))) + x427 * x441 - x43 * x45 + x430 * x646 - x432 * x501 * x595 + x433 * x567 + x433 * x630 + 2 * x433 * x658 - x434 * x449 + x436 * x438 + x436 * x540 + x436 * x568 + x438 * x440 + 10 * x44 * x57 + x440 * x540 + x440 * x568 + x441 * x531 + x441 * x532 + x441 * x628 + x441 * x657 - x442 * x447 - x444 * x46 - x444 * x488 - x446 * x46 - x446 * x488 - x448 * x54 + 16 * x451 * x501 * x576 - x453 * x499 * x730 - x454 * x455 + x454 * x464 + x454 * x479 + x454 * x731 + x454 * x733 + x454 * x734 - x455 * x456 + x456 * x464 + x456 * x479 + x456 * x731 + x456 * x733 + x456 * x734 + x458 * x459 + x458 * x460 + x459 * x462 - x46 * x47 + x460 * x462 + x466 * x467 + x466 * x685 + x467 * x480 - x468 * x50 + x472 * x515 + x472 * x575 + x472 * x586 + x472 * x677 + x472 * x689 + x474 * x475 + x475 * x592 + x475 * x593 - x477 * x63 - x478 * x63 + x480 * x685 + x481 * x53 + x481 * x55 - x482 * x483 - x482 * x696 + x483 * x484 + x484 * x696 + x484 * x697 - x485 * x491 + x485 * x625 * x713 + x486 * x504 * x711 + x49 * x50 - x491 * x492 - x494 * x94 - x496 * x613 - x496 * x95 + x497 * x626 - x498 * x707 - 42 * x499 * x644 - x501 * x502 * x57 + x504 * x506 + x504 * x507 + x504 * x620 + x504 * x710 + x506 * x535 + x507 * x535 - x508 * x510 - x508 * x541 + x51 * x52 - x51 * x83 - x510 * x712 - x511 * x512 - x511 * x604 - x512 * x542 + x514 * x526 + x516 * x585 - x517 * x650 * x651 - x521 * x613 * x614 - 28 * x527 * x572 - x527 * x708 * x96 + x535 * x620 + x535 * x682 * xi + x535 * x710 + x536 * x539 - x537 * x551 - x539 * x647 + x54 * x559 * x79 - x541 * x712 - x542 * x604 - x542 * x608 - x543 * x624 - 14 * x543 * x651 - 22 * x543 * x660 - x548 * x94 - 12 * x549 * x644 - x549 * x662 - x551 * x552 - x553 * x554 + x560 * x561 + x560 * x592 + x560 * x593 + x560 * x692 + x561 * x563 - x562 * x88 + x563 * x592 + x563 * x593 + x563 * x692 - x564 * x565 - x564 * x652 + x564 * x656 - x565 * x570 - x567 * x629 - x570 * x652 + x570 * x656 - x571 * x572 - x571 * x573 + x572 * x664 + x573 * x664 + x575 * x577 + x577 * x586 + x577 * x677 + x577 * x689 - x578 * x579 - x578 * x594 + x58 * x59 - x601 * x680 + x602 * x680 - x603 * x662 - x604 * x607 - x605 * x622 + x605 * x665 + x605 * x670 - x607 * x608 - x609 * x611 + x61 * x62 - x614 * x615 * x99 - x614 * x623 * x94 + x617 * x619 * x79 - x617 * x678 - x619 * x697 * x715 + x621 * x682 - x622 * x681 - x627 * x707 - x629 * x630 - x63 * x65 - x63 * x67 - x636 * x85 + x640 * x641 - x640 * x643 - x640 * x663 + x640 * x709 + x641 * x642 - x642 * x643 - x642 * x663 + x642 * x709 - 22 * x648 * x649 + x651 * x714 + x657 * x659 + x660 * x714 + x665 * x681 + x670 * x681 + x674 * x675 + x675 * x679 + 24 * x679 * x681 - x696 * x703 - x697 * x703 + x698 * x699 - x698 * x705 + x699 * x700 - x700 * x705 + x701 * x702 + x719 * x721 - x76 * x77 - x76 * x78 + x84 * x86 + x86 * x89 - x87 * x88 - x93 * x95

	
    nothing
end




# this is another coupled equation, for Bd and Gd. the notation used is
#
# ( A11 d_uu Bd + A12 d_uu Gd + B d_u Bd + B2 d_u Gd + C11 Bd + C12 Gd ) = -S1
# ( A21 d_uu Bd + A22 d_uu Gd + B21 d_u Bd + B22 d_u Gd + C21 Bd + C22 Gd ) = -S2

function BdGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,       G      ,        S      ,    Fx     ,    Fy     ,  Sd,
        Bp    ,       Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,       Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,       G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,       G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,       Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,       Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,       G_xx   ,        S_xx   ,
        B_yy  ,       G_yy   ,        S_yy   ,
                       G_xy   ,        S_xy
    ) = vars

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")


x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = 8 * x2
x4 = u ^ 2
x5 = S * x0
x6 = u * xi
x7 = x6 + 1
x8 = x5 + x7
x9 = x8 ^ 3
x10 = x4 * x9
x11 = x8 ^ 2
x12 = G * x0
x13 = tanh(x12)
x14 = 3 * u
x15 = G * x14
x16 = -Sp * x4 + 4 * x5 + 2 * x6 + 1
x17 = x9 * (B * x14 - Bp)
x18 = sech(x12)
x19 = -Fyp
x20 = Fy * u
x21 = x19 + x20
x22 = Sh * u
x23 = 2 * x4
x24 = S * x23
x25 = x24 * xi_y
x26 = 2 * x5
x27 = x26 - 1
x28 = x0 * (Fy * x27 + x22 + x25)
x29 = 2 * x0
x30 = x11 * x2
x31 = x30 * (Sd * x29 + x7 ^ 2)
x32 = u * xi_y
x33 = 5 * x4
x34 = x8 * (-2 * Fyph + u * (Fy ^ 2 * x33 + 2 * Fyh + Fyp ^ 2 - 4 * Fyp * x20) - 2 * x32 * (-x19 - 2 * x20))
x35 = 2 * u
x36 = x18 * x35
x37 = -Fxp
x38 = Fx * u
x39 = x37 + x38
x40 = St * u
x41 = x24 * xi_x
x42 = 2 * x38
x43 = u * xi_x
x44 = (4 * x39 * x4 * (Fx * x27 + x40 + x41) - x8 * (-2 * Fxpt + u * (Fx ^ 2 * x33 + Fxp ^ 2 - 4 * Fxp * x38 + 2 * Fxt) - 2 * x43 * (-x37 - x42))) * exp(2 * x1)
x45 = u ^ 4
x46 = Fxh * u
x47 = 3 * G
x48 = Fyt * u
x49 = Fx * Fyp
x50 = x4 * x47
x51 = Fxp * Fy
x52 = 4 * x2
x53 = sinh(x12)
x54 = u * x53
x55 = Fx * x45
x56 = cosh(x12)
x57 = Bp * x0
x58 = 3 * B
x59 = x45 * x58
x60 = u ^ 5 * x58
x61 = 3 * x1
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x10 * x3
BB[1,2] = 0
CC[1,1] = -u * x11 * x3 * (x13 * x4 * x8 * (-Gp + x15) + x16)
CC[1,2] = -x0 * x13 * x17 * x3
SS[1] = -8 * x18 * x21 * x28 + x18 * x45 * x52 * x8 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp - Gh * x38 + Gt * x20 + x32 * (-Fx * Gp + 3 * Fxp * G) + x46 * x47 - x47 * x48 - x49 * x50 + x50 * x51 + xi_x * (-Fyp * x15 + Gp * x20)) - x31 * (12 * B * u - 4 * Bp) + x34 * x36 + x36 * x44
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x10 * x52
CC[2,1] = x17 * x2 * x29 * sinh(2 * x12)
CC[2,2] = -4 * u * x16 * x30
SS[2] = -x0 * x52 * x56 * (Fx * Sh * x4 - Fxp * x22 + 4 * Fy * S * x55 + Fy * St * x4 - Fy * x42 - Fyp * x40 + x21 * x41 + x25 * x39 - x26 * x49 - x26 * x51 + x49 + x51) + x2 * x35 * x56 * x8 * (-Bh * Fxp * x0 + Bh * x55 - Bt * Fy * x45 + Bt * Fyp * x0 + 5 * Fx * Fy * x0 + Fxh * x57 - Fxh * x59 + Fxp * Fyp * u - Fxph - Fypt - Fyt * x57 + Fyt * x59 - x23 * x49 - x23 * x51 + x32 * (Bp * Fx * x0 - Fxp * x61 - Fxp + x42) - x43 * (-Fyp * x61 + Fyp + x20 * (Bp * x4 - 2)) + x46 + x48 + x49 * x60 - x51 * x60) + 4 * x21 * x28 * x53 + x31 * (-6 * G * u + 2 * Gp) - x34 * x54 + x44 * x54

    
    
    nothing
end



function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,S0_xx, S0_yy, S0_xy, S0_txx, S0_tyy, S0_txy,u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")
    @tilde_inner("Sp")
    @tilde_inner("Gp")
    @tilde_inner("Bp")

    @hat_inner("Bp")
    @hat_inner("Gp")
    @hat_inner("Sp")
    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")


x0 = S0_t * u
x1 = S0 * xi
x2 = u ^ 3
x3 = S * x2
x4 = S0 + u * x1 + x0 + x3
x5 = x4 ^ 3
x6 = u ^ 2
x7 = B * x2
x8 = exp(x7)
x9 = 8 * x8
x10 = x6 * x9
x11 = x4 ^ 2
x12 = Sp * u
x13 = S * x6
x14 = S0_t + x1
x15 = 1 / u
x16 = S0 * (x15 + xi) + S0_t + x13
x17 = x16 ^ 3 * x9
x18 = u ^ 4
x19 = 1 / x18
x20 = 4 * S0
x21 = -Sp
x22 = S * u
x23 = -x21 - 2 * x22
x24 = S0 + x23 * x6
x25 = x11 * x24
x26 = S0_x ^ 2
x27 = S0_y ^ 2
x28 = S0 * S0_xx
x29 = S0 * S0_yy
x30 = S0 ^ 2
x31 = S0_t ^ 2
x32 = S0 ^ 4
x33 = S0 ^ 3
x34 = 2 * S0_t
x35 = 1 / x33
x36 = x24 * x35
x37 = x16 ^ 2
x38 = x37 * x8
x39 = S0 * S0_ty
x40 = 2 * x39
x41 = S * x18
x42 = u * xi
x43 = 2 * x3
x44 = x0 * (S0 - x43) + x30 * (x42 + 1)
x45 = Sh * u
x46 = Fy * S0
x47 = 2 * xi_y
x48 = x30 * x6
x49 = S0_y * x44 + x40 * x41 + x48 * (Fy * x43 + x13 * x47 + x45 - x46)
x50 = 1 / x32
x51 = 1 / x6
x52 = G * x2
x53 = cosh(x52)
x54 = 4 * x53
x55 = x50 * x51 * x54
x56 = S0 * S0_tx
x57 = 2 * x56
x58 = Fx * S0
x59 = 2 * xi_x
x60 = S0_x * x44 + x41 * x57 + x48 * (Fx * x43 + St * u + x13 * x59 - x58)
x61 = sinh(x52)
x62 = x50 * x61
x63 = -Fyp
x64 = Fy * u
x65 = x63 + x64
x66 = -x65
x67 = -Fxp
x68 = Fx * u
x69 = x67 + x68
x70 = 1 / x30
x71 = S0_t * S0_y
x72 = x39 - x71
x73 = x70 * x72
x74 = x64 + xi_y
x75 = x73 + x74
x76 = S0_t * S0_x
x77 = x56 - x76
x78 = x70 * x77
x79 = x68 + xi_x
x80 = x78 + x79
x81 = 2 * u
x82 = 3 * u
x83 = G * x82
x84 = 3 * G
x85 = S0_y * x0
x86 = Fy * x6
x87 = x30 * (Gh + x83 * xi_y + x84 * x86) + x39 * x83 - x84 * x85
x88 = S0_x * x0
x89 = Fx * x6
x90 = x30 * (Gt + x83 * xi_x + x84 * x89) + x56 * x83 - x84 * x88
x91 = u ^ 5
x92 = 2 * x91
x93 = x50 * x92
x94 = 2 * x64
x95 = S0_ty * x6
x96 = 2 * x12
x97 = -Spp + x96
x98 = -x97
x99 = 1 / S0
x100 = x98 * x99
x101 = S0_t * x6
x102 = x70 * (x101 * x97 + x30)
x103 = x27 * x34
x104 = u * x33
x105 = x46 * x6
x106 = Fy ^ 2
x107 = x2 * x33
x108 = 2 * S0_y
x109 = Fyh * x104 - S0_t * x29 + S0_tyy * x30 + x103 - x105 * x71 + x106 * x107 + x33 * x86 * xi_y + x33 * xi_yy + x39 * (x105 - x108)
x110 = x21 + 3 * x22
x111 = -x110
x112 = x111 * x2
x113 = -x112 * x34 + x30
x114 = u * x111
x115 = B * x82
x116 = 3 * B
x117 = B * x6
x118 = 3 * x117
x119 = x115 * x39 - x116 * x85 + x30 * (Bh + Fy * x118 + x115 * xi_y)
x120 = Gp * x32
x121 = u * x120
x122 = u * x32
x123 = Fxp * x122
x124 = Fyp * x122
x125 = S0_x * x71
x126 = 6 * G
x127 = S0_xy * x0
x128 = 6 * u
x129 = G * x128
x130 = S0_ty * S0_x
x131 = 6 * Gh
x132 = x30 * x88
x133 = 6 * Gt
x134 = x30 * x85
x135 = Fx * Gp
x136 = x33 * x95
x137 = x32 * x6
x138 = x135 * x137
x139 = Fy * Gp
x140 = x137 * x139
x141 = 24 * G
x142 = x141 * x39 * x6
x143 = B * Fy
x144 = x143 * x91
x145 = Gt * x32
x146 = 3 * x145
x147 = S0_ty * x33
x148 = B * x18
x149 = x148 * xi_y
x150 = Fx * x91
x151 = x150 * x32
x152 = Bh * x84
x153 = x18 * x32
x154 = x153 * xi_x
x155 = 2 * x32
x156 = x155 * x2
x157 = x32 * x52
x158 = 3 * x157
x159 = Fy * x32
x160 = x159 * x52
x161 = 3 * x160
x162 = Fxp * x84
x163 = x137 * x162
x164 = Fyp * x84
x165 = x137 * x164
x166 = Gp - x83
x167 = x166 * x33
x168 = 2 * S0_txy
x169 = 3 * Gh
x170 = x30 * x76
x171 = x148 * x170
x172 = x30 * x84
x173 = 30 * x52
x174 = Fx * x71
x175 = x174 * x30
x176 = Fy * x76
x177 = x176 * x30
x178 = x48 * x76
x179 = x141 * xi_y
x180 = x141 * xi_x
x181 = x48 * x71
x182 = -2 * S0 * x166
x183 = u * xi_y
x184 = 30 * G
x185 = x6 * x84
x186 = x50 * x8
x187 = S0_xy * x155
x188 = S0 ^ 5
x189 = x188 * x6
x190 = Fx * x2
x191 = Sp * x153
x192 = x2 * xi_x
x193 = Fy * x188
x194 = Fx * x18
x195 = 4 * Sh
x196 = S * x32
x197 = x196 * x92
x198 = 2 * S0_x
x199 = 4 * St
x200 = x159 * x91
x201 = x18 * x199
x202 = 4 * xi_x
x203 = x153 * xi_y
x204 = Sp * x150
x205 = 12 * x91
x206 = S * x205
x207 = S0 * Sp
x208 = u ^ 7
x209 = Fx * S
x210 = u ^ 6
x211 = 14 * x210
x212 = x196 * xi_y
x213 = S * xi_x
x214 = 4 * S0_t
x215 = x30 * x41
x216 = x33 * xi_x
x217 = x170 * x18
x218 = Sp * x30
x219 = x2 * x218
x220 = x30 * x71
x221 = S * x211
x222 = 12 * x13
x223 = 2 * x30
x224 = x63 + x94
x225 = -x224
x226 = -Bp + x115
x227 = -S0 * x226
x228 = 2 * x6
x229 = S0_ty ^ 2 * x48
x230 = 12 * x117
x231 = x27 * x31
x232 = 4 * B
x233 = Bh + x232 * x86
x234 = x183 * x232 + x233
x235 = xi_y ^ 2
x236 = x128 * xi_y
x237 = x32 * x83
x238 = Gh * x32
x239 = Gp * S0
x240 = Gp * x30
x241 = S0 * x0 * x126
x242 = x0 * x172
x243 = x185 * x32
x244 = x131 * x32
x245 = Bh * x126
x246 = 18 * G
x247 = x153 * x246
x248 = 12 * G
x249 = x137 * x248
x250 = x248 * x6
x251 = x210 * x32
x252 = 36 * G
x253 = x143 * x252
x254 = x246 * x32
x255 = x208 * x254
x256 = B * x91
x257 = x254 * x256
x258 = x246 * x256
x259 = 3 * x7
x260 = x259 - 2
x261 = 6 * x7
x262 = x50 * x6 * x60
x263 = x26 * x34
x264 = x58 * x6
x265 = Fx ^ 2
x266 = x115 * x56 - x116 * x88 + x30 * (Bt + Fx * x118 + x115 * xi_x)
x267 = S0_t * S0_xx
x268 = x133 * x32
x269 = x148 * xi_x
x270 = Bt * x126
x271 = Fx * Fxp
x272 = Fx * xi_x
x273 = x2 * x265
x274 = xi_x ^ 2
x275 = x26 * x31
x276 = Fx * G
x277 = 36 * x276
x278 = B * x277
x279 = B * x265
x280 = G * xi_x
x281 = S0_tx ^ 2
x282 = x259 + 2
x283 = x126 * x282
x284 = 2 * x2
x285 = Bp * x6
x286 = 12 * Bt
x287 = 54 * x256
x288 = Bp * x194
x289 = 12 * Gt
x290 = x210 * x289
x291 = 24 * x148
x292 = G ^ 2
x293 = 36 * x292
x294 = u ^ 8 * x293
x295 = 18 * x292
x296 = x208 * x295
x297 = x70 * x76
x298 = Fx * x297
x299 = x26 * x35
x300 = x275 * x50
x301 = x210 * x295
ABCS[1] = 0
ABCS[2] = -x10 * x5
ABCS[3] = -x10 * x11 * (-x12 + 3 * x13 + x14)

 
if u == 0.0 
	println("entered A u zero")
	ABCS[4] = -4*S0^3*Sp
else
	ABCS[4] = x6 * (S0 * x19 * x5 * x9 + S0_t * x17 + u * x37 * x50 * (-x228 * x61 * (B * x106 * x255 + Bh * Gh * x156 + Fy * x173 * x220 + Fyh * x121 - Fyh * x243 + Fyp * x161 + Gh * x124 - Gs * x32 - S0_t * S0_yy * x240 + S0_ty * (S0 * S0_y * (-x101 * x248 * x260 + x182) + x104 * (6 * Bh * x52 + Fyp * x83 + x131 * (x7 - 1) + x183 * x248 * x260 + x64 * (Gp + x129 * (x261 - 5)))) + S0_tyy * x167 + S0_yy * x242 + x103 * x239 + x106 * x120 * x2 - x106 * x247 + x120 * xi_yy + x126 * x229 * x260 + x131 * x134 - x131 * x148 * x220 - x139 * x48 * x71 + x140 * xi_y + x144 * x244 + x149 * x244 - 30 * x160 * xi_y - x164 * x181 + x165 * xi_y + x179 * x181 - x18 * x220 * x245 + x200 * x245 + x203 * x245 - x210 * x220 * x253 - x220 * x252 * x256 * xi_y - x231 * x250 + x231 * x258 - x235 * x249 + x235 * x257 - x236 * x238 - x237 * xi_yy - 7 * x238 * x86 - x241 * x27 + x251 * x253 * xi_y) + x53 * (x109 * x227 * x228 + x119 ^ 2 * x92 + x119 * x2 * x223 * x66 + x122 * x66 ^ 2 + x223 * (u * x224 * x39 + x225 * x85 - x30 * (Fyph - u * (Fyh + x224 * x64) + x183 * x225)) - x228 * (12 * B * x229 + x128 * x39 * (-x232 * x85 + x234 * x30) + x230 * x231 - 6 * x234 * x30 * x85 + x32 * (Bs + x230 * x235 + x233 * x236 + 6 * x86 * (Bh + x143 * x228))) + x87 ^ 2 * x92)) + x1 * x17 - x11 * x186 * x53 * x81 * (3 * B * Fx * Gh * x32 * x91 + 3 * B * Gh * x18 * x32 * xi_x + 3 * B * Gt * S0_t * S0_y * x18 * x30 + 3 * Bh * G * S0_t * S0_x * x18 * x30 - Bh * Gt * x2 * x32 + 3 * Bt * Fy * G * x32 * x91 + 3 * Bt * G * S0_ty * x18 * x33 + 3 * Bt * G * x18 * x32 * xi_y + Bt * Gh * x2 * x32 - Bt * x172 * x18 * x71 + 36 * Fx * Fy * G * x18 * x32 - Fx * Fyp * x158 + 30 * Fx * G * S0_ty * x2 * x33 + 30 * Fx * G * x2 * x32 * xi_y + 7 * Fx * Gh * x32 * x6 + Fx * Gp * S0_t * S0_y * x30 * x6 + 3 * Fxh * G * x32 * x6 - Fxh * x121 + 3 * Fxp * G * S0_t * S0_y * x30 * x6 - Fxp * x161 + 30 * Fy * G * x2 * x32 * xi_x + Fy * Gp * S0_t * S0_x * x30 * x6 + 7 * Fy * Gt * x32 * x6 - Fy * x135 * x156 + 3 * Fyp * G * S0_t * S0_x * x30 * x6 + 3 * Fyt * G * x32 * x6 - Fyt * x121 + 12 * G * S0 * S0_t * S0_x * S0_y * u + 24 * G * S0_ty * x33 * x6 * xi_x + 24 * G * S0_x * S0_y * x31 * x6 + 6 * G * u * x32 * xi_xy + 24 * G * x32 * x6 * xi_x * xi_y + 2 * Gc * x32 + 6 * Gh * u * x32 * xi_x - Gh * x123 + 2 * Gp * S0_t * S0_xy * x30 + 2 * Gp * S0_ty * S0_x * x30 - Gp * x125 * x20 + 6 * Gt * S0_ty * u * x33 + 6 * Gt * u * x32 * xi_y - Gt * x124 - 3 * Gt * x147 * x148 + S0 * S0_tx * (S0_y * (-x101 * x141 - x182) + u * x30 * (-u * (Bh * x185 + x139 + x164 - x184 * x64) + x141 * x183 + x169 * (x7 + 2)) + x142) - 2 * x120 * xi_xy - x126 * x127 * x30 - x129 * x130 * x30 - x131 * x132 - x133 * x134 - x135 * x136 - x136 * x162 - x138 * xi_y - x140 * xi_x - x142 * x76 - x144 * x146 - x146 * x149 - x151 * x152 - x152 * x154 - x163 * xi_y - x165 * xi_x - x167 * x168 - x169 * x171 - x173 * x175 - x173 * x177 - x178 * x179 - x180 * x181) + x14 * x25 * x9 / x2 + 4 * x15 * x16 * x186 * (-x2 * x53 * (x49 * x90 + x60 * x87) + x61 * (Fx * Fy * Sp * x155 * x210 - Fx * x211 * x212 + Fxh * x189 + Fxh * x191 - Fxh * x197 + Fyt * x189 + Fyt * x191 - Fyt * x197 - 4 * S * x153 * xi_xy - 8 * S0 * x125 * x41 + S0_ty * x190 * x32 - S0_ty * x198 * x219 - S0_ty * x206 * x216 - S0_x * S0_y * x206 * x31 + S0_xy * x214 * x215 - S0_xy * x219 * x34 - Sc * x156 - Sh * x153 * x202 + Sp * x151 * xi_y + Sp * x156 * xi_xy + Sp * x200 * xi_x + x107 * x168 * x23 - x107 * x174 - x107 * x176 + x108 * x32 * x89 + 4 * x125 * x2 * x207 - 2 * x127 * x33 + 4 * x130 * x215 - x147 * x201 + x147 * x204 - x147 * x209 * x211 - x151 * x195 - 16 * x159 * x208 * x209 - x159 * x211 * x213 + x170 * x206 * xi_y + x175 * x221 - x176 * x218 * x91 + x177 * x221 - x187 * x42 - x187 + x188 * x190 * xi_y + x192 * x193 + 2 * x193 * x194 + x195 * x217 + x198 * x32 * x86 - x199 * x200 - x199 * x203 + x2 * x56 * (S0_y * (S0_t * x222 + x20 * x22 - 2 * x207) - x222 * x39 + x30 * (Fy * (S0 + x6 * (-x21 - 14 * x22)) - x222 * xi_y - 4 * x45)) + x201 * x220 - x204 * x220 - x205 * x212 * xi_x + x205 * x213 * x220 + x206 * x39 * x76)) - 12 * x16 ^ 4 * x8 + x16 * (4 * x49 * x50 * x6 * x61 * x87 - x54 * (-S0 * xi_yy - S0_tyy - S0_y * x47 - S0_yy * x15 - S0_yy * xi + x109 * x24 * x35 + x119 * x49 * x50 * x6 + x50 * (x30 * x74 + x72) * (S0_y * x113 + x112 * x40 + x48 * (Sph + x114 * x47 - x81 * (Sh + x110 * x64))) - x6 * (-Spp * x75 ^ 2 + Ss + (Sph + Spp * x75) * (x47 + 2 * x73 + x94)) + x75 * (S0_y * x102 + x100 * x95 + x6 * (Sph + u * (Fy * Spp - Fy * x96 - 2 * Sh) + x98 * xi_y)))) + x19 * x20 * x25 * x8 + 4 * x36 * x38 * (x26 + x27 - x28 - x29 + x30 * x31 + x32 * xi ^ 2 + x33 * x34 * xi) + x38 * x61 * x81 * (Fxph + Fxpp * x75 + Fypp * x80 + Fypt + u * x66 * x69 - u * (Fxh + Fxp * x75) - u * (Fyp * x80 + Fyt) - x75 * (Fxpp + x69 * x81) - x80 * (Fypp + x65 * x81) - x87 * x90 * x93) - x49 ^ 2 * x55 + x49 * x51 * x60 * x62 * x9 + (x16 * (4 * x262 * x61 * x90 + x54 * (S0 * xi_xx + S0_txx + S0_x * x59 + S0_xx * x15 + S0_xx * xi + x262 * x266 - x36 * (Fxt * x104 - S0_t * x28 + S0_txx * x30 + x107 * x265 + x216 * x89 + x263 - x264 * x76 + x33 * xi_xx + x56 * (-x198 + x264)) - x50 * (x30 * x79 + x77) * (S0_x * x113 + x112 * x57 + x48 * (Spt + x114 * x59 - x81 * (St + x110 * x68))) + x6 * (Sb - Spp * x80 ^ 2 + (Spp * x80 + Spt) * (x59 + 2 * x68 + 2 * x78)) - x80 * (S0_tx * x100 * x6 + S0_x * x102 + x6 * (Spt + u * (Fx * Spp - Fx * x96 - 2 * St) + x98 * xi_x)))) + x37 * (u * x53 * (Bb * x228 - Bp * Fxt * x284 - Bp * x265 * x92 - Fxp * u * x59 + Fxp * x0 * x198 * x70 - 2 * Fxpt + 6 * Fxt * x148 + Fxt * x81 - G * x290 * x297 + Gt ^ 2 * x92 + 12 * S0_t * x299 * x7 + S0_tx * x35 * x81 * (u * x198 * (-S0_t * (x230 + 9 * x292 * x91) + x227) + x30 * (u * (Fx * (-x285 + x301 + 27 * x7 + 2) + x128 * (Bt + Gt * x52)) + x67 + xi_x * (x301 + 24 * x7))) + S0_txx * x226 * x228 * x99 + S0_xx * x285 * x34 * x70 + u ^ 9 * x265 * x295 + u * x69 ^ 2 + 6 * x18 * x281 * x70 * (3 * x2 * x292 + x232) + x192 * x286 + x194 * x286 - x2 * x286 * x297 + x202 * x89 + x208 * x276 * x289 - x208 * x293 * x297 * xi_x + 30 * x210 * x279 - x214 * x285 * x299 - x228 * x271 - x261 * x267 * x70 + x261 * xi_xx + x266 ^ 2 * x93 + x266 * x284 * x69 * x70 - 48 * x269 * x297 + x272 * x287 + x272 * x294 + 4 * x273 + x274 * x291 + x274 * x296 + x280 * x290 - 2 * x285 * xi_xx - x287 * x298 + 2 * x288 * x297 - x288 * x59 + x291 * x300 - x294 * x298 + x296 * x300 - 4 * x297 * x89) + x284 * x62 * (Bt * Gt * x156 - Fx * x170 * x173 + Fx * x256 * x268 - Fxt * x121 + Fxt * x243 + Gb * x32 - Gt * x123 + S0_tx * (S0 * x198 * (S0 * x166 - x101 * x283) + x104 * (12 * u * x280 * x282 + u * (-x135 + x148 * x277 - x162 + x184 * x68 + x270 * x6) + x133 * (x7 + 1))) - S0_txx * x167 - S0_xx * x242 - x120 * x273 - x120 * xi_xx + x128 * x145 * xi_x - x132 * x133 - x133 * x171 + x135 * x48 * x76 - x138 * xi_x + 7 * x145 * x89 + x151 * x270 + x154 * x270 + 30 * x157 * x272 - x158 * x271 + x162 * x178 - x163 * xi_x - x170 * x210 * x278 - 36 * x170 * x256 * x280 - x178 * x180 - x217 * x270 + x237 * xi_xx - x239 * x263 + x240 * x267 + x241 * x26 + x247 * x265 + x249 * x274 + x250 * x275 + x251 * x278 * xi_x + x255 * x279 + x257 * x274 + x258 * x275 + x268 * x269 + x281 * x283 * x48)) - x55 * x60 ^ 2) * exp(2 * x7))
	
end

    nothing
    
end
