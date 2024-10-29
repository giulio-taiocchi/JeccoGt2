
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
    #return esc( :($fs = $f_yy + (Fy * u + xi_y) * ( -2*($fp_y) + (Fy * u + xi_y)* ($fpp) )) )
    return esc( :($fs = $f_yy + (Fy * u + xi_y) * ( -2*($fp_y) - (Fy * u + xi_y)* ($fpp) )) )
end

macro cross_inner(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    #return esc( :($fc = $f_xy  - (Fx * u + xi_x) * ($fp_y) -
    #              (Fy * u + xi_y) * ( $fp_x -(Fx * u + xi_x) * ($fpp) ) ) )
    return esc( :($fc = $f_xy  - (Fx * u + xi_x) * ($fp_y) -
                  (Fy * u + xi_y) * ( $fp_x -(Fx * u + xi_x) * ($fpp) ) ) )
end



# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    ( u, xi, B, Bp, G, Gp) = vars

 
  
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
        u, xi, xi_x, xi_y,
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
        u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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

    
    @tilde_inner("Sp")
    @tilde_inner("Fxp")
    @tilde_inner("Fyp")
    @tilde_inner("Bp")
    @tilde_inner("Sp")
    @tilde_inner("Gp")

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

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")

x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = S * x0
x4 = u * xi
x5 = x4 + 1
x6 = x3 + x5
x7 = x6 ^ 4
x8 = x2 * x7
x9 = 2 * x0
x10 = u ^ 2
x11 = 4 * x2
x12 = Sp * x11
x13 = Sd * x2
x14 = 8 * x13
x15 = 16 * Sp * x2
x16 = G * x0
x17 = cosh(x16)
x18 = sinh(x16)
x19 = Gs * x18
x20 = 4 * x17
x21 = Ss * x20
x22 = u ^ 5
x23 = u ^ 4
x24 = u ^ 6
x25 = x17 * x23
x26 = 2 * Bp * Fyh
x27 = x17 * x9
x28 = Bp * xi_yy
x29 = x0 * x17
x30 = 2 * Fy
x31 = Fyp * x30
x32 = 8 * x22
x33 = Gh * x18
x34 = Fy * x33
x35 = 24 * x22
x36 = Fy * x17
x37 = 6 * G
x38 = x18 * x37
x39 = Fyh * x22
x40 = S * x20
x41 = x2 * x23
x42 = Gd * x37
x43 = x38 * xi_yy
x44 = 4 * x24
x45 = 12 * xi_y
x46 = x33 * x45
x47 = S * x24
x48 = 4 * x19
x49 = x23 * xi
x50 = 8 * x2
x51 = Sp * x50
x52 = 8 * xi_yy
x53 = x0 * x14
x54 = 16 * x13
x55 = 16 * Sh
x56 = Bh ^ 2
x57 = x17 * x24
x58 = 2 * x57
x59 = Fy ^ 2
x60 = Fyp ^ 2
x61 = Gh ^ 2
x62 = S ^ 2
x63 = S ^ 3
x64 = u ^ 7
x65 = S ^ 4
x66 = u ^ 10
x67 = u ^ 8
x68 = x17 * x67
x69 = Fy * x68
x70 = 12 * Bh
x71 = B * x70
x72 = B * x17
x73 = x45 * x64
x74 = u ^ 9
x75 = x17 * xi_y
x76 = Bp * x75
x77 = Bp * Fyh
x78 = x20 * x24
x79 = S * x78
x80 = x20 * x49
x81 = Fy * Fyp
x82 = G * Gh
x83 = G * x18
x84 = 12 * x83
x85 = x67 * x84
x86 = Fy * Sh
x87 = x24 * x83
x88 = Fy * xi_y
x89 = 28 * S
x90 = 20 * x34
x91 = x24 * xi
x92 = S * x57
x93 = x57 * xi
x94 = S * x85
x95 = Fyh * xi
x96 = 8 * x92
x97 = 24 * x2
x98 = x64 * x97
x99 = G * Gd
x100 = x98 * x99
x101 = x22 * xi
x102 = x97 * x99
x103 = G * x17
x104 = Gh * x103
x105 = S * x64
x106 = x84 * xi_yy
x107 = 4 * x33
x108 = S * x74
x109 = x33 * xi_y
x110 = 32 * x109
x111 = x64 * xi
x112 = S * x111
x113 = Sp * x54
x114 = Sp * x97
x115 = S * x114
x116 = S * x17
x117 = xi_y ^ 2
x118 = x108 * x20
x119 = x111 * x20
x120 = Bp * x59
x121 = xi ^ 2
x122 = x59 * x64
x123 = 42 * x116
x124 = 8 * x62 * x68
x125 = x117 * x35
x126 = x62 * x74
x127 = 2 * x19
x128 = x121 * x22
x129 = xi ^ 3
x130 = xi ^ 4
x131 = x24 * x62
x132 = x12 * x23
x133 = 24 * B
x134 = Bh * x133
x135 = u ^ 11
x136 = x116 * x135
x137 = Fy * x136
x138 = x74 * xi
x139 = x138 * x36
x140 = Bh * xi_y
x141 = x133 * x140 * xi
x142 = Bp * x88
x143 = x40 * x67
x144 = x111 * x40
x145 = S * x135
x146 = 24 * x104
x147 = Fy * x146
x148 = x84 * x86
x149 = x83 * x88
x150 = 32 * x108
x151 = x116 * x88
x152 = x108 * x84
x153 = 72 * x99
x154 = S * x2
x155 = x67 * xi
x156 = x154 * x155
x157 = S * x66
x158 = x146 * xi_y
x159 = x68 * xi
x160 = Sh * x45
x161 = S * x155
x162 = B ^ 2
x163 = 36 * x162
x164 = x36 * xi_y
x165 = x164 * x74
x166 = u ^ 12
x167 = x66 * xi
x168 = x167 * x40
x169 = x62 * x66
x170 = 2 * x17
x171 = x170 * x28
x172 = x17 * x31
x173 = G ^ 2
x174 = 36 * x173
x175 = x135 * x62
x176 = x121 * x64
x177 = 12 * x176
x178 = 72 * x83
x179 = x59 * x83
x180 = 48 * x155
x181 = S * x68
x182 = 30 * x181
x183 = Fyh * x38
x184 = 36 * x62
x185 = x2 * x66
x186 = x184 * x185
x187 = u ^ 13
x188 = x187 * x63
x189 = x2 * x42
x190 = u ^ 16
x191 = x190 * x65
x192 = x121 * x2
x193 = x192 * x99
x194 = 36 * x24
x195 = x130 * x67
x196 = S * x117
x197 = x117 * xi
x198 = x121 * x24
x199 = u ^ 15
x200 = x47 * x51
x201 = Sp * x14
x202 = x111 * x51
x203 = x116 * x166
x204 = x203 * xi
x205 = S * x166 * xi
x206 = x135 * xi
x207 = 18 * x162
x208 = x17 * x59
x209 = x208 * x66
x210 = x117 * x68
x211 = 2 * x56
x212 = x166 * x62
x213 = x17 * x212
x214 = x121 * x68
x215 = 18 * x173
x216 = 15 * x62
x217 = 2 * x61
x218 = u ^ 14
x219 = x218 * x62
x220 = x36 * x71
x221 = x121 * x66
x222 = x187 * x62
x223 = Bh * x45 * x72
x224 = x121 * x74
x225 = 72 * x162
x226 = x151 * x166
x227 = x164 * x167
x228 = x30 * x76
x229 = 12 * x219
x230 = Fy * x104
x231 = 12 * x221
x232 = x121 * x67
x233 = 72 * x173
x234 = 84 * S
x235 = x129 * x157
x236 = x2 * x62
x237 = x206 * x236
x238 = x218 * x63 * xi
x239 = x104 * x45
x240 = x11 * x18
x241 = x10 * x240
x242 = x116 * x163
x243 = x187 * x59
x244 = x163 * xi
x245 = x135 * x208
x246 = x117 * x135
x247 = x117 * x17
x248 = x247 * x74
x249 = 2 * x214
x250 = 54 * x222
x251 = 30 * x224
x252 = x116 * x174
x253 = x174 * xi
x254 = x62 * x83
x255 = x151 * x187 * xi
x256 = exp(2 * x1)
x257 = Bb * x256
x258 = Gb * x256
x259 = Sb * x256
x260 = 2 * Gt
x261 = Bh * x2
x262 = x260 * x261
x263 = x18 * x41
x264 = Fx * Fy
x265 = x0 * x240
x266 = Fx * xi_y
x267 = Fxh * Gp
x268 = 2 * x2 * x25
x269 = x132 * x18
x270 = x265 * xi
x271 = 2 * x10
x272 = x18 * x2
x273 = Fxp * Fyp
x274 = Fyt * Gp
x275 = Gp * xi_xy
x276 = 8 * Sh
x277 = x164 * x184 * x199
x278 = x121 * x135 * x164
x279 = x116 * x218 * x59
x280 = x117 * x203
x281 = x166 * x184
x282 = 6 * B
x283 = x17 * x282
x284 = Fxt * x256
x285 = x22 * x284
x286 = x207 * x208
x287 = x190 * x62
x288 = x121 * x166
x289 = x219 * x247
x290 = x207 * x66
x291 = x121 * x247
x292 = Fx * x256
x293 = x292 * x32
x294 = Bt * x17
x295 = Bt * x256
x296 = Gt * x18
x297 = x296 * x44
x298 = Fxp * x292
x299 = St * x292
x300 = x208 * x215
x301 = x215 * x66
x302 = x18 * x258
x303 = St * x256
x304 = Gt * x2
x305 = x283 * x304
x306 = Fx * x37
x307 = x11 * x17
x308 = Bh * Gt
x309 = Gp * x264
x310 = x18 * x264
x311 = x310 * x51
x312 = 2 * x22
x313 = Fx * x2
x314 = Gp * x75
x315 = x12 * x18 * x266
x316 = x11 * x267
x317 = x116 * x64
x318 = x101 * x307
x319 = Fxh * x12 * x18
x320 = x11 * x274
x321 = Fyt * x12 * x18
x322 = x275 * x50
x323 = 8 * xi
x324 = Bt ^ 2
x325 = x256 * x58
x326 = Fx ^ 2
x327 = Fxp ^ 2
x328 = Gt ^ 2
x329 = x17 ^ 2
x330 = Bd * x329
x331 = x282 * x330
x332 = Bt * x68
x333 = 12 * B
x334 = x292 * x333
x335 = x296 * x67
x336 = x299 * x333
x337 = Bt * x292
x338 = x292 * x89
x339 = 20 * x292
x340 = x295 * x296
x341 = St * x295
x342 = 12 * Gt * x292
x343 = x284 * xi
x344 = 12 * x343
x345 = 4 * x296 * x303
x346 = x304 * x333
x347 = B * x304 * x45
x348 = x103 * x70
x349 = Fx * Gp * x36
x350 = Gp * x266
x351 = x11 * x350
x352 = x181 * xi
x353 = x256 * x326
x354 = x353 * x64
x355 = 18 * x354
x356 = x170 * x257
x357 = x118 * x256
x358 = x119 * x256
x359 = 2 * x302
x360 = B * x330
x361 = x360 * x98
x362 = x360 * x97
x363 = x170 * x236 * x66
x364 = x192 * x58
x365 = 2 * x273
x366 = x275 * x307
x367 = x133 * x337
x368 = 24 * x138
x369 = x337 * x72
x370 = x292 * x296
x371 = x133 * x370
x372 = x138 * x299
x373 = 24 * x83
x374 = x150 * x292 * xi
x375 = Gt * x103 * x292
x376 = 24 * x375
x377 = 36 * B
x378 = x353 * x83
x379 = 72 * B
x380 = x116 * x353
x381 = 48 * B
x382 = x353 * xi
x383 = x283 * x284
x384 = x175 * x339
x385 = x177 * x292
x386 = 4 * x340
x387 = x256 * x324
x388 = x170 * x298
x389 = x284 * x38
x390 = x256 * x328
x391 = x330 * x379
x392 = Gt * x282 * x36
x393 = x218 * x236
x394 = x192 * x66
x395 = Bh * x17 * x306
x396 = x170 * x350
x397 = x17 * x353
x398 = 2 * x213
x399 = x2 * x331
x400 = x192 * x360
x401 = x296 * x334
x402 = S * x379
x403 = x206 * x378
x404 = x103 * x342
x405 = x250 * x353
x406 = x187 * x353
x407 = x135 * x397
x408 = x218 * x380
x409 = x207 * x397
x410 = x215 * x397
x411 = 3 * G
x412 = u * x411
x413 = x18 * x6
x414 = u * x6
x415 = 2 * x4
x416 = x1 * x6
x417 = 3 * x416
x418 = 24 * x131 + 36 * x3 * x5 + 12 * x5 ^ 2
x419 = x6 ^ 2
x420 = 9 * x0
x421 = x173 * x24
x422 = 9 * x421
x423 = 3 * x421 + 2
x424 = 6 * S * x423 * x5 + x0 * x62 * (x422 + 4) + x420 * (G * x4 + G) ^ 2
x425 = 18 * x421
x426 = x16 * x419
x427 = 5 * x3 + 3 * x4 + x417 + 3
x428 = Gp * x6
x429 = 17 * x3 + 9 * x4 + 7
ABCS[1] = x8 * x9
ABCS[2] = 8 * x10 * x8
ABCS[3] = u * x11 * x7
ABCS[4] = 12 * B * Fx * Gh * S * x135 * x17 * x2 + 12 * B * Fx * Gh * S * x166 * x17 * x2 * xi + 6 * B * Fx * Gh * x121 * x17 * x2 * x66 + 6 * B * Fx * Gh * x17 * x2 * x218 * x62 + 6 * B * Fx * Gh * x17 * x2 * x67 + 12 * B * Fx * Gh * x17 * x2 * x74 * xi + 144 * B * Fy * G * S * x166 * x18 * xi_y + 144 * B * Fy * G * S * x18 * x187 * xi * xi_y + 72 * B * Fy * G * x121 * x135 * x18 * xi_y + 72 * B * Fy * G * x18 * x199 * x62 * xi_y + 144 * B * Fy * G * x18 * x66 * xi * xi_y + 72 * B * Fy * G * x18 * x74 * xi_y + 24 * B * Fy * Gh * S * x135 * x18 + 24 * B * Fy * Gh * S * x166 * x18 * xi + 12 * B * Fy * Gh * x121 * x18 * x66 + 12 * B * Fy * Gh * x18 * x218 * x62 + 12 * B * Fy * Gh * x18 * x67 + 24 * B * Fy * Gh * x18 * x74 * xi + 12 * B * Fy * S * Sh * x135 * x17 + 156 * B * Fy * S * x17 * x66 * xi * xi_y + 144 * B * Fy * S * x17 * x74 * xi_y + 12 * B * Fy * Sh * x17 * x67 + 12 * B * Fy * Sh * x17 * x74 * xi + 54 * B * Fy * x121 * x17 * x67 * xi_y + 102 * B * Fy * x166 * x17 * x62 * xi_y + 42 * B * Fy * x17 * x24 * xi_y + 96 * B * Fy * x17 * x64 * xi * xi_y + 12 * B * Fyh * S * x17 * x67 + 12 * B * Fyh * S * x17 * x74 * xi + 6 * B * Fyh * x121 * x17 * x64 + 6 * B * Fyh * x135 * x17 * x62 + 6 * B * Fyh * x17 * x22 + 12 * B * Fyh * x17 * x24 * xi + 72 * B * G * S * x117 * x135 * x18 + 72 * B * G * S * x117 * x166 * x18 * xi + 72 * B * G * S * x18 * x187 * x59 + 72 * B * G * S * x18 * x218 * x59 * xi + 36 * B * G * x117 * x121 * x18 * x66 + 36 * B * G * x117 * x18 * x218 * x62 + 36 * B * G * x117 * x18 * x67 + 72 * B * G * x117 * x18 * x74 * xi + 36 * B * G * x121 * x166 * x18 * x59 + 72 * B * G * x135 * x18 * x59 * xi + 36 * B * G * x18 * x190 * x59 * x62 + 36 * B * G * x18 * x59 * x66 + 24 * B * Gh * S * x135 * x18 * xi * xi_y + 24 * B * Gh * S * x18 * x66 * xi_y + 12 * B * Gh * x121 * x18 * x74 * xi_y + 12 * B * Gh * x18 * x187 * x62 * xi_y + 12 * B * Gh * x18 * x64 * xi_y + 24 * B * Gh * x18 * x67 * xi * xi_y - B * Gt * x116 * x185 * x45 + 12 * B * S * Sh * x17 * x66 * xi_y + 72 * B * S * x117 * x17 * x67 + 72 * B * S * x117 * x17 * x74 * xi + 84 * B * S * x135 * x17 * x59 * xi + 72 * B * S * x17 * x59 * x66 + 12 * B * S * x17 * x64 * xi_yy + 12 * B * S * x17 * x67 * xi * xi_yy + 12 * B * Sh * x17 * x64 * xi_y + 12 * B * Sh * x17 * x67 * xi * xi_y + 24 * B * x117 * x121 * x17 * x64 + 48 * B * x117 * x135 * x17 * x62 + 24 * B * x117 * x17 * x22 + 48 * B * x117 * x17 * x24 * xi + 6 * B * x121 * x17 * x24 * xi_yy + 30 * B * x121 * x17 * x59 * x74 + 54 * B * x17 * x187 * x59 * x62 + 12 * B * x17 * x22 * xi * xi_yy + 6 * B * x17 * x23 * xi_yy + 18 * B * x17 * x59 * x64 + 48 * B * x17 * x59 * x67 * xi + 6 * B * x17 * x62 * x66 * xi_yy - B * x184 * x190 * x378 - 84 * B * x206 * x380 + 24 * Bd * Bp * S * x121 * x2 * x329 * x67 + 8 * Bd * Bp * S * x129 * x2 * x329 * x74 + 8 * Bd * Bp * S * x2 * x24 * x329 + 24 * Bd * Bp * S * x2 * x329 * x64 * xi + 2 * Bd * Bp * x0 * x2 * x329 + 12 * Bd * Bp * x121 * x135 * x2 * x329 * x62 + 12 * Bd * Bp * x121 * x2 * x22 * x329 + 8 * Bd * Bp * x129 * x2 * x24 * x329 + 2 * Bd * Bp * x130 * x2 * x329 * x64 + 8 * Bd * Bp * x166 * x2 * x329 * x63 + 8 * Bd * Bp * x187 * x2 * x329 * x63 * xi + 2 * Bd * Bp * x199 * x2 * x329 * x65 + 8 * Bd * Bp * x2 * x23 * x329 * xi + 24 * Bd * Bp * x2 * x329 * x62 * x66 * xi + 12 * Bd * Bp * x2 * x329 * x62 * x74 + 24 * Bh * Fy * G * S * x135 * x18 + 24 * Bh * Fy * G * S * x166 * x18 * xi + 12 * Bh * Fy * G * x121 * x18 * x66 + 12 * Bh * Fy * G * x18 * x218 * x62 + 12 * Bh * Fy * G * x18 * x67 + 24 * Bh * Fy * G * x18 * x74 * xi + 28 * Bh * Fy * S * x17 * x67 + 32 * Bh * Fy * S * x17 * x74 * xi + 12 * Bh * Fy * x121 * x17 * x64 + 20 * Bh * Fy * x135 * x17 * x62 + 8 * Bh * Fy * x17 * x22 + 20 * Bh * Fy * x17 * x24 * xi + 24 * Bh * G * S * x135 * x18 * xi * xi_y + 24 * Bh * G * S * x18 * x66 * xi_y + 12 * Bh * G * x121 * x18 * x74 * xi_y + 12 * Bh * G * x18 * x187 * x62 * xi_y + 12 * Bh * G * x18 * x64 * xi_y + 24 * Bh * G * x18 * x67 * xi * xi_y + 8 * Bh * Gh * S * x18 * x66 * xi + 8 * Bh * Gh * S * x18 * x74 + 4 * Bh * Gh * x121 * x18 * x67 + 4 * Bh * Gh * x166 * x18 * x62 + 4 * Bh * Gh * x18 * x24 + 8 * Bh * Gh * x18 * x64 * xi + 4 * Bh * S * Sh * x17 * x74 + 32 * Bh * S * x17 * x64 * xi_y + 32 * Bh * S * x17 * x67 * xi * xi_y + 4 * Bh * Sh * x17 * x24 + 4 * Bh * Sh * x17 * x64 * xi + 12 * Bh * x121 * x17 * x24 * xi_y + 24 * Bh * x17 * x22 * xi * xi_y + 12 * Bh * x17 * x23 * xi_y + 20 * Bh * x17 * x62 * x66 * xi_y - Bh * x72 * x73 + 4 * Bp * Fxt * S * x17 * x256 * x64 + 4 * Bp * Fxt * S * x17 * x256 * x67 * xi + 2 * Bp * Fxt * x121 * x17 * x24 * x256 + 4 * Bp * Fxt * x17 * x22 * x256 * xi + 2 * Bp * Fxt * x17 * x23 * x256 + 2 * Bp * Fxt * x17 * x256 * x62 * x66 + 4 * Bp * S * x17 * x256 * x326 * x66 * xi + 4 * Bp * S * x17 * x256 * x326 * x74 + 2 * Bp * x121 * x17 * x256 * x326 * x67 - Bp * x143 * x95 + 2 * Bp * x166 * x17 * x256 * x326 * x62 + 2 * Bp * x17 * x24 * x256 * x326 + 4 * Bp * x17 * x256 * x326 * x64 * xi - Bp * x20 * x39 * xi - 2 * Bp * x208 * x212 + 4 * Bs * S * x17 * x24 + 4 * Bs * S * x17 * x64 * xi + 2 * Bs * x0 * x17 + 2 * Bs * x121 * x17 * x22 + 4 * Bs * x17 * x23 * xi + 2 * Bs * x17 * x62 * x74 + 12 * Bt * Fy * G * S * x135 * x17 * x2 + 12 * Bt * Fy * G * S * x166 * x17 * x2 * xi + 6 * Bt * Fy * G * x121 * x17 * x2 * x66 + 6 * Bt * Fy * G * x17 * x2 * x218 * x62 + 6 * Bt * Fy * G * x17 * x2 * x67 + 12 * Bt * Fy * G * x17 * x2 * x74 * xi + 12 * Bt * G * S * x135 * x17 * x2 * xi * xi_y + 12 * Bt * G * S * x17 * x2 * x66 * xi_y + 6 * Bt * G * x121 * x17 * x2 * x74 * xi_y + 6 * Bt * G * x17 * x187 * x2 * x62 * xi_y + 6 * Bt * G * x17 * x2 * x64 * xi_y + 12 * Bt * G * x17 * x2 * x67 * xi * xi_y + 4 * Bt * Gh * S * x17 * x2 * x66 * xi + 4 * Bt * Gh * S * x17 * x2 * x74 + 2 * Bt * Gh * x121 * x17 * x2 * x67 + 2 * Bt * Gh * x166 * x17 * x2 * x62 + 2 * Bt * Gh * x17 * x2 * x24 + 4 * Bt * Gh * x17 * x2 * x64 * xi - Bt * St * x256 * x78 - Bt * x339 * x93 + 168 * Fx * Fy * G * S * x135 * x17 * x2 * xi + 144 * Fx * Fy * G * S * x17 * x2 * x66 + 60 * Fx * Fy * G * x121 * x17 * x2 * x74 + 108 * Fx * Fy * G * x17 * x187 * x2 * x62 + 36 * Fx * Fy * G * x17 * x2 * x64 + 96 * Fx * Fy * G * x17 * x2 * x67 * xi + 72 * Fx * Fy * S * x173 * x18 * x187 * x2 + 72 * Fx * Fy * S * x173 * x18 * x2 * x218 * xi + 84 * Fx * Fy * S * x18 * x2 * x64 + 60 * Fx * Fy * S * x18 * x2 * x67 * xi + 36 * Fx * Fy * x121 * x166 * x173 * x18 * x2 + 72 * Fx * Fy * x135 * x173 * x18 * x2 * xi + 36 * Fx * Fy * x173 * x18 * x190 * x2 * x62 + 36 * Fx * Fy * x173 * x18 * x2 * x66 + 30 * Fx * Fy * x18 * x2 * x62 * x66 + 4 * Fx * Fyp * S * x18 * x2 * x24 + 4 * Fx * Fyp * S * x18 * x2 * x64 * xi + 2 * Fx * Fyp * x0 * x18 * x2 + 2 * Fx * Fyp * x121 * x18 * x2 * x22 + 4 * Fx * Fyp * x18 * x2 * x23 * xi + 2 * Fx * Fyp * x18 * x2 * x62 * x74 + 24 * Fx * G * Gh * S * x135 * x18 * x2 + 24 * Fx * G * Gh * S * x166 * x18 * x2 * xi + 12 * Fx * G * Gh * x121 * x18 * x2 * x66 + 12 * Fx * G * Gh * x18 * x2 * x218 * x62 + 12 * Fx * G * Gh * x18 * x2 * x67 + 24 * Fx * G * Gh * x18 * x2 * x74 * xi + 12 * Fx * G * S * Sh * x135 * x17 * x2 + 156 * Fx * G * S * x17 * x2 * x66 * xi * xi_y + 144 * Fx * G * S * x17 * x2 * x74 * xi_y + 12 * Fx * G * Sh * x17 * x2 * x67 + 12 * Fx * G * Sh * x17 * x2 * x74 * xi + 54 * Fx * G * x121 * x17 * x2 * x67 * xi_y + 102 * Fx * G * x166 * x17 * x2 * x62 * xi_y + 42 * Fx * G * x17 * x2 * x24 * xi_y + 96 * Fx * G * x17 * x2 * x64 * xi * xi_y + 28 * Fx * Gh * S * x17 * x2 * x67 + 32 * Fx * Gh * S * x17 * x2 * x74 * xi + 12 * Fx * Gh * x121 * x17 * x2 * x64 + 20 * Fx * Gh * x135 * x17 * x2 * x62 + 8 * Fx * Gh * x17 * x2 * x22 + 20 * Fx * Gh * x17 * x2 * x24 * xi + 72 * Fx * S * x166 * x173 * x18 * x2 * xi_y + 72 * Fx * S * x173 * x18 * x187 * x2 * xi * xi_y + 68 * Fx * S * x18 * x2 * x24 * xi_y + 56 * Fx * S * x18 * x2 * x64 * xi * xi_y + 24 * Fx * Sh * x18 * x2 * x22 + 16 * Fx * Sh * x18 * x2 * x24 * xi - Fx * x108 * x11 * x314 * xi + 36 * Fx * x121 * x135 * x173 * x18 * x2 * xi_y - Fx * x135 * x154 * x348 - Fx * x154 * x166 * x348 * xi + 36 * Fx * x173 * x18 * x199 * x2 * x62 * xi_y + 72 * Fx * x173 * x18 * x2 * x66 * xi * xi_y + 36 * Fx * x173 * x18 * x2 * x74 * xi_y - Fx * x18 * x192 * x24 * x30 + 24 * Fx * x18 * x2 * x62 * x74 * xi_y + 12 * Fxh * G * S * x17 * x2 * x67 + 12 * Fxh * G * S * x17 * x2 * x74 * xi + 6 * Fxh * G * x121 * x17 * x2 * x64 + 6 * Fxh * G * x135 * x17 * x2 * x62 + 6 * Fxh * G * x17 * x2 * x22 + 12 * Fxh * G * x17 * x2 * x24 * xi + 4 * Fxh * S * x18 * x2 * x22 + 8 * Fxh * S * x18 * x2 * x24 * xi + 8 * Fxh * x18 * x2 * x62 * x67 - Fxh * x241 - Fxh * x269 - Fxh * x270 + 4 * Fxp * Fy * S * x18 * x2 * x24 + 4 * Fxp * Fy * S * x18 * x2 * x64 * xi + 2 * Fxp * Fy * x0 * x18 * x2 + 2 * Fxp * Fy * x121 * x18 * x2 * x22 + 4 * Fxp * Fy * x18 * x2 * x23 * xi + 2 * Fxp * Fy * x18 * x2 * x62 * x74 + 4 * Fxt * Gp * S * x18 * x256 * x64 + 4 * Fxt * Gp * S * x18 * x256 * x67 * xi + 2 * Fxt * Gp * x121 * x18 * x24 * x256 + 4 * Fxt * Gp * x18 * x22 * x256 * xi + 2 * Fxt * Gp * x18 * x23 * x256 + 2 * Fxt * Gp * x18 * x256 * x62 * x66 + 4 * Fxt * S * Sp * x17 * x256 * x64 + 4 * Fxt * Sp * x17 * x22 * x256 * xi + 4 * Fxt * Sp * x17 * x23 * x256 + 4 * Fxt * x0 * x17 * x256 * xi + 4 * Fxt * x10 * x17 * x256 + 24 * Fy * G * Gt * S * x135 * x18 * x2 + 24 * Fy * G * Gt * S * x166 * x18 * x2 * xi + 12 * Fy * G * Gt * x121 * x18 * x2 * x66 + 12 * Fy * G * Gt * x18 * x2 * x218 * x62 + 12 * Fy * G * Gt * x18 * x2 * x67 + 24 * Fy * G * Gt * x18 * x2 * x74 * xi + 12 * Fy * G * S * St * x135 * x17 * x2 + 12 * Fy * G * St * x17 * x2 * x67 + 12 * Fy * G * St * x17 * x2 * x74 * xi + 4 * Fy * Gp * S * x18 * x67 * xi_y + 4 * Fy * Gp * S * x18 * x74 * xi * xi_y + 2 * Fy * Gp * x121 * x18 * x64 * xi_y + 2 * Fy * Gp * x135 * x18 * x62 * xi_y + 2 * Fy * Gp * x18 * x22 * xi_y + 4 * Fy * Gp * x18 * x24 * xi * xi_y + 28 * Fy * Gt * S * x17 * x2 * x67 + 32 * Fy * Gt * S * x17 * x2 * x74 * xi + 12 * Fy * Gt * x121 * x17 * x2 * x64 + 20 * Fy * Gt * x135 * x17 * x2 * x62 + 8 * Fy * Gt * x17 * x2 * x22 + 20 * Fy * Gt * x17 * x2 * x24 * xi + 4 * Fy * S * Sp * x17 * x67 * xi_y + 4 * Fy * Sp * x17 * x22 * xi_y + 4 * Fy * Sp * x17 * x24 * xi * xi_y + 24 * Fy * St * x18 * x2 * x22 + 16 * Fy * St * x18 * x2 * x24 * xi + 4 * Fy * x0 * x17 * xi_y - Fy * x134 * x204 + 4 * Fy * x17 * x23 * xi * xi_y - Fy * x204 * x304 * x333 - Fy * x55 * x93 + 4 * Fyh * Gp * S * x18 * x64 + 4 * Fyh * Gp * S * x18 * x67 * xi + 2 * Fyh * Gp * x121 * x18 * x24 + 4 * Fyh * Gp * x18 * x22 * xi + 2 * Fyh * Gp * x18 * x23 + 2 * Fyh * Gp * x18 * x62 * x66 + 4 * Fyh * S * Sp * x17 * x64 + 4 * Fyh * Sp * x17 * x22 * xi + 4 * Fyh * Sp * x17 * x23 + 4 * Fyh * x0 * x17 * xi + 4 * Fyh * x10 * x17 - Fyh * x124 - Fyh * x94 + 12 * Fyt * G * S * x17 * x2 * x67 + 12 * Fyt * G * S * x17 * x2 * x74 * xi + 6 * Fyt * G * x121 * x17 * x2 * x64 + 6 * Fyt * G * x135 * x17 * x2 * x62 + 6 * Fyt * G * x17 * x2 * x22 + 12 * Fyt * G * x17 * x2 * x24 * xi + 4 * Fyt * S * x18 * x2 * x22 + 8 * Fyt * S * x18 * x2 * x24 * xi + 8 * Fyt * x18 * x2 * x62 * x67 - Fyt * x241 - Fyt * x269 - Fyt * x270 + 24 * G * Gt * S * x135 * x18 * x2 * xi * xi_y + 24 * G * Gt * S * x18 * x2 * x66 * xi_y + 12 * G * Gt * x121 * x18 * x2 * x74 * xi_y + 12 * G * Gt * x18 * x187 * x2 * x62 * xi_y + 12 * G * Gt * x18 * x2 * x64 * xi_y + 24 * G * Gt * x18 * x2 * x67 * xi * xi_y + 12 * G * S * St * x17 * x2 * x66 * xi_y + 24 * G * S * x17 * x2 * x64 * xi_xy + 24 * G * S * x17 * x2 * x67 * xi * xi_xy + 12 * G * St * x17 * x2 * x64 * xi_y + 12 * G * St * x17 * x2 * x67 * xi * xi_y + 12 * G * x121 * x17 * x2 * x24 * xi_xy + 24 * G * x17 * x2 * x22 * xi * xi_xy + 12 * G * x17 * x2 * x23 * xi_xy + 12 * G * x17 * x2 * x62 * x66 * xi_xy - G * x342 * x68 + 8 * Gc * S * x17 * x2 * x24 + 8 * Gc * S * x17 * x2 * x64 * xi + 4 * Gc * x0 * x17 * x2 + 4 * Gc * x121 * x17 * x2 * x22 + 8 * Gc * x17 * x2 * x23 * xi + 4 * Gc * x17 * x2 * x62 * x74 + 24 * Gd * Gp * S * x121 * x2 * x67 + 8 * Gd * Gp * S * x129 * x2 * x74 + 8 * Gd * Gp * S * x2 * x24 + 24 * Gd * Gp * S * x2 * x64 * xi + 2 * Gd * Gp * x0 * x2 + 12 * Gd * Gp * x121 * x135 * x2 * x62 + 12 * Gd * Gp * x121 * x2 * x22 + 8 * Gd * Gp * x129 * x2 * x24 + 2 * Gd * Gp * x130 * x2 * x64 + 8 * Gd * Gp * x166 * x2 * x63 + 8 * Gd * Gp * x187 * x2 * x63 * xi + 2 * Gd * Gp * x199 * x2 * x65 + 8 * Gd * Gp * x2 * x23 * xi + 24 * Gd * Gp * x2 * x62 * x66 * xi + 12 * Gd * Gp * x2 * x62 * x74 + 8 * Gh * Gt * S * x18 * x2 * x66 * xi + 8 * Gh * Gt * S * x18 * x2 * x74 + 4 * Gh * Gt * x121 * x18 * x2 * x67 + 4 * Gh * Gt * x166 * x18 * x2 * x62 + 4 * Gh * Gt * x18 * x2 * x24 + 8 * Gh * Gt * x18 * x2 * x64 * xi + 4 * Gh * S * St * x17 * x2 * x74 + 4 * Gh * St * x17 * x2 * x24 + 4 * Gh * St * x17 * x2 * x64 * xi + 4 * Gp * S * x18 * x24 * xi_yy + 4 * Gp * S * x18 * x256 * x326 * x66 * xi + 4 * Gp * S * x18 * x256 * x326 * x74 + 4 * Gp * S * x18 * x59 * x66 * xi + 4 * Gp * S * x18 * x59 * x74 + 4 * Gp * S * x18 * x64 * xi * xi_yy + 2 * Gp * x0 * x18 * xi_yy + 2 * Gp * x121 * x18 * x22 * xi_yy + 2 * Gp * x121 * x18 * x256 * x326 * x67 + 2 * Gp * x121 * x18 * x59 * x67 + 2 * Gp * x166 * x18 * x256 * x326 * x62 + 2 * Gp * x166 * x18 * x59 * x62 + 4 * Gp * x18 * x23 * xi * xi_yy + 2 * Gp * x18 * x24 * x256 * x326 + 2 * Gp * x18 * x24 * x59 + 4 * Gp * x18 * x256 * x326 * x64 * xi + 4 * Gp * x18 * x59 * x64 * xi + 2 * Gp * x18 * x62 * x74 * xi_yy + 4 * Gt * S * Sh * x17 * x2 * x74 + 32 * Gt * S * x17 * x2 * x64 * xi_y + 32 * Gt * S * x17 * x2 * x67 * xi * xi_y + 4 * Gt * Sh * x17 * x2 * x24 + 4 * Gt * Sh * x17 * x2 * x64 * xi + 12 * Gt * x121 * x17 * x2 * x24 * xi_y + 24 * Gt * x17 * x2 * x22 * xi * xi_y + 12 * Gt * x17 * x2 * x23 * xi_y + 20 * Gt * x17 * x2 * x62 * x66 * xi_y - Gt * x187 * x236 * x283 * xi_y + 8 * S * Sc * x18 * x2 * x24 + 16 * S * Sd * x121 * x2 * x24 + 16 * S * Sd * x2 * x22 * xi + 4 * S * Sp * x17 * x24 * xi_yy + 4 * S * Sp * x17 * x256 * x326 * x74 + 4 * S * Sp * x17 * x59 * x74 + 16 * S * u * x2 + 72 * S * x0 * x121 * x2 + 56 * S * x10 * x2 * xi - S * x100 + 40 * S * x129 * x2 * x23 + 8 * S * x130 * x2 * x22 - S * x158 * x206 + 2 * S * x17 * x22 * x256 * x327 + 2 * S * x17 * x22 * x60 + 2 * S * x17 * x24 * x256 * x327 * xi + 2 * S * x17 * x24 * x60 * xi + 16 * S * x18 * x2 * x22 * xi * xi_xy + 16 * S * x18 * x2 * x23 * xi_xy - S * x22 * x240 * x273 - S * x25 * x52 - S * x315 * x67 - S * x361 + 8 * Sc * x0 * x18 * x2 + 8 * Sc * x18 * x2 * x23 * xi + 24 * Sd * x2 * x62 * x64 + 32 * Sd * x2 * x62 * x67 * xi + 16 * Sd * x2 * x63 * x66 + 4 * Sh ^ 2 * x17 * x24 - Sh * x107 * x108 - Sh * x107 * x111 - Sh * x33 * x44 - Sh * x35 * x36 - Sh * x73 * x83 + 4 * Sp * x0 * x17 * xi_yy + 4 * Sp * x17 * x23 * xi * xi_yy + 4 * Sp * x17 * x24 * x256 * x326 + 4 * Sp * x17 * x24 * x59 + 4 * Sp * x17 * x256 * x326 * x64 * xi + 4 * Sp * x17 * x59 * x64 * xi - Sp * x263 * x323 * xi_xy - Sp * x53 + 4 * St ^ 2 * x17 * x24 * x256 + 16 * St * x18 * x2 * x22 * xi * xi_y + 16 * St * x18 * x2 * x23 * xi_y - St * x24 * x272 * x276 - u * x14 - x0 * x129 * x15 + 2 * x0 * x17 * x256 * x327 * xi + 2 * x0 * x17 * x60 * xi - x0 * x18 * x51 * xi_xy - x0 * x20 * x259 - x0 * x21 + 2 * x0 * x256 * x6 * xi_xx * (x17 * (-3 * B * x414 + Bp * x6 - 4 * S * u + 2 * Sp) + x413 * (Gp - x412)) - x10 * x114 * x121 + x10 * x17 * x256 * x327 + x10 * x17 * x60 - x10 * x54 * xi - x100 * x129 - x101 * x102 - x101 * x106 - 12 * x101 * x264 * x272 - x101 * x319 - x101 * x321 - x101 * x362 - x101 * x55 * x75 - x102 * x188 - x102 * x235 - x102 * x238 - x104 * x73 - x105 * x106 - x105 * x110 - x105 * x319 - x105 * x321 - x106 * x161 - 144 * x108 * x149 - x108 * x178 * x197 - x108 * x192 * x391 - 72 * x108 * x192 * x99 - x108 * x307 * x308 - x108 * x311 - 8 * x108 * x340 - x108 * x344 * x72 - x108 * x345 - x108 * x349 * x50 - 20 * x109 * x169 - x109 * x35 * xi - x11 * x116 * x308 * x66 * xi - x11 * x212 * x349 - x11 * x214 * x309 - x11 * x275 * x29 - x11 * x309 * x57 - x110 * x161 - x111 * x116 * x322 - 96 * x111 * x149 - 56 * x111 * x151 - x111 * x307 * x308 - 8 * x111 * x340 - x111 * x345 - x111 * x349 * x50 - x112 * x113 - 4 * x112 * x302 - x112 * x48 - x113 * x47 - x113 * x49 - x115 * x128 - x115 * x49 - x116 * x125 - x116 * x133 * x140 * x66 - x116 * x167 * x309 * x50 - x116 * x206 * x347 - x116 * x32 * xi * xi_yy - x117 * x124 - 24 * x117 * x176 * x83 - x118 * x120 - x118 * x142 * xi - x118 * x341 - x118 * x56 - x118 * x61 - x119 * x120 - x119 * x341 - x119 * x56 - x119 * x61 - x12 * x131 - x12 * x232 * x62 - x12 - x120 * x168 - x120 * x249 - x120 * x58 + x121 * x17 * x23 * x256 * x327 + x121 * x17 * x23 * x60 + x121 * x17 * x24 * x256 * x326 + x121 * x17 * x24 * x59 + 68 * x121 * x2 * x24 * x62 + 8 * x121 * x2 * x63 * x74 - x121 * x263 * x365 - x121 * x53 - x121 * x58 * x77 - x122 * x123 - 18 * x122 * x83 - x123 * x354 - x124 * x284 - x125 * x83 - x126 * x127 - 24 * x126 * x164 - x126 * x171 - x126 * x172 - x126 * x201 - x126 * x356 - x126 * x359 - x126 * x366 - x126 * x388 - x127 * x128 - x128 * x171 - x128 * x172 - x128 * x201 - x128 * x356 - x128 * x359 - x128 * x366 - x128 * x388 + 16 * x129 * x2 * x62 * x64 - x129 * x200 - x129 * x361 - x130 * x132 - x133 * x145 * x370 - x134 * x137 - x134 * x139 - x135 * x236 * x396 - x136 * x141 - x136 * x336 - x136 * x367 - x137 * x346 - x138 * x147 - x138 * x148 - x138 * x313 * x348 - x138 * x337 * x373 - x138 * x371 - x139 * x346 - x141 * x68 - x142 * x143 - x142 * x78 * xi - x144 * x257 - x144 * x28 - x144 * x298 - x144 * x81 - x145 * x147 - x145 * x148 - x145 * x299 * x84 - x145 * x337 * x373 - x145 * x376 - x147 * x205 - 156 * x149 * x157 * xi - 102 * x149 * x212 - 54 * x149 * x232 - x15 * x4 - x150 * x34 * xi - x152 * x343 - x152 * x95 - x153 * x156 - x153 * x237 - x155 * x160 * x83 - x156 * x391 - x157 * x158 - x157 * x160 * x83 - x157 * x178 * x353 - x157 * x178 * x59 - x157 * x323 * x340 - x159 * x347 - 24 * x159 * x82 * xi_y - x162 * x277 - x163 * x165 - x163 * x278 - x165 * x174 - x168 * x387 - x168 * x390 - x168 * x56 - x168 * x61 - x169 * x17 * x26 - x169 * x43 + 6 * x17 * x22 * x256 * x326 * xi + 6 * x17 * x22 * x59 * xi + 9 * x17 * x23 * x256 * x326 + 9 * x17 * x23 * x59 + x17 * x256 * x327 * x62 * x67 - x17 * x299 * x35 - x17 * x52 * x62 * x64 + x17 * x60 * x62 * x67 - x173 * x277 - x174 * x278 - x175 * x183 - x175 * x228 - x175 * x383 - x175 * x389 - x175 * x90 - x176 * x183 - x176 * x228 - x176 * x383 - x176 * x389 - x177 * x34 - x178 * x196 * x67 - x179 * x180 - x179 * x206 * x234 - x179 * x250 - x179 * x251 + 16 * x18 * x2 * x62 * x64 * xi_xy - x18 * x200 * xi_xy - x18 * x236 * x365 * x67 - x18 * x258 * x9 - x180 * x378 - x181 * x284 * x333 - x181 * x351 - x182 * x382 - x182 * x59 * xi - x186 * x360 - x186 * x99 - x187 * x378 * x402 - x188 * x362 - x189 * x191 - x189 * x195 - x19 * x9 - x191 * x399 - x192 * x396 * x64 - x193 * x194 - x193 * x281 - x194 * x400 - x195 * x399 - 24 * x196 * x93 - 48 * x197 * x87 - x198 * x43 - x198 * x46 + 88 * x2 * x22 * x62 * xi + 36 * x2 * x23 * x62 + 24 * x2 * x63 * x64 + 32 * x2 * x63 * x67 * xi + 4 * x2 * x65 * x66 - x2 * x9 * xi_x * (x17 * (x2 * (Fx * (78 * B * x112 + 51 * B * x126 + 27 * B * x128 - Bp * x10 * x419 - 2 * Sp * x10 * x6 + 21 * x1 + x108 * x163 + x108 * x174 + x111 * x163 + x111 * x174 + 12 * x131 + x157 * x244 + x157 * x253 + x207 * x212 + x207 * x232 + x207 * x24 + x212 * x215 + x215 * x232 + 34 * x3 + x379 * x47 + x381 * x49 - x415 + x425 + x49 * x89 - 2) + 2 * u * (Bt * x427 * x6 + 3 * Gt * x426 + St * (4 * x4 + x417 + 4))) - x414 * (Gh * (10 * x3 + 6 * x4 + x417 + 6) - u * (-Fy * x412 * x429 + Fy * x428 + x10 * x411 * (Bh * x3 + Bh * x4 + Bh - 2 * Sh)))) + x18 * (u * (-6 * Gh * x426 + x2 * x6 * (u * (Fx * x412 * (12 * x416 + x429) - Fx * x428 + x10 * x37 * (Bt * x3 + Bt * x4 + Bt + St)) + x260 * x427) - x276 * x5) + x30 * (Sp * x10 * x6 - 3 * x131 * x423 - x3 * (14 * x4 + x425 * x5 + 17) - x5 * (x422 * x5 - 1))) - x271 * xi_y * (x103 * x418 + x18 * x424)) - x202 * x310 - x202 * x62 - x204 * x367 - x205 * x337 * x373 - x205 * x371 - x205 * x376 - x207 * x209 - x207 * x210 - x207 * x289 - x209 * x215 - x209 * x216 - x21 * x47 - x21 * x49 - x210 * x215 - x211 * x213 - x211 * x214 - x212 * x386 - x213 * x217 - x213 * x262 - x214 * x217 - x214 * x262 - x215 * x289 - x216 * x397 * x66 - 12 * x218 * x254 * x337 - x218 * x382 * x402 * x83 - x219 * x220 - x219 * x401 - x219 * x404 - x22 * x30 * x76 - x22 * x315 - x220 * x221 - x221 * x337 * x84 - x221 * x401 - x221 * x404 - x222 * x223 - x222 * x239 - x223 * x224 - x224 * x239 - x224 * x305 * xi_y - x225 * x226 - x225 * x227 - x225 * x255 - x226 * x233 - x227 * x233 - x229 * x230 - x229 * x369 - x23 * x43 - x23 * x46 - x230 * x231 - x231 * x369 - x232 * x386 - x233 * x255 - x234 * x403 - x235 * x362 - x237 * x391 - x238 * x362 - x24 * x311 - x240 * x266 * x49 - x240 * x273 * x47 * xi - x242 * x243 - x242 * x246 - x242 * x406 - x243 * x252 - x244 * x245 - x244 * x248 - x244 * x279 - x244 * x280 - x244 * x407 - x244 * x408 - x245 * x253 - x246 * x252 - 48 * x246 * x254 - x248 * x253 - x249 * x387 - x249 * x390 - x25 * x26 - x25 * x322 * xi - x25 * x55 * xi_y - x251 * x353 * x72 - x251 * x378 - x252 * x406 - x253 * x279 - x253 * x280 - x253 * x407 - x253 * x408 - x256 * x312 * xi_x ^ 2 * (x17 * (B * x418 + x162 * x419 * x420 + x424) + x37 * x413 * (4 * x3 + x415 + x417 + 2)) - x257 * x27 - x257 * x79 - x257 * x80 - x259 * x79 - x259 * x80 - x261 * x306 * x68 - x262 * x57 - 18 * x263 * x264 - x265 * x266 - x267 * x268 - x267 * x318 - x267 * x363 - x267 * x364 - x268 * x274 - x27 * x28 - x27 * x298 - x270 * x273 - x271 * x272 * x273 - x274 * x318 - x274 * x363 - x274 * x364 - x28 * x79 - x28 * x80 - x281 * x400 - x282 * x304 * x69 - x283 * x285 - x284 * x333 * x93 - x284 * x94 - x285 * x38 - x285 * x40 - x286 * x287 - x286 * x288 - x287 * x300 - x287 * x409 - x287 * x410 - x288 * x300 - x288 * x377 * x378 - x288 * x409 - x288 * x410 - x29 * x31 - x290 * x291 - x290 * x397 - x291 * x301 - x293 * x294 - x293 * x296 - x294 * x374 - x294 * x384 - x294 * x385 - x295 * x297 - x296 * x339 * x91 - x296 * x374 - x296 * x384 - x296 * x385 - x297 * x303 - x298 * x79 - x298 * x80 - x299 * x85 - 16 * x299 * x93 - x3 * x51 - x301 * x397 - 4 * x302 * x47 - 4 * x302 * x49 - x305 * x64 * xi_y - x312 * x313 * x314 - x315 * x91 - x316 * x317 - x316 * x352 - x317 * x320 - x32 * x34 - x320 * x352 - x322 * x92 - x324 * x325 - x324 * x357 - x324 * x358 - x325 * x328 - x328 * x357 - x328 * x358 - x331 * x41 - x332 * x334 - x332 * x338 - x334 * x335 - x335 * x338 - x336 * x68 - x337 * x85 - x34 * x67 * x89 - x343 * x96 - x344 * x87 - x351 * x93 - x355 * x72 - x355 * x83 - x368 * x369 - x368 * x375 - 12 * x372 * x72 - x372 * x84 - x377 * x378 * x66 - x379 * x380 * x66 - x379 * x403 - x38 * x39 - x381 * x382 * x68 - x387 * x398 - x39 * x40 - x390 * x398 - x392 * x393 - x392 * x394 - x393 * x395 - x394 * x395 - x40 * x64 * x77 - x405 * x72 - x405 * x83 - x41 * x42 - x47 * x48 - x48 * x49 - x56 * x58 - x58 * x61 - x69 * x71 - 12 * x69 * x82 - x79 * x81 - x80 * x81 - x85 * x86 - 42 * x87 * x88 - 12 * x87 * x95 - 68 * x88 * x92 - x90 * x91 - x95 * x96

    nothing
end
