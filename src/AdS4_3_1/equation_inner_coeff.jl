
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
x21 = 16 * Sp
x22 = Bs * x17
x23 = x16 * x5
x24 = 2 * x23
x25 = sinh(x15)
x26 = x0 * x25
x27 = 2 * x26
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
x43 = Bh * x42
x44 = Gh * x25
x45 = 4 * x44
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
x59 = 6 * G
x60 = x25 * x59
x61 = Fyh * x32
x62 = 2 * x35
x63 = Gp * x25
x64 = Fyh * x63
x65 = 8 * S
x66 = Fyp * x44
x67 = Fyp * xi_y
x68 = S * x35
x69 = x28 * x68
x70 = 4 * xi
x71 = x23 * x70
x72 = x60 * xi_yy
x73 = x42 * x45
x74 = x44 * x48
x75 = S * x42
x76 = Gs * x25
x77 = 4 * x76
x78 = Sp * x8
x79 = x7 * xi
x80 = 16 * xi_y
x81 = Sh * x80
x82 = x13 * x16
x83 = Bh ^ 2
x84 = 2 * x42
x85 = x16 * x84
x86 = Fy ^ 2
x87 = 3 * x36
x88 = Gh ^ 2
x89 = S ^ 2
x90 = x7 * x89
x91 = S ^ 3
x92 = x7 * x91
x93 = u ^ 7
x94 = 48 * x93
x95 = u ^ 10
x96 = 12 * x7
x97 = x95 * x96
x98 = x16 * x38
x99 = u ^ 8
x100 = 12 * B
x101 = x100 * x99
x102 = x16 * xi_y
x103 = Bh * x102
x104 = x100 * x93
x105 = x42 * x56
x106 = x16 * x31
x107 = Fy * x44
x108 = x16 * x58
x109 = Fy * x42
x110 = 48 * x109
x111 = B * x102
x112 = Fyh * x16
x113 = S * x99
x114 = x100 * x113
x115 = x42 * xi
x116 = x100 * x112
x117 = x31 * x33
x118 = x44 * xi_y
x119 = x16 * xi_yy
x120 = Sh * x102
x121 = x33 * xi
x122 = x121 * xi_yy
x123 = 12 * G
x124 = x123 * x25
x125 = x124 * x99
x126 = 32 * x113
x127 = 24 * x115
x128 = S * x93
x129 = Bh * Fyp
x130 = x129 * x28
x131 = x32 * xi
x132 = Bh * xi_y
x133 = x124 * x93
x134 = Bh * x44
x135 = u ^ 9
x136 = S * x135
x137 = 8 * x136
x138 = x93 * xi
x139 = 8 * x138
x140 = x136 * x28
x141 = S * x16
x142 = x141 * x93
x143 = 32 * x142
x144 = x138 * x28
x145 = x33 * xi_y
x146 = 2 * Fy
x147 = Bp * x146
x148 = Fyh * x128
x149 = Bp * x28
x150 = x149 * xi
x151 = Bp * xi_yy
x152 = S * x144
x153 = 8 * x56
x154 = x16 * x75
x155 = x36 * xi
x156 = Gh * x123
x157 = Fy * x16
x158 = x157 * x99
x159 = G * x25
x160 = x63 * xi_y
x161 = x146 * x160
x162 = Fy * x102
x163 = 76 * x75
x164 = 16 * x115
x165 = Fy * xi_y
x166 = Fyh * x124
x167 = 4 * x128
x168 = x32 * x70
x169 = x75 * xi
x170 = 12 * x169
x171 = x13 * xi
x172 = x32 * x60
x173 = Fyp * x45
x174 = S * x32
x175 = x174 * x28
x176 = x67 * xi
x177 = x175 * xi
x178 = x102 * x93
x179 = S * xi_yy
x180 = Sh * xi_y
x181 = x124 * xi_yy
x182 = Sh * x45
x183 = x63 * xi_yy
x184 = 4 * x75
x185 = 4 * x54
x186 = x128 * x70
x187 = Sp * x20
x188 = x42 * x82
x189 = 24 * B
x190 = x16 * x189
x191 = x86 * x93
x192 = xi_y ^ 2
x193 = x192 * x32
x194 = Bp * x85
x195 = x135 * x89
x196 = xi ^ 2
x197 = x196 * x32
x198 = 24 * x159
x199 = x63 * x86
x200 = 54 * x141
x201 = 6 * x121
x202 = x89 * x99
x203 = x16 * x202
x204 = 10 * x203
x205 = x196 * x36
x206 = 2 * x205
x207 = x17 * x174
x208 = x0 * x17
x209 = x30 * xi
x210 = x89 * x93
x211 = Fyph * x17
x212 = x0 * x196
x213 = 2 * x195
x214 = 2 * x197
x215 = x196 * x7
x216 = xi ^ 3
x217 = x216 * x7
x218 = xi ^ 4
x219 = x174 * x8
x220 = 24 * x141
x221 = x32 * x89
x222 = x99 * xi
x223 = x187 * x196
x224 = x14 * x35
x225 = u ^ 11
x226 = S * x225
x227 = x190 * x38
x228 = x135 * xi
x229 = S * x95
x230 = x132 * x190
x231 = x16 * x56
x232 = x100 * x136
x233 = x138 * x231
x234 = x159 * x165
x235 = 72 * B
x236 = x135 * x235
x237 = x107 * x189
x238 = x100 * x141
x239 = x225 * x238
x240 = Fy * x111
x241 = 156 * x136
x242 = 108 * x240
x243 = x136 * xi
x244 = x16 * x67
x245 = x100 * x115
x246 = x118 * x189
x247 = x238 * x95
x248 = x198 * x38
x249 = 36 * x243
x250 = x113 * xi
x251 = x132 * x198
x252 = x95 * xi
x253 = x252 * x65
x254 = x126 * xi
x255 = x113 * x165
x256 = x165 * xi
x257 = x113 * x150
x258 = x124 * x56
x259 = Fy * x225
x260 = G * Gh
x261 = x220 * x260
x262 = 24 * x260
x263 = x157 * x228
x264 = x124 * x58
x265 = 108 * x234
x266 = x63 * x70
x267 = x138 * x141
x268 = x113 * x70
x269 = x124 * x67
x270 = x261 * xi_y
x271 = x102 * x222
x272 = x124 * x180
x273 = 2 * u
x274 = x25 * x7
x275 = x273 * x274
x276 = 36 * B
x277 = x159 * x86
x278 = x277 * x95
x279 = x16 * x86
x280 = x279 * x95
x281 = B * S
x282 = 84 * x281
x283 = B * x279
x284 = 60 * x222
x285 = x225 * x89
x286 = x112 * x31
x287 = x196 * x93
x288 = x159 * x192
x289 = x276 * x288
x290 = x16 * x192
x291 = B * x290
x292 = 72 * x113
x293 = x89 * x95
x294 = x119 * x31
x295 = 48 * x291
x296 = x196 * x42
x297 = B ^ 2
x298 = 36 * x297
x299 = x135 * x162
x300 = 22 * x285
x301 = 14 * x287
x302 = x17 * x293
x303 = u ^ 12
x304 = x303 * x89
x305 = Bh * x45
x306 = x196 * x99
x307 = 20 * x293
x308 = x16 * x296
x309 = x28 * xi
x310 = x229 * x309
x311 = Bp * x140
x312 = x191 * xi
x313 = x293 * x52
x314 = x194 * x196
x315 = x52 * xi_yy
x316 = x28 * x56
x317 = G ^ 2
x318 = 36 * x317
x319 = 4 * x136
x320 = x136 * x82
x321 = 42 * x250
x322 = Fyh * x60
x323 = 2 * x293
x324 = x196 * x84
x325 = S * x85
x326 = 48 * x288
x327 = x75 * x78
x328 = x303 * xi
x329 = S * x328
x330 = x226 * xi
x331 = 144 * x234
x332 = x281 * x331
x333 = x162 * x252
x334 = x100 * x243
x335 = x229 * xi
x336 = 18 * x297
x337 = x290 * x99
x338 = x17 * x83
x339 = exp(2 * x6)
x340 = Fxpt * x339
x341 = 18 * x317
x342 = 21 * x89
x343 = 5 * x296
x344 = x17 * x88
x345 = x90 * x93
x346 = 16 * x345
x347 = x135 * x196
x348 = u ^ 14
x349 = x348 * x89
x350 = x100 * x98
x351 = x196 * x95
x352 = u ^ 13
x353 = x352 * x89
x354 = x100 * x103
x355 = x106 * x56
x356 = x100 * x107
x357 = 60 * x306
x358 = 72 * x281
x359 = x277 * x358
x360 = x225 * x277
x361 = x235 * xi
x362 = x225 * x279
x363 = 96 * xi
x364 = x281 * x363
x365 = x106 * x67
x366 = x288 * x358
x367 = x100 * x118
x368 = 72 * x243
x369 = 72 * x297
x370 = x141 * x303
x371 = x165 * x370
x372 = x124 * x38
x373 = x124 * x132
x374 = x102 * x147
x375 = x150 * x229
x376 = x56 * x60
x377 = x156 * x157
x378 = 72 * x317
x379 = x60 * x67
x380 = x102 * x156
x381 = 2 * x5
x382 = x274 * x381
x383 = x29 * x7
x384 = Sc * x8
x385 = 60 * x353
x386 = 36 * x347
x387 = x192 * x287
x388 = S * x298
x389 = x279 * x352
x390 = x298 * xi
x391 = x225 * x290
x392 = x135 * x290
x393 = x52 * x86
x394 = S * x318
x395 = x318 * xi
x396 = 2 * x199
x397 = u ^ 15
x398 = x397 * x89
x399 = x234 * x235
x400 = x196 * x225
x401 = x348 * xi
x402 = x141 * x256 * x352
x403 = Bb * x339
x404 = Fxt * x339
x405 = Gb * x339
x406 = x29 * x339
x407 = 2 * Gt
x408 = x407 * x7
x409 = Gh * x7
x410 = Bt * x409
x411 = x274 * x35
x412 = Fx * Fy
x413 = Fx * x7
x414 = 4 * Fyp
x415 = x26 * x414
x416 = Fx * x32
x417 = 10 * x416
x418 = x16 * x417
x419 = x20 * x25
x420 = Fxh * x7
x421 = x33 * x59
x422 = Gp * x420
x423 = 2 * x36
x424 = Fxh * x25
x425 = Fy * x7
x426 = 4 * Fxp
x427 = Fxp * x274
x428 = Fxp * x36
x429 = Fxp * xi_y
x430 = x274 * x68
x431 = 4 * x430
x432 = Fxph * x25
x433 = 4 * x79
x434 = x433 * x5
x435 = x28 * x425
x436 = Gt * x7
x437 = Fy * x436
x438 = Spt * x8
x439 = Fy * x25
x440 = x438 * x439
x441 = Fypt * x25
x442 = Fyt * x7
x443 = Gp * x442
x444 = Fyt * x25
x445 = x36 * x96
x446 = G * xi_xy
x447 = Gc * x8
x448 = St * x409
x449 = Gp * xi_xy
x450 = Gpt * xi_y
x451 = Sh * x436
x452 = Gt * xi_y
x453 = x25 * x384
x454 = 16 * xi_xy
x455 = 8 * Sh
x456 = x438 * xi_y
x457 = St * x80
x458 = u ^ 16
x459 = x458 * x89
x460 = x276 * x277
x461 = x196 * x303
x462 = x162 * x298
x463 = x388 * xi
x464 = x279 * x348
x465 = x290 * x303
x466 = x162 * x318
x467 = x394 * xi
x468 = Fxp ^ 2 * x339
x469 = x279 * x336
x470 = x290 * x349
x471 = x336 * x95
x472 = x196 * x290
x473 = x36 * x404
x474 = Bt * x339
x475 = 2 * x339
x476 = x25 * x339
x477 = Gt * x476
x478 = Bt * x477
x479 = 4 * x42
x480 = x339 * x46
x481 = Bt * St
x482 = Fx * Fxp
x483 = x339 * x416
x484 = St * x16
x485 = Fxp * x476
x486 = x407 * x485
x487 = x404 * x63
x488 = x279 * x341
x489 = x341 * x95
x490 = x25 * x405
x491 = St * x477
x492 = Sb * x339
x493 = Fx * x99
x494 = x16 * x99
x495 = Bh * x59
x496 = x28 * x79
x497 = x496 * x93
x498 = Gt * x497
x499 = Bt * x59
x500 = Bt * Gh
x501 = G * x413
x502 = Fy * x46
x503 = Gp * x413
x504 = Gpp * x416
x505 = x274 * x412
x506 = Fx * x25
x507 = x109 * x506
x508 = x416 * x439
x509 = Spp * x8
x510 = x96 * xi
x511 = x16 * x59
x512 = x75 * x8
x513 = Fyp * x506
x514 = x54 * x8
x515 = G * x96
x516 = Fx * x16
x517 = x516 * x99
x518 = Sh * x515
x519 = x102 * x501
x520 = x16 * x409
x521 = 32 * Fx
x522 = x113 * x521
x523 = x115 * x20
x524 = Gh * x516
x525 = 2 * Fx
x526 = Gp * x525
x527 = x526 * x7
x528 = Gpp * x28
x529 = x413 * x528
x530 = x274 * xi_y
x531 = Fx * x530
x532 = 16 * x79
x533 = x25 * xi_y
x534 = x416 * x533
x535 = x506 * xi_y
x536 = x35 * x535
x537 = x16 * x515
x538 = Fxh * x113
x539 = x115 * x537
x540 = x128 * x28
x541 = Fxh * Gp
x542 = x32 * x496
x543 = x14 * x424
x544 = x169 * x96
x545 = Fxp * x7
x546 = Fxp * x439
x547 = x174 * x427
x548 = Fxp * x28
x549 = Gh * x79
x550 = x548 * x549
x551 = 4 * xi_y
x552 = x174 * x433
x553 = x159 * x96
x554 = Gt * x99
x555 = St * x515
x556 = Fy * Gpt
x557 = x556 * x8
x558 = x16 * x437
x559 = Gt * x157
x560 = Fyp * Gt
x561 = Fyt * x537
x562 = Fyt * Gp
x563 = x14 * x444
x564 = x452 * x553
x565 = x20 * x446
x566 = x136 * x8
x567 = Gt * x44
x568 = x567 * x8
x569 = x449 * x8
x570 = x450 * x8
x571 = x452 * x7
x572 = x20 * x452
x573 = x25 * xi_xy
x574 = x25 * x456
x575 = x25 * x79
x576 = Bt ^ 2
x577 = x339 * x85
x578 = Fx ^ 2 * x339
x579 = Gt ^ 2
x580 = x90 * x99
x581 = 10 * x580
x582 = x215 * x62
x583 = 2 * x345
x584 = x215 * x27
x585 = Gc * x28
x586 = x135 * x90
x587 = x215 * x32
x588 = x474 * x516
x589 = x339 * x482
x590 = x42 * x589
x591 = Fx * x101
x592 = x339 * x484
x593 = x16 * x404
x594 = x149 * x404
x595 = Fx * x125
x596 = x16 * x474
x597 = x474 * x548
x598 = x140 * x339
x599 = x144 * x339
x600 = 8 * x339
x601 = Fx * xi
x602 = Gt * x123
x603 = x339 * x602
x604 = St * x339
x605 = Fx * x477
x606 = Gt * x485
x607 = x124 * x404
x608 = Fx * x239
x609 = x100 * x135
x610 = x559 * x79
x611 = Gt * x102
x612 = x611 * x79
x613 = Fx * x225
x614 = x141 * x613
x615 = Bh * x515
x616 = x135 * x601
x617 = x229 * x496
x618 = Bt * x515
x619 = x141 * x259
x620 = Bt * x141
x621 = G * x97
x622 = x621 * xi_y
x623 = Fy * x95
x624 = 168 * x141
x625 = G * x79
x626 = 120 * x412
x627 = Gp * x412
x628 = x601 * x93
x629 = x157 * x8
x630 = Gpp * x8
x631 = x113 * x630
x632 = x16 * x412
x633 = x630 * x632
x634 = x25 * x412
x635 = x509 * x634
x636 = x634 * x79
x637 = x136 * x515
x638 = x601 * x8
x639 = Fx * x226
x640 = x20 * x44
x641 = G * x640
x642 = G * x135
x643 = x16 * x601
x644 = 108 * Fx
x645 = 36 * Fx
x646 = Fx * xi_y
x647 = x46 * x646
x648 = Gp * x79
x649 = x14 * x535
x650 = x535 * x79
x651 = Fxp * x157
x652 = S * x138
x653 = x184 * x575
x654 = Fxp * x102 * x515
x655 = Fy * Gt
x656 = x159 * x655
x657 = x20 * x656
x658 = x113 * x496
x659 = x159 * x572
x660 = x578 * x93
x661 = x17 * x403
x662 = x578 * x63
x663 = x468 * xi
x664 = x17 * x340
x665 = x303 * x90
x666 = Bh * x407
x667 = x17 * x500
x668 = x215 * x99
x669 = x317 * x645
x670 = x90 * x95
x671 = x215 * x25
x672 = x416 * x671
x673 = x225 * x90
x674 = 22 * x673
x675 = 14 * Fx
x676 = x287 * x675
x677 = x511 * x673
x678 = x287 * x511
x679 = x17 * x670
x680 = x215 * x85
x681 = Fxp * x25
x682 = Fxp * Gh
x683 = x25 * x429
x684 = x16 * x670
x685 = Fyp * x407
x686 = x16 * x215
x687 = x28 * x449
x688 = x28 * x450
x689 = x190 * x474
x690 = Fxp * x339
x691 = x643 * x690
x692 = x189 * x477
x693 = x198 * x474
x694 = x243 * x645
x695 = Fxp * x474
x696 = x601 * x690
x697 = Gt * x339
x698 = G * x220 * x697
x699 = G * x16
x700 = x124 * x604
x701 = x238 * x303
x702 = G * x510
x703 = x225 * xi_y
x704 = x141 * x621
x705 = x303 * x601
x706 = S * x705
x707 = x159 * x578
x708 = x276 * x707
x709 = x16 * x578
x710 = x709 * x95
x711 = B * x709
x712 = x106 * x404
x713 = 4 * x478
x714 = x310 * x339
x715 = x404 * x60
x716 = x348 * x90
x717 = x106 * x716
x718 = Fx * Gh
x719 = x215 * x95
x720 = x352 * x90
x721 = x31 * x611
x722 = x135 * x215
x723 = x495 * x516
x724 = x157 * x499
x725 = x102 * x499
x726 = S * x378
x727 = x352 * x726
x728 = x28 * x627
x729 = Fyp * x59
x730 = Fx * x44
x731 = G * x102
x732 = x102 * x673
x733 = x102 * x287
x734 = Fxp * x59
x735 = x17 * x339
x736 = x576 * x735
x737 = x579 * x735
x738 = x100 * x588
x739 = x106 * x589
x740 = x100 * x605
x741 = x358 * x707
x742 = x225 * x709
x743 = Fx * x124 * x474
x744 = x589 * x60
x745 = x516 * x603
x746 = x141 * x578
x747 = x352 * x746
x748 = x52 * x578
x749 = 2 * x662
x750 = x439 * x669
x751 = x348 * x746
x752 = x336 * x709
x753 = x341 * x709
x754 = 3 * G
x755 = u * x754
x756 = x25 * x4
x757 = 2 * x4
x758 = 3 * x6
x759 = x4 * x758
x760 = (G * x2 + G) ^ 2
x761 = 9 * x0
x762 = x317 * x42
x763 = 9 * x762 + 4
x764 = 6 * x3 * (3 * x762 + 2)
x765 = Gpp * x4
x766 = -Spp * x757
x767 = 9 * x32
x768 = 6 * x15
x769 = u * x196
x770 = x1 * xi
x771 = x762 * xi
x772 = 9 * x317
x773 = 5 * x1 + 3 * x2 + x759 + 3
x774 = 6 * x4 * x6
x775 = 9 * x1 + 5 * x2 + 4
x776 = u * x59
x777 = x525 * xi
x778 = Bp * Fx
x779 = Fx * x5
x780 = Bp * x777
x781 = 2 * Sp
x782 = B * Fx
x783 = x113 * x645
x784 = x297 * x645
x785 = Gt * x59
x786 = Fx * x336
x787 = Fx * x341
ABCS[1] = 0
ABCS[2] = -x4 ^ 3 * x9
ABCS[3] = -x10 * x9 * (-x11 + 3 * x12 + xi)
ABCS[4] = -B * x252 * x331 + Bh * Gt * x617 + Bh * x140 * x436 - Bh * x308 * x48 + Bh * x498 + Bh * x537 * x616 - Bp * x144 * x578 + Bp * x256 * x46 - Bs * x152 - Bs * x53 - Bs * x55 - Bt * Fxp * x196 * x577 - Bt * Fy * x370 * x702 - Bt * x428 * x475 + Fx * x164 * x592 - Fx * x549 * x701 - Fxh * x243 * x537 + Fxh * x382 - Fxh * x539 - Fxh * x677 + Fxp * Fy * x704 * xi - Fxp * Fyp * x653 - Fxp * x415 * x79 + Fxph * x275 + Fxph * x431 + Fxph * x584 - 4 * Fy * x113 * x160 - Fy * x13 * x145 + Fy * x261 * x328 + Fy * x287 * x529 - Fy * x553 * x554 + Fyh * x204 + Fyh * x206 - Fyh * x24 + Fyh * x257 + Fyh * x313 + Fyh * x314 + Fyp * x128 * x25 * x638 + Fyp * x36 * x408 + Fyp * x413 * x42 * x511 + Fyp * x436 * x540 + Fyp * x516 * x637 + Fyp * x537 * x628 - Fyp * x582 * x681 + Fyp * x601 * x704 - Fyph * x177 - Fyph * x18 - Fyph * x69 - Fyph * x71 + Fypt * x275 + Fypt * x431 + Fypt * x584 + Fyt * x382 - Fyt * x539 - Fyt * x677 - Gc * x383 - Gh * St * x497 + Gp * x141 * x623 * x638 - Gp * x27 * xi_yy + Gp * x496 * x538 + Gp * x628 * x629 + Gpp * x215 * x647 + Gpt * x102 * x512 + Gpt * x215 * x502 + Gpt * x250 * x629 + Gpt * x35 * x435 + Gs * x27 + Gt * x106 * x215 * x623 - Gt * x45 * x665 - S ^ 4 * x97 - S * x104 * x119 + 84 * S * x278 + S * x360 * x363 + Sb * x406 - Sh ^ 2 * x46 - Sh * x416 * x419 - Sh * x42 * x506 * x532 - Sh * x498 + Sh * x73 + Ss * x29 + Ss * x53 + Ss * x55 - St ^ 2 * x480 - St * x109 * x25 * x532 - St * x141 * x622 + St * x274 * x42 * x455 - St * x419 * x57 + x0 * x21 * x217 - x0 * x22 - x0 * x339 * x757 * xi_xx * (x16 * (-3 * B * u * x4 + Bp * x4 + 2 * Sp - 4 * x19) + x756 * (Gp - x755)) - 96 * x1 * x215 + x1 * x78 - x100 * x108 * x228 - x100 * x120 * x222 - x100 * x122 + x100 * x233 - x101 * x107 - x101 * x108 + x101 * x588 + x101 * x612 + x101 * x98 + x102 * x504 * x8 * xi + x102 * x601 * x631 + x103 * x104 - x103 * x254 - x103 * x307 - x104 * x118 - x104 * x120 - x104 * x691 + x105 * x106 - x105 * x60 - x106 * x409 * x493 + x106 * x437 * x99 - x106 * x590 - x106 * x718 * x719 + x107 * x126 + x107 * x127 + x107 * x249 + x107 * x300 + x107 * x301 + x108 * x164 - x109 * x266 * xi_y + x109 * x511 * x545 - x110 * x111 + x110 * x159 * xi_y - x112 * x114 + x112 * x170 + x113 * x166 - x113 * x266 * x404 - x113 * x269 + x113 * x28 * x503 * xi_y - x113 * x309 * x695 + x113 * x550 - x113 * x561 + x113 * x607 + x113 * x635 - 84 * x113 * x636 + x113 * x649 + x113 * x654 - x114 * x119 * xi + x114 * x244 + x114 * x593 - x115 * x116 + x115 * x166 - x115 * x269 - x115 * x295 + x115 * x326 + x115 * x607 + x115 * x633 + x115 * x635 + x115 * x649 + x115 * x654 - x116 * x243 + x117 * x404 + x117 * x67 + 32 * x118 * x128 + 24 * x118 * x131 + x118 * x254 + x118 * x307 + 8 * x119 * x210 - 80 * x12 * x79 - x121 * x13 * x404 - 24 * x121 * x132 + x121 * x557 - x121 * x565 - x121 * x572 + x121 * x81 + x122 * x65 - x123 * x684 * xi_xy - x123 * x716 * x730 - x124 * x136 * x589 - x124 * x229 * x696 - x125 * x38 + x125 * x58 - x126 * x558 - x126 * x612 - x126 * x98 + x127 * x588 + x127 * x605 - x127 * x98 + x128 * x130 - x128 * x173 - x128 * x404 * x82 + x128 * x409 * x548 + x128 * x440 - 108 * x128 * x505 + x128 * x509 * x535 + x128 * x543 + x128 * x563 - x128 * x594 - x128 * x597 - 64 * x128 * x650 + x129 * x196 * x85 + x129 * x302 - x13 * x473 - x13 * x49 + x130 * x131 + x130 * x250 - x131 * x173 + x131 * x181 + x131 * x543 + x131 * x563 - x131 * x594 - x131 * x597 - x132 * x133 - x132 * x143 + x133 * x179 + x133 * x180 - x133 * x696 - x134 * x137 - x134 * x139 - x134 * x253 - x135 * x518 * x643 - x135 * x530 * x669 - x136 * x16 * x549 * x645 - x136 * x165 * x266 + x136 * x182 - x136 * x258 - 36 * x136 * x610 + x136 * x634 * x78 + x137 * x478 + x138 * x182 - x138 * x242 - x138 * x258 + x138 * x265 + x138 * x515 * x651 - x138 * x568 - x138 * x578 * x82 + x139 * x478 + x14 * x196 * x202 + x14 * x42 * x89 + x14 * x534 + x14 - x140 * x410 - x140 * x448 - x140 * x451 - x140 * x47 + x140 * x646 * x648 + x140 * x83 + x140 * x88 - 192 * x141 * x225 * x412 * x625 + x142 * x557 - x142 * x565 + x142 * x630 * x646 - x143 * x571 - x144 * x47 + x144 * x83 + x144 * x88 + x145 * x147 + x145 * x527 + x148 * x149 - x148 * x82 + x149 * x165 * x243 + x149 * x255 + x149 * x312 + x150 * x61 + x151 * x152 + x151 * x53 + x151 * x55 + x152 * x403 - x153 * x154 - x153 * x155 - x154 * x447 - x154 * x482 * x600 + x154 * x569 - x155 * x447 + x155 * x569 + x155 * x570 + x156 * x158 + x156 * x178 - x157 * x501 * x94 + x157 * x665 * x734 + x158 * x215 * x734 - x158 * x555 - x16 * x232 * x589 - x16 * x250 * x565 + x16 * x408 * x43 + x16 * x566 * x627 + x16 * x665 * x666 - x161 * x285 - x161 * x287 - x161 * x32 + x162 * x163 + 28 * x162 * x195 - x163 * x531 + x165 * x197 * x28 + 64 * x165 * x267 + x165 * x55 + x166 * x243 - x167 * x487 - x167 * x606 - x167 * x64 - x168 * x487 - x168 * x606 - x168 * x64 + 24 * x169 * x290 - x17 * x202 * x67 + x170 * x593 - x171 * x34 - x171 * x37 + x172 * x404 - x172 * x67 - x173 * x250 + x174 * x223 - x174 * x454 * x575 - x175 * x67 - x176 * x29 - x176 * x53 - x177 * x340 + x178 * x31 * x436 - x178 * x499 * x7 - x178 * x555 - x178 * x625 * x644 - x179 * x188 - x18 * x340 - x180 * x247 + x181 * x250 - x183 * x184 - x183 * x185 - x183 * x186 - x183 * x213 - x183 * x214 + x184 * x490 + x185 * x490 + x186 * x490 + x186 * x76 + x187 * x68 * xi - x188 * x256 - x188 * x578 - x188 * x86 + x189 * x226 * x605 - x19 * x20 - x190 * x191 - x190 * x193 - x190 * x387 + x190 * x660 + x191 * x198 + x191 * x200 - x191 * x266 + 8 * x192 * x203 + x193 * x198 + x193 * x220 - x194 * x578 + x194 * x86 - x195 * x22 - x195 * x28 * x589 + x195 * x315 - x195 * x316 + x195 * x661 - x196 * x483 * x548 - x196 * x621 * x730 - x196 * x656 * x97 - x197 * x22 + x197 * x315 - x197 * x316 + x197 * x661 + x198 * x387 + x198 * x660 - x199 * x229 * x70 - x199 * x319 - x199 * x84 + x2 * x21 * x7 + x200 * x660 + x201 * x578 + x201 * x86 + x203 * x30 + x203 * x468 + x204 * x404 + x205 * x30 + x205 * x468 + x206 * x404 - x206 * x67 + x207 * x30 + x207 * x468 + x208 * x209 + x208 * x403 + x208 * x663 + x209 * x325 - x210 * x211 - x210 * x664 + x210 * x78 * xi - x211 * x212 - x212 * x664 + x213 * x490 + x213 * x76 + x214 * x490 + x214 * x76 - x215 * x303 * x750 - x215 * x45 * x554 - 60 * x215 * x493 * x731 + x215 * x494 * x666 - 10 * x215 * x507 + x215 * x517 * x729 + x216 * x327 - x216 * x346 - 48 * x217 * x68 - x218 * x219 + x218 * x224 - x219 * x424 - x219 * x444 - 136 * x221 * x79 + x222 * x230 - x222 * x246 - x222 * x251 + x222 * x272 - x222 * x659 - 56 * x222 * x92 + x223 * x5 + x224 * x424 + x224 * x444 + x225 * x270 * xi + x225 * x361 * x707 - x225 * x366 - x225 * x378 * x636 + x226 * x227 - x226 * x237 - x226 * x248 + x226 * x264 + x226 * x363 * x707 - x226 * x657 + x227 * x228 + x227 * x329 - x228 * x237 - x228 * x248 + x228 * x264 - x228 * x657 + x229 * x230 - x229 * x246 - x229 * x251 - x229 * x266 * x578 + x229 * x272 - x229 * x659 + 84 * x229 * x707 + x23 * x30 + x23 * x468 + x230 * x330 + x231 * x232 - x233 * x65 - x234 * x236 + x234 * x241 + 168 * x234 * x335 + x234 * x357 - x236 * x288 * xi - x237 * x329 + x239 * x437 + x239 * x452 * x79 - x239 * x58 - x24 * x404 - x24 * x67 - x240 * x241 - x240 * x357 - x241 * x519 - x242 * x304 - x243 * x269 - x243 * x561 + x243 * x607 + x243 * x633 + x243 * x654 + x244 * x245 + x244 * x334 + x245 * x593 - x246 * x330 + x247 * x56 * xi + x247 * x571 - x247 * x696 - x248 * x329 - x249 * x98 + x25 * x438 * x57 * xi - x251 * x330 + x253 * x478 - x255 * x82 - x257 * x404 - x258 * x335 + x259 * x261 - x26 * x384 + x26 * x425 * x426 + x26 * x429 * x433 + x26 * x456 + x26 * x78 * xi_xy + x262 * x263 + x262 * x271 - x263 * x555 - x263 * x618 + x265 * x304 - x266 * x660 - x267 * x447 + x267 * x569 + x267 * x570 - x268 * x606 - x268 * x64 + x27 * x405 + x270 * x95 - x271 * x555 - x271 * x618 - x274 * x623 * x669 - x276 * x278 + x276 * x347 * x709 + x277 * x284 + x277 * x385 + x277 * x386 + x279 * x321 + x279 * x343 + x28 * x556 * x670 - x280 * x282 + x280 * x336 + x280 * x341 + x280 * x342 - 168 * x281 * x333 + x282 * x710 - x283 * x284 - x283 * x385 - x283 * x386 + x284 * x707 + x284 * x711 - x285 * x286 - x285 * x295 + x285 * x322 + x285 * x326 + x285 * x365 + x285 * x374 - x285 * x379 + x285 * x712 + x285 * x715 - x286 * x287 + x287 * x322 + x287 * x365 + x287 * x374 - x287 * x379 + x287 * x712 + x287 * x715 + x288 * x292 + x288 * x368 - x289 * x349 - x289 * x351 - x289 * x99 - x29 * x56 - x291 * x292 - x291 * x368 - x293 * x294 - x293 * x486 + x293 * x72 - x294 * x296 - x296 * x486 - x296 * x611 * x96 + x296 * x72 + x296 * x74 - 92 * x296 * x90 + x298 * x299 + x298 * x747 + x299 * x318 + x300 * x588 + x300 * x605 - x300 * x98 - x301 * x558 - x301 * x98 - x302 * x695 - x303 * x332 - x303 * x531 * x726 - x304 * x305 + x304 * x338 + x304 * x344 + x304 * x355 - x304 * x376 + x304 * x393 - x304 * x396 + x304 * x713 + x304 * x736 + x304 * x737 - x304 * x739 - x304 * x744 - x304 * x748 - x304 * x749 - x305 * x306 + x306 * x338 + x306 * x344 + x306 * x355 - x306 * x376 + x306 * x393 - x306 * x396 + x306 * x713 + x306 * x736 + x306 * x737 - x306 * x739 - x306 * x744 - x306 * x748 - x306 * x749 - x308 * x446 * x96 - x31 * x34 - x31 * x37 + x310 * x83 + x310 * x88 - x311 * x578 + x311 * x86 - x312 * x82 - x313 * x404 - x314 * x404 + x318 * x747 + x319 * x491 - x319 * x662 - x32 * x457 * x575 + x32 * x475 * xi_x ^ 2 * (x16 * (S * x764 + x0 * x763 * x89 + x10 * x297 * x761 + x100 * (3 * x1 * x3 + x3 ^ 2 + x84 * x89) + x760 * x761) + x59 * x756 * (4 * x1 + 2 * x2 + x759 + 2)) + x32 * x550 - x320 * x578 - x320 * x86 + x321 * x709 - x323 * x487 - x323 * x64 - x323 * x66 - x324 * x487 - x324 * x64 - x324 * x66 + x325 * x663 + x327 * x573 - x328 * x366 - x329 * x657 + x33 * x404 * x65 + 24 * x33 * x58 - x330 * x659 - x332 * x352 * xi + x333 * x369 + x333 * x378 + x334 * x593 - x335 * x568 + x336 * x337 + x336 * x470 + x337 * x341 + x34 * x65 - x340 * x69 - x340 * x71 + x341 * x470 + x342 * x710 + x343 * x709 - x346 * x573 + x347 * x354 - x347 * x367 - x347 * x373 + x347 * x380 - x347 * x564 - x347 * x8 * x91 - x348 * x636 * x726 + x349 * x350 - x349 * x356 - x349 * x372 + x349 * x377 + x349 * x738 + x349 * x740 + x349 * x743 + x349 * x745 + x35 * x440 - x35 * x486 + x35 * x529 * xi_y + x35 * x72 + x35 * x74 - 60 * x35 * x90 + x350 * x351 - x351 * x356 - x351 * x372 + x351 * x377 + x351 * x738 + x351 * x740 + x351 * x743 + x351 * x745 - x352 * x359 + x352 * x741 + x353 * x354 - x353 * x367 - x353 * x373 + x353 * x380 - x359 * x401 + x36 * x81 - x360 * x361 - x362 * x364 + x362 * x390 + x362 * x395 + x364 * x742 + x369 * x371 + x369 * x402 + x37 * x65 + x370 * x601 * x615 + x371 * x378 - x375 * x578 + x375 * x86 + x378 * x402 - x378 * x650 * x95 - x38 * x39 + x381 * x7 * xi_x * (x16 * (-x4 * x5 * (-Fy * (u * (Gp * x4 - x775 * x776) + 2 * x765) + Gh * (10 * x1 + 6 * x2 + x759 + 6) - x755 * (Fyp * x4 + x5 * (Bh * x1 + Bh * x2 + Bh - 2 * Sh))) + x7 * (-Fxp * x10 * (x758 + 1) + x5 * (-Bp * x525 * x68 + Bt * x757 * x773 + 38 * Fx * x12 - Fx * x68 * x781 + Gt * x768 + St * (8 * x2 + x774 + 8) - u * x778 - x11 * x525 + x115 * x784 + 84 * x169 * x782 - x174 * x780 + x189 * x779 + x195 * x785 + 30 * x196 * x35 * x782 + x197 * x785 + 54 * x202 * x782 - x210 * x778 - x212 * x778 + x221 * x675 + x243 * x669 + x243 * x784 + 78 * x281 * x416 + x285 * x786 + x285 * x787 + x287 * x786 + x287 * x787 + x297 * x783 + x317 * x783 + x336 * x416 + x341 * x416 - x5 * x780 + x521 * x770 + x525 * x769 + x54 * x602 + 54 * x6 * x601 + x602 * x652 + x602 * x75 + x645 * x771 + x777 - x779 * x781 * xi))) + x25 * (Fyp * x10 + x5 * (-Gh * x10 * x768 - x146 * (-x11 * x4 + x113 * x341 + 19 * x12 + 7 * x221 + x243 * x341 + x285 * x772 + x287 * x772 + x317 * x767 + x766 + x769 + 16 * x770 + 18 * x771 + xi) - x3 * x455 + x4 * x7 * (u * (-Fx * Gp * x4 + Fx * x776 * (x774 + x775) - Fxp * x4 * x754 + x5 * x59 * (Bt * x1 + Bt * x2 + Bt + St)) + x407 * x773))) - x273 * xi_y * (x16 * x4 * (x123 * x5 * (2 * x1 + x3) - x765) + x25 * (x12 * x764 + x221 * x763 + x760 * x767 + x766))) + x382 * x429 + x383 * x449 + x383 * x450 + x385 * x707 + x385 * x711 + x386 * x707 + x388 * x389 + x388 * x391 + x389 * x394 - x39 * x437 + x390 * x392 + x390 * x742 + x390 * x751 + x391 * x394 + x392 * x395 + x395 * x742 + x395 * x751 - x397 * x533 * x669 * x90 - x398 * x399 + x398 * x462 + x398 * x466 - x399 * x400 + x40 * x41 - x40 * x48 + x400 * x462 + x400 * x466 + x401 * x741 + x403 * x53 + x403 * x55 - x406 * x482 - x409 * x418 + 2 * x409 * x428 - x409 * x608 - x41 * x427 * x5 - x41 * x580 * x681 - x410 * x85 + 6 * x411 * x412 - x411 * x457 + x412 * x528 * x673 - 72 * x412 * x642 * x686 + x413 * x415 + x413 * x494 * x495 + x414 * x506 * x586 - x414 * x547 + x414 * x672 + x417 * x477 + x418 * x474 - 48 * x42 * x519 + x42 * x685 * x686 - x420 * x421 - x420 * x678 + x421 * x429 * x7 - x421 * x442 + x422 * x423 + x422 * x540 + x423 * x443 - x424 * x544 - x424 * x581 - x424 * x582 - x425 * x494 * x499 + x426 * x439 * x586 + x426 * x57 * x671 - x428 * x600 * x601 + x429 * x653 - x43 * x45 - x430 * x454 + x432 * x434 + x432 * x552 + x432 * x583 - x433 * x536 + x434 * x441 + x435 * x504 - x436 * x73 - x439 * x602 * x716 + x439 * x628 * x78 - x44 * x493 * x515 + 10 * x44 * x57 + x441 * x552 + x441 * x583 - x442 * x678 + x443 * x540 - x444 * x544 - x444 * x581 - x444 * x582 - x445 * x446 - x445 * x452 - x448 * x46 - x451 * x46 - x453 * x54 - x453 * x75 - x458 * x750 * x90 - x459 * x460 + x459 * x469 + x459 * x488 + x459 * x708 + x459 * x752 + x459 * x753 - x46 * x47 - x460 * x461 + x461 * x469 + x461 * x488 + x461 * x708 + x461 * x752 + x461 * x753 + x463 * x464 + x463 * x465 + x464 * x467 + x465 * x467 + x471 * x472 + x471 * x709 + x472 * x489 - x473 * x50 + x474 * x595 + x477 * x522 + x477 * x591 + x477 * x676 + x477 * x694 + x478 * x479 + x479 * x491 + x480 * x481 + x481 * x598 + x481 * x599 + 24 * x483 * x484 - x487 * x62 + x489 * x709 + x49 * x50 + x491 * x70 * x93 + x492 * x53 + x492 * x55 - x494 * x625 * x626 - x497 * x500 - x500 * x617 - x501 * x623 * x624 + x502 * x503 - x505 * x727 + x507 * x78 + x508 * x509 - x508 * x510 + x509 * x534 * xi + x509 * x536 + x51 * x52 - x51 * x82 + x512 * x513 + x512 * x546 + x513 * x514 + x514 * x546 + x516 * x665 * x729 - x517 * x518 + x517 * x603 - x518 * x614 - x520 * x522 - x520 * x676 + x522 * x596 - x523 * x524 - x523 * x559 - x524 * x609 * x79 - x524 * x674 + x526 * x732 + x527 * x733 + x528 * x646 * x670 - x533 * x602 * x720 - 28 * x535 * x586 - x537 * x538 + x54 * x573 * x78 + x54 * x574 + x54 * x77 + x541 * x542 + x541 * x679 + x541 * x680 + x542 * x560 + x542 * x562 + x545 * x59 * x733 + x546 * x652 * x8 + x547 * x551 - x551 * x672 - x555 * x619 - x559 * x674 + x560 * x658 + x562 * x658 + x562 * x679 + x562 * x680 - x564 * x93 - x566 * x567 + x574 * x75 + x576 * x577 + x576 * x598 + x576 * x599 + x576 * x714 + x577 * x579 - x578 * x87 + x579 * x598 + x579 * x599 + x579 * x714 + 2 * x580 * x683 + x582 * x683 - x585 * x586 - x585 * x587 + x586 * x687 + x586 * x688 + x587 * x687 + x587 * x688 - x590 * x60 + x591 * x592 + x592 * x601 * x609 + x595 * x604 + x596 * x676 + x596 * x694 + x60 * x61 - x601 * x640 * x642 + x604 * x608 + x609 * x610 - 20 * x611 * x670 + x613 * x698 + x614 * x615 + x616 * x689 + x616 * x692 + x616 * x693 + 24 * x616 * x697 * x699 + x616 * x700 - x618 * x619 - x62 * x64 - x62 * x66 - x620 * x622 - x620 * x702 * x703 - x624 * x625 * x646 * x95 - x626 * x699 * x720 + x631 * x632 - 42 * x634 * x670 + x637 * x651 - x639 * x641 + x639 * x689 + x639 * x693 + x639 * x700 - x641 * x706 - x644 * x665 * x731 + x647 * x648 - x65 * x691 * x93 - x650 * x727 + x655 * x701 * x79 + x655 * x717 - x662 * x84 - x665 * x667 + x665 * x728 - x667 * x668 + x668 * x728 - x669 * x671 * x703 + x679 * x682 + x680 * x682 + x684 * x685 + x689 * x706 + x692 * x706 + x693 * x706 + x698 * x705 + x708 * x95 + x716 * x723 - x716 * x724 - x717 * x718 + x719 * x723 - x719 * x724 + x720 * x721 - x720 * x725 + x721 * x722 - x722 * x725 + x732 * x734 + x75 * x77 + x83 * x85 + x85 * x88 - x86 * x87 - x92 * x94

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

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

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
x37 = G * x18
x38 = 6 * x37
x39 = Fyh * x22
x40 = S * x20
x41 = 6 * G
x42 = x2 * x23
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
x55 = x25 * xi_y
x56 = 16 * Sh
x57 = Bh ^ 2
x58 = x17 * x24
x59 = 2 * x58
x60 = Fy ^ 2
x61 = Fyp ^ 2
x62 = Gh ^ 2
x63 = S ^ 2
x64 = S ^ 3
x65 = u ^ 7
x66 = S ^ 4
x67 = u ^ 10
x68 = 12 * B
x69 = Bh * x68
x70 = u ^ 8
x71 = x17 * x70
x72 = Fy * x71
x73 = x17 * xi_y
x74 = x65 * x73
x75 = u ^ 9
x76 = Bp * x73
x77 = Bp * Fyh
x78 = x20 * x24
x79 = S * x78
x80 = x20 * x49
x81 = Fy * Fyp
x82 = G * Gh
x83 = 12 * x82
x84 = Fy * Sh
x85 = 12 * x37
x86 = x70 * x85
x87 = x24 * x37
x88 = Fy * xi_y
x89 = S * x70
x90 = 28 * x89
x91 = 20 * x34
x92 = x24 * xi
x93 = S * x58
x94 = x58 * xi
x95 = x85 * x89
x96 = Fyh * xi
x97 = 12 * x87
x98 = 8 * x93
x99 = 24 * x2
x100 = x65 * x99
x101 = G * Gd
x102 = x100 * x101
x103 = x22 * xi
x104 = x101 * x99
x105 = x17 * x65
x106 = x45 * x82
x107 = S * x65
x108 = x85 * xi_yy
x109 = x37 * x65
x110 = 4 * x33
x111 = S * x75
x112 = x33 * xi_y
x113 = 32 * x112
x114 = x65 * xi
x115 = S * x114
x116 = Sp * x54
x117 = S * x49
x118 = Sp * x99
x119 = S * x17
x120 = x32 * xi
x121 = xi_y ^ 2
x122 = x111 * x20
x123 = x114 * x20
x124 = Bp * x60
x125 = xi ^ 2
x126 = x60 * x65
x127 = 42 * x119
x128 = 8 * x63 * x71
x129 = x121 * x35
x130 = x63 * x75
x131 = 2 * x19
x132 = x125 * x22
x133 = xi ^ 3
x134 = xi ^ 4
x135 = x24 * x63
x136 = x12 * x23
x137 = u ^ 11
x138 = x119 * x137
x139 = 24 * B
x140 = Bh * Fy * x139
x141 = x75 * xi
x142 = x139 * x141
x143 = Bh * xi_y
x144 = x119 * x139 * x143
x145 = Bp * x88
x146 = x40 * x70
x147 = x114 * x40
x148 = 24 * x82
x149 = Fy * x138
x150 = x141 * x36
x151 = x84 * x85
x152 = S * x137
x153 = x37 * x88
x154 = x111 * xi
x155 = 32 * x154
x156 = x119 * x88
x157 = x111 * x85
x158 = 72 * x101
x159 = x70 * xi
x160 = S * x159 * x2
x161 = x148 * xi_y
x162 = x71 * xi
x163 = S * x67
x164 = Sh * x45
x165 = S * x159
x166 = x36 * xi_y
x167 = B ^ 2
x168 = 36 * x167
x169 = x168 * x75
x170 = u ^ 12
x171 = x67 * xi
x172 = x171 * x40
x173 = x63 * x67
x174 = 2 * x17
x175 = x174 * x28
x176 = x17 * x31
x177 = G ^ 2
x178 = 36 * x177
x179 = x178 * x75
x180 = x137 * x63
x181 = x125 * x65
x182 = 12 * x181
x183 = 72 * x37
x184 = x163 * x183
x185 = x37 * x60
x186 = 48 * x159
x187 = x60 * xi
x188 = S * x71
x189 = 30 * x188
x190 = Fyh * x38
x191 = 36 * x63
x192 = x2 * x67
x193 = x191 * x192
x194 = u ^ 13
x195 = x194 * x64
x196 = Gd * x41
x197 = u ^ 16
x198 = x197 * x2 * x66
x199 = x125 * x2
x200 = x101 * x199
x201 = 36 * x24
x202 = x134 * x2 * x70
x203 = x121 * x183
x204 = x121 * xi
x205 = x125 * x24
x206 = u ^ 15
x207 = x47 * x51
x208 = Sp * x14
x209 = x114 * x51
x210 = x119 * xi
x211 = x170 * x210
x212 = x137 * xi
x213 = x163 * xi
x214 = x119 * x212
x215 = 18 * x167
x216 = x17 * x60
x217 = x216 * x67
x218 = x121 * x71
x219 = 2 * x57
x220 = x170 * x63
x221 = x17 * x220
x222 = x125 * x71
x223 = 18 * x177
x224 = 15 * x63
x225 = 2 * x62
x226 = u ^ 14
x227 = x226 * x63
x228 = x36 * x69
x229 = x125 * x67
x230 = x194 * x63
x231 = x69 * x73
x232 = x125 * x75
x233 = 72 * x167
x234 = x156 * x170
x235 = x166 * x171
x236 = x36 * x83
x237 = x125 * x70
x238 = 72 * x177
x239 = 84 * S
x240 = x133 * x163
x241 = x2 * x63
x242 = x212 * x241
x243 = x226 * x64 * xi
x244 = x106 * x17
x245 = x11 * x18
x246 = x10 * x245
x247 = x119 * x168
x248 = x194 * x60
x249 = x212 * x216
x250 = x121 * x137
x251 = x17 * x204
x252 = 2 * x222
x253 = 54 * x230
x254 = 30 * x232
x255 = x119 * x178
x256 = x37 * x63
x257 = x121 * x125
x258 = x156 * x194 * xi
x259 = exp(2 * x1)
x260 = Bb * x259
x261 = Gb * x259
x262 = Sb * x259
x263 = 2 * Gt
x264 = Bh * x2
x265 = x263 * x264
x266 = x18 * x42
x267 = Fx * Fy
x268 = x0 * x245
x269 = Fx * xi_y
x270 = 2 * Gp * x25
x271 = x136 * x18
x272 = x268 * xi
x273 = x18 * x2
x274 = Fxp * Fyp
x275 = 2 * x274
x276 = Fy * Gpt
x277 = x11 * x276
x278 = 8 * x266
x279 = Fy * Spt
x280 = x11 * x29
x281 = Gp * xi_xy
x282 = Gpt * xi_y
x283 = 8 * Sh
x284 = x0 * x18
x285 = Spt * xi_y
x286 = x166 * x191 * x206
x287 = x125 * x137 * x166
x288 = x187 * x226
x289 = x170 * x204
x290 = x170 * x191
x291 = 6 * B
x292 = x17 * x291
x293 = Fxt * x259
x294 = x22 * x293
x295 = x215 * x216
x296 = x197 * x63
x297 = x125 * x170
x298 = x17 * x215
x299 = x121 * x227
x300 = x257 * x67
x301 = Fx * x259
x302 = x301 * x32
x303 = Bt * x17
x304 = Bt * x259
x305 = Gt * x18
x306 = x305 * x44
x307 = Fxp * x301
x308 = St * x301
x309 = x216 * x223
x310 = x17 * x223
x311 = x18 * x261
x312 = St * x259
x313 = Gt * x2
x314 = x313 * xi_y
x315 = x292 * x314
x316 = x11 * x17
x317 = Bh * Gt
x318 = Gp * x267
x319 = Fx * x36
x320 = Gpp * x11
x321 = x18 * x267
x322 = x321 * x51
x323 = x267 * x273
x324 = Fx * x2
x325 = 2 * xi_y
x326 = Gp * x325
x327 = Fx * x320
x328 = x12 * x18
x329 = x269 * x328
x330 = Spp * x269
x331 = Fxh * Gp
x332 = x119 * x65
x333 = x11 * x332
x334 = x103 * x316
x335 = Fxh * x328
x336 = x276 * x50
x337 = Gpt * xi
x338 = x18 * x50
x339 = x107 * x338
x340 = x120 * x273
x341 = Fyt * Gp
x342 = Fyt * x328
x343 = x281 * x50
x344 = x282 * x50
x345 = x278 * xi
x346 = Bt ^ 2
x347 = x259 * x59
x348 = Fx ^ 2
x349 = Fxp ^ 2
x350 = Gt ^ 2
x351 = x17 ^ 2
x352 = Bd * x351
x353 = x291 * x352
x354 = Bt * x301
x355 = x68 * x71
x356 = x301 * x305
x357 = x304 * x305
x358 = St * x304
x359 = 12 * G
x360 = Gt * x301
x361 = x293 * xi
x362 = 4 * x305 * x312
x363 = x313 * x68
x364 = Gt * xi_y
x365 = x314 * x68
x366 = 12 * Bh * G * x324
x367 = x141 * x17
x368 = x11 * x67
x369 = Gp * x319
x370 = Gpp * x50
x371 = x267 * x370
x372 = Spp * x267 * x338
x373 = Gp * x11 * x269
x374 = x269 * x370
x375 = Gpp * x73
x376 = x188 * xi
x377 = x11 * x376
x378 = x114 * x119
x379 = B * x17
x380 = x259 * x348
x381 = x380 * x65
x382 = 18 * x381
x383 = x174 * x260
x384 = x122 * x259
x385 = x123 * x259
x386 = 2 * x311
x387 = B * x352
x388 = x100 * x387
x389 = x387 * x99
x390 = x174 * x241 * x67
x391 = x199 * x59
x392 = x368 * x63
x393 = x125 * x58
x394 = x281 * x316
x395 = x11 * x73
x396 = Gpt * x395
x397 = x139 * x354
x398 = x141 * x354
x399 = 24 * x37
x400 = x354 * x399
x401 = x301 * x303
x402 = 24 * G * x360
x403 = x308 * x85
x404 = x380 * x67
x405 = 36 * B
x406 = 72 * B
x407 = 48 * B
x408 = x380 * xi
x409 = x292 * x293
x410 = 20 * x180
x411 = 4 * x357
x412 = x259 * x346
x413 = x174 * x307
x414 = x37 * x380
x415 = x293 * x38
x416 = x259 * x350
x417 = x352 * x406
x418 = Gt * x291 * x36
x419 = x226 * x241
x420 = x199 * x67
x421 = Bh * x17 * x41
x422 = Fx * x326
x423 = S * x170 * xi
x424 = 2 * x221
x425 = x199 * x387
x426 = x401 * x68
x427 = x356 * x68
x428 = S * x406
x429 = x212 * x414
x430 = x17 * x359 * x360
x431 = x253 * x380
x432 = x194 * x380
x433 = x17 * x212 * x380
x434 = x226 * x408
x435 = x298 * x380
x436 = x310 * x380
x437 = 3 * G * u
x438 = x18 * x6
x439 = u * x6
x440 = 2 * x4
x441 = x1 * x6
x442 = 3 * x441
x443 = x6 ^ 2
x444 = 9 * x0
x445 = x177 * x24
x446 = 9 * x445
x447 = 6 * S * x5 * (3 * x445 + 2) + x0 * x63 * (x446 + 4) + x444 * (G * x4 + G) ^ 2
x448 = G * x10
x449 = 17 * x3
x450 = 9 * x177
x451 = 5 * x3 + 3 * x4 + x442 + 3
x452 = 9 * x4 + x449 + 7
x453 = 3 * x16
x454 = 36 * x111
x455 = 36 * x114
x456 = 36 * x213
ABCS[1] = x8 * x9
ABCS[2] = 8 * x10 * x8
ABCS[3] = u * x11 * x7
ABCS[4] = 12 * B * Fx * Gh * S * x137 * x17 * x2 + 12 * B * Fx * Gh * S * x17 * x170 * x2 * xi + 6 * B * Fx * Gh * x125 * x17 * x2 * x67 + 6 * B * Fx * Gh * x17 * x2 * x226 * x63 + 6 * B * Fx * Gh * x17 * x2 * x70 + 12 * B * Fx * Gh * x17 * x2 * x75 * xi + 144 * B * Fy * G * S * x170 * x18 * xi_y + 144 * B * Fy * G * S * x18 * x194 * xi * xi_y + 72 * B * Fy * G * x125 * x137 * x18 * xi_y + 72 * B * Fy * G * x18 * x206 * x63 * xi_y + 144 * B * Fy * G * x18 * x67 * xi * xi_y + 72 * B * Fy * G * x18 * x75 * xi_y + 24 * B * Fy * Gh * S * x137 * x18 + 24 * B * Fy * Gh * S * x170 * x18 * xi + 12 * B * Fy * Gh * x125 * x18 * x67 + 12 * B * Fy * Gh * x18 * x226 * x63 + 12 * B * Fy * Gh * x18 * x70 + 24 * B * Fy * Gh * x18 * x75 * xi + 12 * B * Fy * S * Sh * x137 * x17 + 156 * B * Fy * S * x17 * x67 * xi * xi_y + 144 * B * Fy * S * x17 * x75 * xi_y + 12 * B * Fy * Sh * x17 * x70 + 12 * B * Fy * Sh * x17 * x75 * xi + 54 * B * Fy * x125 * x17 * x70 * xi_y + 102 * B * Fy * x17 * x170 * x63 * xi_y + 42 * B * Fy * x17 * x24 * xi_y + 96 * B * Fy * x17 * x65 * xi * xi_y + 12 * B * Fyh * S * x17 * x70 + 12 * B * Fyh * S * x17 * x75 * xi + 6 * B * Fyh * x125 * x17 * x65 + 6 * B * Fyh * x137 * x17 * x63 + 6 * B * Fyh * x17 * x22 + 12 * B * Fyh * x17 * x24 * xi + 72 * B * G * S * x121 * x137 * x18 + 72 * B * G * S * x121 * x170 * x18 * xi + 72 * B * G * S * x18 * x194 * x60 + 72 * B * G * S * x18 * x226 * x60 * xi + 36 * B * G * x121 * x125 * x18 * x67 + 36 * B * G * x121 * x18 * x226 * x63 + 36 * B * G * x121 * x18 * x70 + 72 * B * G * x121 * x18 * x75 * xi + 36 * B * G * x125 * x170 * x18 * x60 + 72 * B * G * x137 * x18 * x60 * xi + 36 * B * G * x18 * x197 * x60 * x63 + 36 * B * G * x18 * x60 * x67 + 24 * B * Gh * S * x137 * x18 * xi * xi_y + 24 * B * Gh * S * x18 * x67 * xi_y + 12 * B * Gh * x125 * x18 * x75 * xi_y + 12 * B * Gh * x18 * x194 * x63 * xi_y + 12 * B * Gh * x18 * x65 * xi_y + 24 * B * Gh * x18 * x70 * xi * xi_y + 12 * B * S * Sh * x17 * x67 * xi_y + 72 * B * S * x121 * x17 * x70 + 72 * B * S * x121 * x17 * x75 * xi + 84 * B * S * x137 * x17 * x60 * xi + 72 * B * S * x17 * x60 * x67 + 12 * B * S * x17 * x65 * xi_yy + 12 * B * S * x17 * x70 * xi * xi_yy + 12 * B * Sh * x17 * x65 * xi_y + 12 * B * Sh * x17 * x70 * xi * xi_y + 24 * B * x121 * x125 * x17 * x65 + 48 * B * x121 * x137 * x17 * x63 + 24 * B * x121 * x17 * x22 + 48 * B * x121 * x17 * x24 * xi + 6 * B * x125 * x17 * x24 * xi_yy + 30 * B * x125 * x17 * x60 * x75 + 54 * B * x17 * x194 * x60 * x63 + 12 * B * x17 * x22 * xi * xi_yy + 6 * B * x17 * x23 * xi_yy + 18 * B * x17 * x60 * x65 + 48 * B * x17 * x60 * x70 * xi + 6 * B * x17 * x63 * x67 * xi_yy - B * x191 * x197 * x414 - 84 * B * x214 * x380 + 24 * Bd * Bp * S * x125 * x2 * x351 * x70 + 8 * Bd * Bp * S * x133 * x2 * x351 * x75 + 8 * Bd * Bp * S * x2 * x24 * x351 + 24 * Bd * Bp * S * x2 * x351 * x65 * xi + 2 * Bd * Bp * x0 * x2 * x351 + 12 * Bd * Bp * x125 * x137 * x2 * x351 * x63 + 12 * Bd * Bp * x125 * x2 * x22 * x351 + 8 * Bd * Bp * x133 * x2 * x24 * x351 + 2 * Bd * Bp * x134 * x2 * x351 * x65 + 8 * Bd * Bp * x170 * x2 * x351 * x64 + 8 * Bd * Bp * x194 * x2 * x351 * x64 * xi + 2 * Bd * Bp * x2 * x206 * x351 * x66 + 8 * Bd * Bp * x2 * x23 * x351 * xi + 24 * Bd * Bp * x2 * x351 * x63 * x67 * xi + 12 * Bd * Bp * x2 * x351 * x63 * x75 + 24 * Bh * Fy * G * S * x137 * x18 + 24 * Bh * Fy * G * S * x170 * x18 * xi + 12 * Bh * Fy * G * x125 * x18 * x67 + 12 * Bh * Fy * G * x18 * x226 * x63 + 12 * Bh * Fy * G * x18 * x70 + 24 * Bh * Fy * G * x18 * x75 * xi + 28 * Bh * Fy * S * x17 * x70 + 32 * Bh * Fy * S * x17 * x75 * xi + 12 * Bh * Fy * x125 * x17 * x65 + 20 * Bh * Fy * x137 * x17 * x63 + 8 * Bh * Fy * x17 * x22 + 20 * Bh * Fy * x17 * x24 * xi + 24 * Bh * G * S * x137 * x18 * xi * xi_y + 24 * Bh * G * S * x18 * x67 * xi_y + 12 * Bh * G * x125 * x18 * x75 * xi_y + 12 * Bh * G * x18 * x194 * x63 * xi_y + 12 * Bh * G * x18 * x65 * xi_y + 24 * Bh * G * x18 * x70 * xi * xi_y + 8 * Bh * Gh * S * x18 * x67 * xi + 8 * Bh * Gh * S * x18 * x75 + 4 * Bh * Gh * x125 * x18 * x70 + 4 * Bh * Gh * x170 * x18 * x63 + 4 * Bh * Gh * x18 * x24 + 8 * Bh * Gh * x18 * x65 * xi + 4 * Bh * S * Sh * x17 * x75 + 32 * Bh * S * x17 * x65 * xi_y + 32 * Bh * S * x17 * x70 * xi * xi_y + 4 * Bh * Sh * x17 * x24 + 4 * Bh * Sh * x17 * x65 * xi + 12 * Bh * x125 * x17 * x24 * xi_y - Bh * x142 * x36 + 24 * Bh * x17 * x22 * xi * xi_y + 12 * Bh * x17 * x23 * xi_y + 20 * Bh * x17 * x63 * x67 * xi_y + 4 * Bp * Fxt * S * x17 * x259 * x65 + 4 * Bp * Fxt * S * x17 * x259 * x70 * xi + 2 * Bp * Fxt * x125 * x17 * x24 * x259 + 4 * Bp * Fxt * x17 * x22 * x259 * xi + 2 * Bp * Fxt * x17 * x23 * x259 + 2 * Bp * Fxt * x17 * x259 * x63 * x67 + 4 * Bp * S * x17 * x259 * x348 * x67 * xi + 4 * Bp * S * x17 * x259 * x348 * x75 + 2 * Bp * x125 * x17 * x259 * x348 * x70 - Bp * x125 * x30 * x74 - Bp * x146 * x96 + 2 * Bp * x17 * x170 * x259 * x348 * x63 + 2 * Bp * x17 * x24 * x259 * x348 + 4 * Bp * x17 * x259 * x348 * x65 * xi - Bp * x20 * x39 * xi - 2 * Bp * x216 * x220 + 4 * Bs * S * x17 * x24 + 4 * Bs * S * x17 * x65 * xi + 2 * Bs * x0 * x17 + 2 * Bs * x125 * x17 * x22 + 4 * Bs * x17 * x23 * xi + 2 * Bs * x17 * x63 * x75 + 12 * Bt * Fy * G * S * x137 * x17 * x2 + 12 * Bt * Fy * G * S * x17 * x170 * x2 * xi + 6 * Bt * Fy * G * x125 * x17 * x2 * x67 + 6 * Bt * Fy * G * x17 * x2 * x226 * x63 + 6 * Bt * Fy * G * x17 * x2 * x70 + 12 * Bt * Fy * G * x17 * x2 * x75 * xi + 12 * Bt * G * S * x137 * x17 * x2 * xi * xi_y + 12 * Bt * G * S * x17 * x2 * x67 * xi_y + 6 * Bt * G * x125 * x17 * x2 * x75 * xi_y + 6 * Bt * G * x17 * x194 * x2 * x63 * xi_y + 6 * Bt * G * x17 * x2 * x65 * xi_y + 12 * Bt * G * x17 * x2 * x70 * xi * xi_y + 4 * Bt * Gh * S * x17 * x2 * x67 * xi + 4 * Bt * Gh * S * x17 * x2 * x75 + 2 * Bt * Gh * x125 * x17 * x2 * x70 + 2 * Bt * Gh * x17 * x170 * x2 * x63 + 2 * Bt * Gh * x17 * x2 * x24 + 4 * Bt * Gh * x17 * x2 * x65 * xi - Bt * St * x259 * x78 + 168 * Fx * Fy * G * S * x137 * x17 * x2 * xi + 144 * Fx * Fy * G * S * x17 * x2 * x67 + 60 * Fx * Fy * G * x125 * x17 * x2 * x75 + 108 * Fx * Fy * G * x17 * x194 * x2 * x63 + 36 * Fx * Fy * G * x17 * x2 * x65 + 96 * Fx * Fy * G * x17 * x2 * x70 * xi + 72 * Fx * Fy * S * x177 * x18 * x194 * x2 + 72 * Fx * Fy * S * x177 * x18 * x2 * x226 * xi + 84 * Fx * Fy * S * x18 * x2 * x65 + 60 * Fx * Fy * S * x18 * x2 * x70 * xi + 36 * Fx * Fy * x125 * x170 * x177 * x18 * x2 + 72 * Fx * Fy * x137 * x177 * x18 * x2 * xi + 36 * Fx * Fy * x177 * x18 * x197 * x2 * x63 + 36 * Fx * Fy * x177 * x18 * x2 * x67 + 30 * Fx * Fy * x18 * x2 * x63 * x67 + 4 * Fx * Fyp * S * x18 * x2 * x24 + 4 * Fx * Fyp * S * x18 * x2 * x65 * xi + 2 * Fx * Fyp * x0 * x18 * x2 + 2 * Fx * Fyp * x125 * x18 * x2 * x22 + 4 * Fx * Fyp * x18 * x2 * x23 * xi + 2 * Fx * Fyp * x18 * x2 * x63 * x75 + 24 * Fx * G * Gh * S * x137 * x18 * x2 + 24 * Fx * G * Gh * S * x170 * x18 * x2 * xi + 12 * Fx * G * Gh * x125 * x18 * x2 * x67 + 12 * Fx * G * Gh * x18 * x2 * x226 * x63 + 12 * Fx * G * Gh * x18 * x2 * x70 + 24 * Fx * G * Gh * x18 * x2 * x75 * xi + 12 * Fx * G * S * Sh * x137 * x17 * x2 + 156 * Fx * G * S * x17 * x2 * x67 * xi * xi_y + 144 * Fx * G * S * x17 * x2 * x75 * xi_y + 12 * Fx * G * Sh * x17 * x2 * x70 + 12 * Fx * G * Sh * x17 * x2 * x75 * xi + 54 * Fx * G * x125 * x17 * x2 * x70 * xi_y + 102 * Fx * G * x17 * x170 * x2 * x63 * xi_y + 42 * Fx * G * x17 * x2 * x24 * xi_y + 96 * Fx * G * x17 * x2 * x65 * xi * xi_y + 28 * Fx * Gh * S * x17 * x2 * x70 + 32 * Fx * Gh * S * x17 * x2 * x75 * xi + 12 * Fx * Gh * x125 * x17 * x2 * x65 + 20 * Fx * Gh * x137 * x17 * x2 * x63 + 8 * Fx * Gh * x17 * x2 * x22 + 20 * Fx * Gh * x17 * x2 * x24 * xi - Fx * Gp * x154 * x395 + 72 * Fx * S * x170 * x177 * x18 * x2 * xi_y + 72 * Fx * S * x177 * x18 * x194 * x2 * xi * xi_y + 68 * Fx * S * x18 * x2 * x24 * xi_y + 56 * Fx * S * x18 * x2 * x65 * xi * xi_y + 24 * Fx * Sh * x18 * x2 * x22 + 16 * Fx * Sh * x18 * x2 * x24 * xi + 36 * Fx * x125 * x137 * x177 * x18 * x2 * xi_y + 36 * Fx * x177 * x18 * x2 * x206 * x63 * xi_y + 72 * Fx * x177 * x18 * x2 * x67 * xi * xi_y + 36 * Fx * x177 * x18 * x2 * x75 * xi_y - Fx * x18 * x199 * x24 * x30 + 24 * Fx * x18 * x2 * x63 * x75 * xi_y - Fx * x264 * x41 * x71 - Fx * x375 * x392 - Fx * x419 * x421 - Fx * x420 * x421 + 12 * Fxh * G * S * x17 * x2 * x70 + 12 * Fxh * G * S * x17 * x2 * x75 * xi + 6 * Fxh * G * x125 * x17 * x2 * x65 + 6 * Fxh * G * x137 * x17 * x2 * x63 + 6 * Fxh * G * x17 * x2 * x22 + 12 * Fxh * G * x17 * x2 * x24 * xi + 4 * Fxh * S * x18 * x2 * x22 + 8 * Fxh * S * x18 * x2 * x24 * xi + 8 * Fxh * x18 * x2 * x63 * x70 - Fxh * x2 * x270 - Fxh * x246 - Fxh * x271 - Fxh * x272 + 4 * Fxp * Fy * S * x18 * x2 * x24 + 4 * Fxp * Fy * S * x18 * x2 * x65 * xi + 2 * Fxp * Fy * x0 * x18 * x2 + 2 * Fxp * Fy * x125 * x18 * x2 * x22 + 4 * Fxp * Fy * x18 * x2 * x23 * xi + 2 * Fxp * Fy * x18 * x2 * x63 * x75 + 4 * Fxt * Gp * S * x18 * x259 * x65 + 4 * Fxt * Gp * S * x18 * x259 * x70 * xi + 2 * Fxt * Gp * x125 * x18 * x24 * x259 + 4 * Fxt * Gp * x18 * x22 * x259 * xi + 2 * Fxt * Gp * x18 * x23 * x259 + 2 * Fxt * Gp * x18 * x259 * x63 * x67 + 4 * Fxt * S * Sp * x17 * x259 * x65 + 4 * Fxt * Sp * x17 * x22 * x259 * xi + 4 * Fxt * Sp * x17 * x23 * x259 + 4 * Fxt * x0 * x17 * x259 * xi + 4 * Fxt * x10 * x17 * x259 + 24 * Fy * G * Gt * S * x137 * x18 * x2 + 24 * Fy * G * Gt * S * x170 * x18 * x2 * xi + 12 * Fy * G * Gt * x125 * x18 * x2 * x67 + 12 * Fy * G * Gt * x18 * x2 * x226 * x63 + 12 * Fy * G * Gt * x18 * x2 * x70 + 24 * Fy * G * Gt * x18 * x2 * x75 * xi + 12 * Fy * G * S * St * x137 * x17 * x2 + 12 * Fy * G * St * x17 * x2 * x70 + 12 * Fy * G * St * x17 * x2 * x75 * xi + 4 * Fy * Gp * S * x18 * x70 * xi_y + 4 * Fy * Gp * S * x18 * x75 * xi * xi_y + 2 * Fy * Gp * x125 * x18 * x65 * xi_y + 2 * Fy * Gp * x137 * x18 * x63 * xi_y + 2 * Fy * Gp * x18 * x22 * xi_y + 4 * Fy * Gp * x18 * x24 * xi * xi_y + 28 * Fy * Gt * S * x17 * x2 * x70 + 32 * Fy * Gt * S * x17 * x2 * x75 * xi + 12 * Fy * Gt * x125 * x17 * x2 * x65 + 20 * Fy * Gt * x137 * x17 * x2 * x63 + 8 * Fy * Gt * x17 * x2 * x22 + 20 * Fy * Gt * x17 * x2 * x24 * xi + 4 * Fy * S * Sp * x17 * x70 * xi_y + 4 * Fy * Sp * x17 * x22 * xi_y + 4 * Fy * Sp * x17 * x24 * xi * xi_y + 24 * Fy * St * x18 * x2 * x22 + 16 * Fy * St * x18 * x2 * x24 * xi + 4 * Fy * x0 * x17 * xi_y - Fy * x148 * x211 + 4 * Fy * x17 * x23 * xi * xi_y - Fy * x211 * x313 * x68 - Fy * x56 * x94 + 4 * Fyh * Gp * S * x18 * x65 + 4 * Fyh * Gp * S * x18 * x70 * xi + 2 * Fyh * Gp * x125 * x18 * x24 + 4 * Fyh * Gp * x18 * x22 * xi + 2 * Fyh * Gp * x18 * x23 + 2 * Fyh * Gp * x18 * x63 * x67 + 4 * Fyh * S * Sp * x17 * x65 + 4 * Fyh * Sp * x17 * x22 * xi + 4 * Fyh * Sp * x17 * x23 + 4 * Fyh * x0 * x17 * xi + 4 * Fyh * x10 * x17 - Fyh * x128 - Fyh * x95 + 12 * Fyt * G * S * x17 * x2 * x70 + 12 * Fyt * G * S * x17 * x2 * x75 * xi + 6 * Fyt * G * x125 * x17 * x2 * x65 + 6 * Fyt * G * x137 * x17 * x2 * x63 + 6 * Fyt * G * x17 * x2 * x22 + 12 * Fyt * G * x17 * x2 * x24 * xi + 4 * Fyt * S * x18 * x2 * x22 + 8 * Fyt * S * x18 * x2 * x24 * xi + 8 * Fyt * x18 * x2 * x63 * x70 - Fyt * x2 * x270 - Fyt * x246 - Fyt * x271 - Fyt * x272 + 24 * G * Gt * S * x137 * x18 * x2 * xi * xi_y + 24 * G * Gt * S * x18 * x2 * x67 * xi_y + 12 * G * Gt * x125 * x18 * x2 * x75 * xi_y + 12 * G * Gt * x18 * x194 * x2 * x63 * xi_y + 12 * G * Gt * x18 * x2 * x65 * xi_y + 24 * G * Gt * x18 * x2 * x70 * xi * xi_y + 12 * G * S * St * x17 * x2 * x67 * xi_y + 24 * G * S * x17 * x2 * x65 * xi_xy + 24 * G * S * x17 * x2 * x70 * xi * xi_xy + 12 * G * St * x17 * x2 * x65 * xi_y + 12 * G * St * x17 * x2 * x70 * xi * xi_y + 12 * G * x125 * x17 * x2 * x24 * xi_xy + 24 * G * x17 * x2 * x22 * xi * xi_xy + 12 * G * x17 * x2 * x23 * xi_xy + 12 * G * x17 * x2 * x63 * x67 * xi_xy + 8 * Gc * S * x17 * x2 * x24 + 8 * Gc * S * x17 * x2 * x65 * xi + 4 * Gc * x0 * x17 * x2 + 4 * Gc * x125 * x17 * x2 * x22 + 8 * Gc * x17 * x2 * x23 * xi + 4 * Gc * x17 * x2 * x63 * x75 + 24 * Gd * Gp * S * x125 * x2 * x70 + 8 * Gd * Gp * S * x133 * x2 * x75 + 8 * Gd * Gp * S * x2 * x24 + 24 * Gd * Gp * S * x2 * x65 * xi + 2 * Gd * Gp * x0 * x2 + 12 * Gd * Gp * x125 * x137 * x2 * x63 + 12 * Gd * Gp * x125 * x2 * x22 + 8 * Gd * Gp * x133 * x2 * x24 + 2 * Gd * Gp * x134 * x2 * x65 + 8 * Gd * Gp * x170 * x2 * x64 + 8 * Gd * Gp * x194 * x2 * x64 * xi + 2 * Gd * Gp * x2 * x206 * x66 + 8 * Gd * Gp * x2 * x23 * xi + 24 * Gd * Gp * x2 * x63 * x67 * xi + 12 * Gd * Gp * x2 * x63 * x75 - Gd * x41 * x42 + 8 * Gh * Gt * S * x18 * x2 * x67 * xi + 8 * Gh * Gt * S * x18 * x2 * x75 + 4 * Gh * Gt * x125 * x18 * x2 * x70 + 4 * Gh * Gt * x170 * x18 * x2 * x63 + 4 * Gh * Gt * x18 * x2 * x24 + 8 * Gh * Gt * x18 * x2 * x65 * xi + 4 * Gh * S * St * x17 * x2 * x75 + 4 * Gh * St * x17 * x2 * x24 + 4 * Gh * St * x17 * x2 * x65 * xi + 4 * Gp * S * x18 * x24 * xi_yy + 4 * Gp * S * x18 * x259 * x348 * x67 * xi + 4 * Gp * S * x18 * x259 * x348 * x75 + 4 * Gp * S * x18 * x60 * x67 * xi + 4 * Gp * S * x18 * x60 * x75 + 4 * Gp * S * x18 * x65 * xi * xi_yy + 2 * Gp * x0 * x18 * xi_yy + 2 * Gp * x125 * x18 * x22 * xi_yy + 2 * Gp * x125 * x18 * x259 * x348 * x70 + 2 * Gp * x125 * x18 * x60 * x70 + 2 * Gp * x170 * x18 * x259 * x348 * x63 + 2 * Gp * x170 * x18 * x60 * x63 + 4 * Gp * x18 * x23 * xi * xi_yy + 2 * Gp * x18 * x24 * x259 * x348 + 2 * Gp * x18 * x24 * x60 + 4 * Gp * x18 * x259 * x348 * x65 * xi + 4 * Gp * x18 * x60 * x65 * xi + 2 * Gp * x18 * x63 * x75 * xi_yy - Gpt * x36 * x392 + 4 * Gt * S * Sh * x17 * x2 * x75 + 32 * Gt * S * x17 * x2 * x65 * xi_y + 32 * Gt * S * x17 * x2 * x70 * xi * xi_y + 4 * Gt * Sh * x17 * x2 * x24 + 4 * Gt * Sh * x17 * x2 * x65 * xi + 12 * Gt * x125 * x17 * x2 * x24 * xi_y + 24 * Gt * x17 * x2 * x22 * xi * xi_y + 12 * Gt * x17 * x2 * x23 * xi_y + 20 * Gt * x17 * x2 * x63 * x67 * xi_y + 8 * S * Sc * x18 * x2 * x24 + 16 * S * Sd * x125 * x2 * x24 + 16 * S * Sd * x2 * x22 * xi + 4 * S * Sp * x17 * x24 * xi_yy + 4 * S * Sp * x17 * x259 * x348 * x75 + 4 * S * Sp * x17 * x60 * x75 + 16 * S * u * x2 + 72 * S * x0 * x125 * x2 + 56 * S * x10 * x2 * xi - S * x102 - S * x118 * x132 + 40 * S * x133 * x2 * x23 + 8 * S * x134 * x2 * x22 + 2 * S * x17 * x22 * x259 * x349 + 2 * S * x17 * x22 * x61 + 2 * S * x17 * x24 * x259 * x349 * xi + 2 * S * x17 * x24 * x61 * xi + 16 * S * x18 * x2 * x22 * xi * xi_xy + 16 * S * x18 * x2 * x23 * xi_xy - S * x22 * x245 * x274 - S * x25 * x52 - S * x388 + 8 * Sc * x0 * x18 * x2 + 8 * Sc * x18 * x2 * x23 * xi + 24 * Sd * x2 * x63 * x65 + 32 * Sd * x2 * x63 * x70 * xi + 16 * Sd * x2 * x64 * x67 + 4 * Sh ^ 2 * x17 * x24 - Sh * x109 * x45 - Sh * x110 * x111 - Sh * x110 * x114 - Sh * x33 * x44 - Sh * x35 * x36 + 4 * Sp * x0 * x17 * xi_yy + 4 * Sp * x17 * x23 * xi * xi_yy + 4 * Sp * x17 * x24 * x259 * x348 + 4 * Sp * x17 * x24 * x60 + 4 * Sp * x17 * x259 * x348 * x65 * xi + 4 * Sp * x17 * x60 * x65 * xi - Sp * x345 * xi_xy - Sp * x53 - Spp * x32 * x323 + 4 * St ^ 2 * x17 * x24 * x259 + 16 * St * x18 * x2 * x22 * xi * xi_y + 16 * St * x18 * x2 * x23 * xi_y - St * x24 * x273 * x283 - u * x14 - x0 * x133 * x15 + 2 * x0 * x17 * x259 * x349 * xi + 2 * x0 * x17 * x61 * xi - x0 * x20 * x262 - x0 * x21 + 2 * x0 * x259 * x6 * xi_xx * (x17 * (-3 * B * x439 + Bp * x6 - 4 * S * u + 2 * Sp) + x438 * (Gp - x437)) - x10 * x118 * x125 + x10 * x17 * x259 * x349 + x10 * x17 * x61 - x10 * x273 * x275 - x10 * x54 * xi - 72 * x101 * x111 * x199 - x101 * x193 - x102 * x133 - x103 * x104 - x103 * x108 - 12 * x103 * x323 - x103 * x335 - x103 * x342 - x103 * x389 - x103 * x56 * x73 - x104 * x195 - x104 * x240 - x104 * x243 - x105 * x106 - x105 * x199 * x422 - x107 * x108 - x107 * x113 - x107 * x335 - x107 * x342 - x108 * x165 - 24 * x109 * x257 - x11 * x220 * x369 - x11 * x222 * x318 - x11 * x318 * x58 - 144 * x111 * x153 - x111 * x17 * x361 * x68 - x111 * x199 * x352 * x406 - x111 * x316 * x317 - x111 * x322 - 8 * x111 * x357 - x111 * x362 - x111 * x369 * x50 - 20 * x112 * x173 - x112 * x35 * xi - x113 * x165 - 96 * x114 * x153 - 56 * x114 * x156 - x114 * x316 * x317 - 8 * x114 * x357 - x114 * x362 - x114 * x369 * x50 - x115 * x116 - 4 * x115 * x311 - x115 * x48 - x116 * x47 - x116 * x49 - x117 * x118 - x119 * x120 * xi_yy - x119 * x129 - x119 * x161 * x67 - x119 * x171 * x318 * x50 - x119 * x192 * x364 * x68 - x119 * x404 * x406 - x12 * x135 - x12 * x237 * x63 - x12 - x120 * x324 * x375 - x121 * x128 - x122 * x124 - x122 * x145 * xi - x122 * x358 - x122 * x57 - x122 * x62 - x123 * x124 - x123 * x358 - x123 * x57 - x123 * x62 - x124 * x172 - x124 * x252 - x124 * x59 + x125 * x17 * x23 * x259 * x349 + x125 * x17 * x23 * x61 + x125 * x17 * x24 * x259 * x348 + x125 * x17 * x24 * x60 + 68 * x125 * x2 * x24 * x63 + 8 * x125 * x2 * x64 * x75 - x125 * x266 * x275 - x125 * x53 - x125 * x59 * x77 - x126 * x127 - 18 * x126 * x37 - x127 * x381 - x128 * x293 - x129 * x37 - x130 * x131 - 24 * x130 * x166 - x130 * x175 - x130 * x176 - x130 * x208 - x130 * x383 - x130 * x386 - x130 * x394 - x130 * x396 - x130 * x413 - x131 * x132 - x132 * x175 - x132 * x176 - x132 * x208 - x132 * x383 - x132 * x386 - x132 * x394 - x132 * x396 - x132 * x413 + 16 * x133 * x2 * x63 * x65 - x133 * x207 - x133 * x388 - x134 * x136 - x137 * x17 * x241 * x422 - x138 * x140 - x138 * x308 * x68 - x138 * x366 - x138 * x397 - x138 * x402 - x139 * x143 * x71 * xi - x139 * x152 * x356 - x139 * x356 * x423 - x140 * x211 - x141 * x151 - x141 * x17 * x308 * x68 - x141 * x403 - x142 * x356 - x144 * x212 - x144 * x67 - x145 * x146 - x145 * x78 * xi - x147 * x260 - x147 * x28 - x147 * x307 - x147 * x81 - x148 * x149 - x148 * x150 - x149 * x363 - x15 * x4 - x150 * x363 - x151 * x152 - x152 * x400 - x152 * x403 - 156 * x153 * x213 - 102 * x153 * x220 - 54 * x153 * x237 - x154 * x203 - x154 * x319 * x370 - x155 * x34 - x155 * x356 - x155 * x401 - x157 * x361 - x157 * x96 - x158 * x160 - x158 * x242 - x159 * x164 * x37 - x160 * x417 - x161 * x162 - x161 * x214 - x162 * x365 - x163 * x164 * x37 - x166 * x169 - x166 * x179 - x167 * x286 - x168 * x249 - x168 * x287 - x168 * x433 - x169 * x251 - x17 * x173 * x26 + 6 * x17 * x22 * x259 * x348 * xi - x17 * x22 * x324 * x326 + 6 * x17 * x22 * x60 * xi - x17 * x224 * x404 + 9 * x17 * x23 * x259 * x348 + 9 * x17 * x23 * x60 + x17 * x259 * x349 * x63 * x70 - x17 * x308 * x35 - x17 * x52 * x63 * x65 + x17 * x61 * x63 * x70 - x172 * x412 - x172 * x416 - x172 * x57 - x172 * x62 - x173 * x43 - x177 * x286 - x178 * x249 - x178 * x287 - x178 * x433 - x179 * x251 + 16 * x18 * x2 * x63 * x65 * xi_xy - x18 * x207 * xi_xy - x18 * x241 * x275 * x70 - x18 * x261 * x9 - x180 * x190 - x180 * x30 * x76 - x180 * x327 * x36 - x180 * x409 - x180 * x415 - x180 * x91 - x181 * x190 - x181 * x327 * x36 - x181 * x409 - x181 * x415 - x182 * x34 - x182 * x356 - x182 * x401 - x184 * x380 - x184 * x60 - x185 * x186 - x185 * x212 * x239 - x185 * x253 - x185 * x254 - x186 * x414 - x187 * x189 - x188 * x293 * x68 - 28 * x188 * x354 - x188 * x371 - x188 * x373 - x189 * x408 - x19 * x9 - x193 * x387 - x194 * x241 * x292 * x364 - x194 * x414 * x428 - x195 * x389 - x196 * x198 - x196 * x202 - x198 * x353 + 88 * x2 * x22 * x63 * xi + 36 * x2 * x23 * x63 - x2 * x32 * x337 * x36 + 24 * x2 * x64 * x65 + 32 * x2 * x64 * x70 * xi + 4 * x2 * x66 * x67 - x2 * x9 * xi_x * (x17 * (x2 * (Fx * (78 * B * x115 + 51 * B * x130 + 27 * B * x132 - Bp * x10 * x443 - 2 * Sp * x10 * x6 + 21 * x1 + 28 * x117 + 12 * x135 + x167 * x454 + x167 * x455 + x167 * x456 + x177 * x454 + x177 * x455 + x177 * x456 + x215 * x220 + x215 * x237 + x215 * x24 + x220 * x223 + x223 * x237 + 34 * x3 + x406 * x47 + x407 * x49 - x440 + 18 * x445 - 2) + 2 * u * (Bt * x451 * x6 + Gt * x443 * x453 + St * (4 * x4 + x442 + 4))) - x439 * (-Fy * (2 * Gpp * x6 + u * (Gp * x6 - x437 * x452)) + Gh * (10 * x3 + 6 * x4 + x442 + 6) - x453 * (Bh * x3 + Bh * x4 + Bh - 2 * Sh))) + x18 * (u * (-6 * Gh * x16 * x443 + x2 * x6 * (u * (-Fx * Gp * x6 + Fx * x437 * (12 * x441 + x452) + 6 * x448 * (Bt * x3 + Bt * x4 + Bt + St)) + x263 * x451) - x283 * x5) + x30 * (Sp * x10 * x6 + 2 * Spp * u * x6 + u * xi - x111 * x223 - x114 * x223 - 14 * x117 - 6 * x135 - x213 * x223 - x220 * x450 - x237 * x450 - x446 - x449 + 1)) + x325 * (x17 * x6 * (Gpp * x6 - 12 * x448 * (2 * x3 + x5)) + x18 * (2 * Spp * x6 - x10 * x447))) - x200 * x201 - x200 * x290 - x201 * x425 - x202 * x353 - x203 * x89 - 48 * x204 * x87 - 24 * x204 * x93 - x205 * x43 - x205 * x46 - x209 * x321 - x209 * x63 - x21 * x47 - x21 * x49 - x210 * x317 * x368 - x211 * x366 - x211 * x397 - x211 * x402 - 8 * x213 * x357 - x214 * x365 - x215 * x217 - x215 * x218 - x217 * x223 - x217 * x224 - x218 * x223 - x219 * x221 - x219 * x222 - 2 * x22 * x259 * xi_x ^ 2 * (x17 * (x167 * x443 * x444 + x447 + x68 * (2 * x135 + 3 * x3 * x5 + x5 ^ 2)) + x41 * x438 * (4 * x3 + x440 + x442 + 2)) - x22 * x30 * x76 - x22 * x319 * x320 - x22 * x329 - x220 * x411 - x221 * x225 - x221 * x265 - x222 * x225 - x222 * x265 - 12 * x226 * x256 * x354 - x226 * x37 * x408 * x428 - x227 * x228 - x227 * x236 - x227 * x426 - x227 * x427 - x227 * x430 - x228 * x229 - x229 * x236 - x229 * x354 * x85 - x229 * x426 - x229 * x427 - x229 * x430 - x23 * x43 - x23 * x46 - x230 * x231 - x230 * x244 - x231 * x232 - x232 * x244 - x232 * x315 - x233 * x234 - x233 * x235 - x233 * x258 - x234 * x238 - x235 * x238 - x237 * x411 - x238 * x258 - x239 * x429 - x24 * x322 - x240 * x389 - x242 * x417 - x243 * x389 - x245 * x269 * x49 - x245 * x274 * x47 * xi - x247 * x248 - x247 * x250 - x247 * x288 - x247 * x289 - x247 * x432 - x247 * x434 - x248 * x255 - x25 * x26 - x25 * x277 - x25 * x343 * xi - x250 * x255 - 48 * x250 * x256 - x252 * x412 - x252 * x416 - x254 * x379 * x380 - x254 * x414 - x255 * x288 - x255 * x289 - x255 * x432 - x255 * x434 - x260 * x27 - x260 * x79 - x260 * x80 - x262 * x79 - x262 * x80 - x265 * x58 - 18 * x266 * x267 - x268 * x269 - x269 * x320 * x393 - x27 * x28 - x27 * x307 - x272 * x274 - x277 * x393 - x278 * x279 - x278 * x330 - x279 * x339 - x279 * x340 - x28 * x79 - x28 * x80 - x280 * x281 - x280 * x282 - x284 * x285 * x50 - x284 * x51 * xi_xy - x285 * x338 * x47 - x285 * x345 - x29 * x31 - x290 * x425 - x291 * x313 * x72 - x292 * x294 - x293 * x68 * x94 - x293 * x95 - x294 * x38 - x294 * x40 - x295 * x296 - x295 * x297 - x296 * x309 - x296 * x435 - x296 * x436 - x297 * x309 - x297 * x405 * x414 - x297 * x435 - x297 * x436 - x298 * x299 - x298 * x300 - x298 * x404 - x299 * x310 - x3 * x51 - x300 * x310 - x302 * x303 - x302 * x305 - x304 * x306 - x306 * x312 - x307 * x79 - x307 * x80 - x308 * x355 - x308 * x86 - 16 * x308 * x94 - x310 * x404 - 4 * x311 * x47 - 4 * x311 * x49 - x315 * x65 - x32 * x34 - x327 * x55 - x329 * x89 - x329 * x92 - x330 * x339 - x330 * x340 - x331 * x333 - x331 * x334 - x331 * x377 - x331 * x390 - x331 * x391 - x332 * x336 - x332 * x374 - x333 * x341 - x334 * x341 - x336 * x376 - x337 * x50 * x55 - x34 * x90 - x341 * x377 - x341 * x390 - x341 * x391 - x343 * x378 - x343 * x93 - x344 * x378 - x344 * x93 - x346 * x347 - x346 * x384 - x346 * x385 - x347 * x350 - x350 * x384 - x350 * x385 - x353 * x42 - x354 * x355 - x354 * x86 - 20 * x354 * x94 - x356 * x410 - x356 * x68 * x70 - x356 * x90 - 20 * x356 * x92 - x359 * x360 * x71 - x361 * x97 - x361 * x98 - x366 * x367 - x367 * x402 - x37 * x382 - x37 * x404 * x405 - x37 * x431 - x371 * x94 - x372 * x89 - x372 * x92 - x373 * x94 - x374 * x376 - x379 * x382 - 24 * x379 * x398 - x379 * x431 - x38 * x39 - x39 * x40 - x398 * x399 - x40 * x65 * x77 - x400 * x423 - x401 * x410 - x406 * x429 - x407 * x408 * x71 - x412 * x424 - x416 * x424 - x418 * x419 - x418 * x420 - x47 * x48 - x48 * x49 - x55 * x56 - x57 * x59 - x59 * x62 - x69 * x72 - x69 * x74 - x72 * x83 - x79 * x81 - x80 * x81 - x84 * x86 - 42 * x87 * x88 - 68 * x88 * x93 - x91 * x92 - x96 * x97 - x96 * x98

    nothing
end
