
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

x0 = u ^ 6
x1 = u ^ 5
x2 = u ^ 4
x3 = 3 * u
x4 = (-B * x3 + Bp) ^ 2 * cosh(G * u ^ 3) ^ 2
ABCS[1] = 4 * x0
ABCS[2] = 24 * x1
ABCS[3] = x2 * (9 * G ^ 2 * x0 - 6 * G * Gp * x1 + Gp ^ 2 * x2 + x2 * x4 + 24)
ABCS[4] = x1 * (x4 + (-G * x3 + Gp) ^ 2) * (S0 * (u * xi + 1 ) + S0_t)

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
    
    
x0 = S0_t * u
x1 = u * xi
x2 = S0 * x1
x3 = u ^ 3
x4 = S * x3
x5 = S0 + x2 + x4
x6 = x0 + x5
x7 = x6 ^ 2
x8 = 2 * u
x9 = x7 * x8
x10 = -Bp
x11 = 3 * u
x12 = B * x11
x13 = x10 + x12
x14 = -x13
x15 = u ^ 2
x16 = G * x3
x17 = 2 * x16
x18 = cosh(x17)
x19 = x15 * x18
x20 = Bp * x15
x21 = B * x3
x22 = 3 * x21
x23 = x20 - x22
x24 = 2 * Gp
x25 = 6 * u
x26 = G * x25
x27 = sinh(x17)
x28 = x14 * x27
x29 = x24 - x26 - x28
x30 = exp(x21)
x31 = x15 * x7
x32 = x30 * x31
x33 = 2 * S0
x34 = Spp * u
x35 = Sp * x15
x36 = Gp ^ 2
x37 = u ^ 7
x38 = S * x37
x39 = x36 * x38
x40 = u ^ 4
x41 = x36 * x40
x42 = u ^ 8
x43 = S * x42
x44 = 6 * Gp
x45 = G * x44
x46 = u ^ 5
x47 = G * x46
x48 = x44 * x47
x49 = G ^ 2
x50 = u ^ 9
x51 = 9 * S
x52 = x50 * x51
x53 = u ^ 6
x54 = 9 * x53
x55 = x49 * x54
x56 = S0 * xi
x57 = x46 * x56
x58 = x53 * x56
x59 = x37 * x56
x60 = -Gp
x61 = G * x11
x62 = x60 + x61
x63 = -x62
x64 = x28 * x40
x65 = x63 * x64
x66 = x5 * x65
x67 = Bpp * x5
x68 = Bp ^ 2
x69 = x15 * x68
x70 = x5 * x69
x71 = S * u
x72 = 9 * x71
x73 = S * x40
x74 = 6 * B
x75 = x73 * x74
x76 = 7 * x1
x77 = x1 + 1
x78 = 6 * x21
x79 = x77 * x78
x80 = x76 - x79
x81 = 5 * x1
x82 = -3 * B * x3 * x77 + x81
x83 = cosh(x16) ^ 2
x84 = u * x83
x85 = 4 * S0_t
x86 = Bp * S0
x87 = 2 * Sp
x88 = -x87
x89 = S * x11
x90 = S0_t ^ 2
x91 = Bpp - u * (Bp * (7 - x78) + x12 * (x22 - 5) + x69)
x92 = S ^ 2
x93 = u ^ 11 * x92
x94 = 18 * G
x95 = B * x27
x96 = x94 * x95
x97 = x93 * x96
x98 = S0 ^ 2
x99 = 18 * B * x47
x100 = x27 * x98 * x99
x101 = x50 * x92
x102 = x24 * x27
x103 = Bp * x101 * x102
x104 = x3 * x98
x105 = Bp * x104
x106 = x102 * x105
x107 = u ^ 10 * x92
x108 = x107 * x27
x109 = Gp * x74
x110 = x108 * x109
x111 = x40 * x98
x112 = x109 * x27
x113 = x111 * x112
x114 = 6 * G
x115 = Bp * x114
x116 = x108 * x115
x117 = Bp * x98
x118 = x114 * x40
x119 = x117 * x118 * x27
x120 = B * x42
x121 = S * x120
x122 = S0 * x27
x123 = 36 * G
x124 = x121 * x122 * x123
x125 = x98 * xi
x126 = x125 * x53
x127 = x123 * x126 * x95
x128 = xi ^ 2
x129 = x128 * x98
x130 = x129 * x37
x131 = x130 * x96
x132 = 4 * S0
x133 = S * x53
x134 = x132 * x133
x135 = Bp * x27
x136 = Gp * x134 * x135
x137 = 4 * x125
x138 = x27 * x40
x139 = Bp * Gp * x137 * x138
x140 = Bp * x46
x141 = x102 * x129 * x140
x142 = x129 * x53
x143 = x112 * x142
x144 = x115 * x142 * x27
x145 = Gp * S0
x146 = 12 * x38
x147 = x145 * x146 * x95
x148 = 12 * x46
x149 = x125 * x148
x150 = Gp * x149 * x95
x151 = 12 * G
x152 = x27 * x86
x153 = x151 * x152 * x38
x154 = G * x148
x155 = x125 * x135 * x154
x156 = S * x50
x157 = x27 * x56
x158 = B * x123 * x156 * x157
x159 = x132 * xi
x160 = Bp * x38
x161 = Gp * x159 * x160 * x27
x162 = 12 * x121
x163 = Gp * x157 * x162
x164 = Bp * x151 * x157 * x43
x165 = 4 * x3
x166 = Sp * u
x167 = 24 * S0
x168 = 2 * x36
x169 = 4 * u
x170 = 18 * x49
x171 = x46 * x98
x172 = S * x15
x173 = Gp * x151
x174 = x129 * x46
x175 = 36 * x49
x176 = S0 * x43
x177 = G * Gp
x178 = 24 * x125
x179 = x156 * x56
x180 = -Gp * x178 * x47 + 48 * S0 * x172 - Sp ^ 2 * x165 + Spp * x132 + 4 * Spp * x4 + x101 * x168 + x104 * x168 - x107 * x173 - x111 * x173 + 8 * x125 + x126 * x175 + x129 * x169 + x130 * x170 + x134 * x36 + x137 * x41 - x142 * x173 + x148 * x92 + x159 * x34 + x159 * x39 - x166 * x167 - x167 * x177 * x38 + x168 * x174 + x170 * x171 + x170 * x93 + x175 * x176 + x175 * x179 - 24 * x177 * x43 * x56 - 16 * x35 * x56 + 32 * x4 * x56
x181 = Gpp * x33
x182 = 2 * Gpp
x183 = 2 * x18
x184 = x14 * x183
x185 = x3 * x63
x186 = x184 * x185
x187 = -x20
x188 = x22 + 7
x189 = x22 + 5
x190 = Bpp + u * (-Bp * (x78 + 7) + x12 * x189 + x69)
x191 = 54 * S * x47
x192 = S0 * x15
x193 = x192 * x94
x194 = Gp * Sp * x165
x195 = x182 * x4
x196 = 10 * u * x145
x197 = Sp * x40
x198 = x151 * x197
x199 = 22 * Gp * x73
x200 = x27 * x68
x201 = x133 * x200
x202 = S0 * x200 * x3
x203 = 30 * x16 * x56
x204 = x1 * x181
x205 = B ^ 2
x206 = x205 * x27 * x42 * x51
x207 = 9 * x205
x208 = x122 * x207 * x46
x209 = x15 * x56
x210 = 14 * Gp * x209
x211 = x40 * x68
x212 = x157 * x211
x213 = x38 * x74
x214 = x135 * x213
x215 = x40 * x74
x216 = x152 * x215
x217 = x205 * x54
x218 = x157 * x217
x219 = Bp * x74
x220 = x219 * x27 * x57
x221 = Bpp * x27
x222 = 11 * Bp
x223 = x197 * x74
x224 = 5 * u
x225 = Bp * x133
x226 = Bp * x3
x227 = 9 * B
x228 = B * x46
x229 = S * x228
x230 = B * x56
x231 = 7 * x20
x232 = Bp * x56
x233 = x18 * x5
x234 = 15 * x21
x235 = Gp * x213 + S0 * x221 - S0 * x99 + x109 * x57 + x114 * x160 + x115 * x57 + x118 * x86 - x121 * x94 + 2 * x14 * x185 * x233 + x145 * x215 - x152 * x224 - x157 * x231 + x157 * x234 + x192 * x227 * x27 + x2 * x221 + x221 * x4 - x222 * x27 * x73 - x223 * x27 - x225 * x24 + x226 * x27 * x87 + 27 * x229 * x27 - x230 * x53 * x94 - x232 * x24 * x40 - x24 * x3 * x86
x236 = x30 * x6
x237 = -Sp
x238 = S * x8 + x237
x239 = -x15 * x238
x240 = S0 + x239
x241 = x132 * x240
x242 = x132 + 4 * x239
x243 = 4 * x240 ^ 2
x244 = S0_x * xi
x245 = S_x * x15
x246 = 4 * x240
x247 = x169 * x6
x248 = -4 * x166
x249 = x30 * x9
x250 = S0 * S0_tx - S0_t * S0_x
x251 = 1 / x98
x252 = x243 * x251
x253 = x83 * x9
x254 = 6 * x32
x255 = 6 * x31
x256 = B_x * x83
x257 = Spp + 6 * x172 + x248
x258 = x247 * x257
x259 = x30 * x7
x260 = u * x27
x261 = x259 * x260
x262 = 4 * x6
x263 = x13 * x84
x264 = 3 * x27
x265 = B_y * x259
x266 = 2 * x40
x267 = x266 * x62
x268 = x13 * x83
x269 = x15 * x268
x270 = x262 * x269
x271 = x268 * x40
x272 = x269 * x6
x273 = x266 * x7
x274 = x14 ^ 2
x275 = x274 * xi_x
x276 = x273 * x83
x277 = Gpp + x25 * (G * x8 + x60)
x278 = x249 * x277
x279 = x250 * x251
x280 = Bpp + x25 * (B * x8 + x10)
x281 = x280 * xi_x
x282 = x259 * xi_y
x283 = x138 * x274
x284 = G_y * x259
x285 = x13 * x183
x286 = x285 * x40
x287 = x273 * (x28 + x62)
x288 = x267 * x282
x289 = x261 * x280
x290 = x236 * x29
x291 = x274 * x279
x292 = 2 * x15
x293 = x290 * x292
x294 = S0 * S0_ty - S0_t * S0_y
x295 = x251 * x294
x296 = x15 * x33
x297 = S0_y * xi
x298 = x279 * x280
x299 = x286 * x62
x300 = x287 * x62
x301 = -S0 + x15 * x238
x302 = x301 * xi_x
x303 = x259 * x295
x304 = x267 * x303
x305 = x293 * x301
x306 = -x24 + x26 - x28
x307 = Bpp * x98
x308 = Bpp * u
x309 = Bp * x197
x310 = Bpp * x73
x311 = B * x133
x312 = x38 * x68
x313 = 15 * x40
x314 = x170 * x53
x315 = x33 * xi
x316 = x76 + x79
x317 = x22 * x77 + x81
x318 = u * x18
x319 = 18 * x205
x320 = x246 * x30
x321 = x30 * xi_y
x322 = S0_y * x30
x323 = S_y * x15
x324 = x169 * x236
x325 = x260 * x7
x326 = x273 * x62
x327 = x277 * x9
x328 = x249 * x83
x329 = x40 * x7
x330 = x184 * x329
x331 = x257 * x324
x332 = x138 * x7
x333 = x13 * x326
x334 = 4 * x236
x335 = x306 * x6
x336 = x269 * x334
x337 = x292 * x335
x338 = x266 * x274 * x83
x339 = x280 * x328
x340 = -Gp - x28 + x61
x341 = x301 * x336
AA[1,1] = x9
AA[1,2] = 0
BB[1,1] = x7 * (x14 * x19 + x23 + 8)
BB[1,2] = x29 * x32
CC[1,1] = -x100 - x103 - x106 + x110 + x113 + x116 + x119 - x124 - x127 - x131 - x136 - x139 - x141 + x143 + x144 + x147 + x150 + x153 + x155 - x158 - x161 + x163 + x164 + x180 + x5 * x83 * (2 * u * (S0 * x12 * (-x82 - 3) + x20 * (x71 * (11 - x78) + x88) + x22 * (x87 + x89 * (x21 - 3)) + x70 + x86 * (x80 + 5)) - 2 * x67) + x8 * x90 * (x41 - x48 + x55 - x65 - x84 * x91 + 2) + x85 * (S0 * x41 - S0 * x48 + S0 * x55 + x1 * x33 + x33 + x34 - 4 * x35 + x36 * x57 + x39 + 8 * x4 - x43 * x45 - x45 * x58 + x49 * x52 + 9 * x49 * x59 - x66 + x84 * (u * (-Bp * (S0 * (-x80 - 6) + x15 * (Sp - x72 + x75)) + x12 * (S0 * (-x82 - 4) + x15 * (Sp + x71 * (x22 - 7))) + x70) - x67)) - x97
CC[1,2] = x236 * (-x0 * (x182 - x186 - x190 * x27 + x8 * (Gp * (-x187 - x188) + x61 * (x187 + x189))) - x181 - x191 - x193 - x194 - x195 + x196 + x198 + x199 + x201 + x202 - x203 - x204 + x206 + x208 + x210 + x212 - x214 - x216 + x218 - x220 + x235)
SS[1] = B_y * x264 * x32 + Bp_x * x253 - Bp_y * x261 + G_x * x287 - G_y * x254 + Gp_y * x249 + S0_tx * x242 - S0_tx * x270 + S0_ty * x293 - S0_x * x262 * x263 - S0_x * (x159 + 12 * x172 + x248 + x85) + S0_y * x290 * x8 - S_x * x262 * x271 + S_y * x266 * x290 - Sp_x * x247 + x13 * x138 * x265 + x13 * x256 * x273 - x13 * x288 - x13 * x304 - x132 * x272 * xi_x + x241 * xi_x + x242 * x244 - x243 * xi_x - x244 * x270 + x245 * x246 + 8 * x245 * x6 + x246 * x272 * x279 - x250 * x252 - x253 * x281 - x253 * x298 - x255 * x256 + x258 * x279 + x258 * xi_x - x265 * x267 - x270 * x302 + x275 * x276 + x276 * x291 - x278 * x295 - x278 * xi_y + x279 * x300 + x282 * x283 + x282 * x299 + x283 * x303 + x284 * x286 + x289 * x295 + x289 * xi_y + x290 * x296 * xi_y + x293 * x297 + x295 * x305 + x299 * x303 + x300 * xi_x + x305 * xi_y
AA[2,1] = 0
AA[2,2] = x249
BB[2,1] = -x306 * x31
BB[2,2] = x259 * (x13 * x19 + x187 + x22 + 8)
CC[2,1] = -x6 * (x0 * (x182 + x186 + x27 * x91 + x8 * (-Gp * (x23 + 7) + x61 * (x23 + 5))) + x181 + x191 + x193 + x194 + x195 - x196 - x198 - x199 - x201 - x202 + x203 + x204 - x206 - x208 - x210 - x212 + x214 + x216 - x218 + x220 + x235)
CC[2,2] = x30 * (-B * Bp * x149 + B * x129 * x313 - B * x146 * x86 - 18 * Bp * S * x57 + Bpp * x33 * x4 + Bpp * x53 * x92 - S0 * x223 + 36 * S0 * x229 + 2 * S0_t * (Bpp * x209 - 6 * S0 * x20 + 12 * S0 * x21 + S0 * x211 + S0 * x217 + S0 * x308 + S0 * x314 - 3 * Sp * x228 + x1 * x132 + x132 - x140 * x51 - x145 * x154 + x156 * x170 + x170 * x59 - x173 * x43 - x173 * x58 + x205 * x52 + x207 * x59 - x219 * x43 - x219 * x58 - 7 * x226 * x56 + x230 * x313 + x309 + x310 + 21 * x311 + x312 + x315 * x36 * x46 + x318 * (u * (-Bp * (S0 * (x316 + 6) + x15 * (x237 + x72 + x75)) + x12 * (S0 * (x317 + 4) + x15 * (x188 * x71 + x237)) + x70) + x67) + x33 * x41 + 2 * x34 - 8 * x35 + 2 * x39 + 16 * x4 - x46 * x74 * x86 + x57 * x68 + 2 * x66) - Sp * x213 + Sp * x226 * x33 - Sp * x57 * x74 + u * x90 * (-Gp * x154 - x140 * x74 + x190 * x318 + x211 + x217 - x231 + x234 + x308 + x314 + 2 * x41 + 2 * x65 + 4) + 2 * x1 * x307 + x100 + x101 * x68 + x103 + x104 * x68 - 7 * x105 * x128 + x106 - x107 * x219 - x110 - x113 - x116 - x117 * x215 - x117 * x224 - x119 + 27 * x120 * x92 + x124 - 12 * x125 * x20 + 2 * x125 * x211 + x126 * x319 + x127 + x128 * x15 * x307 + x130 * x207 + x131 + x133 * x33 * x68 + x136 + x139 + x141 - x142 * x219 - x143 - x144 - x147 + x15 * x227 * x98 - x150 - x153 - x155 + x158 + x161 - x162 * x232 - x163 - x164 + x171 * x207 + x174 * x68 + x176 * x319 + x178 * x21 + x179 * x319 + x180 + x207 * x93 - x222 * x37 * x92 + x225 * x87 + x233 * (u * (-Bp * (S0 * (x316 + 5) + x15 * (x71 * (x78 + 11) + x88)) + x12 * (S0 * (x317 + 3) + x15 * (x88 + x89 * (x21 + 3))) + x70) + x67) + x307 + x309 * x315 + x310 * x315 + 42 * x311 * x56 + x312 * x315 - 16 * x73 * x86 + x97)
SS[2] = -B_x * x264 * x31 + B_x * x326 - B_x * x64 * x7 + B_y * x254 * x83 + Bp_x * x325 - Bp_y * x328 - G_x * x255 + G_x * x330 + Gp_x * x9 - S0_tx * x337 + S0_ty * x320 + S0_ty * x336 - S0_x * x335 * x8 + S0_y * x263 * x334 - S_x * x266 * x335 + S_y * x271 * x334 - Sp_y * x324 + x132 * x236 * x269 * xi_y + 8 * x236 * x323 + x240 * x279 * x337 + x241 * x321 - x243 * x321 - x244 * x337 + x246 * x322 * xi - x252 * x294 * x30 + x265 * x266 * x268 + x266 * x284 * x340 + x275 * x332 - x279 * x285 * x329 * x62 - x279 * x327 + x279 * x333 - x281 * x325 + x282 * x338 + x288 * x340 + x291 * x332 + x295 * x331 + x295 * x339 + x295 * x341 - x296 * x335 * xi_x + x297 * x336 - x298 * x325 - x302 * x337 + x303 * x338 + x304 * x340 + x320 * x323 - 4 * x322 * (S0_t - x166 + 3 * x172 + x56) - x327 * xi_x + x330 * x62 * xi_x + x331 * xi_y + x333 * xi_x + x339 * xi_y + x341 * xi_y

A11=AA[1,1]
A12=AA[1,2]
A21=AA[2,1]
A22=AA[2,2]
B11=BB[1,1]
B12=BB[1,2]
B21=BB[2,1]
B22=BB[2,2]
C11=CC[1,1]
C12=CC[1,1]
C21=CC[1,1]
C22=CC[1,1]
S22=SS[2]
S11=SS[1]

#println("Matrix C11 = $C11 and the source is $S0")
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

x0 = S0 * xi
x1 = u ^ 3
x2 = S * x1
x3 = S0 + S0_t * u + u * x0 + x2
x4 = B * x1
x5 = exp(x4)
x6 = u ^ 2
x7 = 8 * x6
x8 = x5 * x7
x9 = S * x6
x10 = S0_t + x0
x11 = 1 / u
x12 = S0 * x11 + x10 + x9
x13 = x12 ^ 3
x14 = 8 * S0
x15 = x5 * x6
x16 = 4 * S0
x17 = Sp * x6
x18 = S0 + x17 - 2 * x2
x19 = x12 ^ 2
x20 = x19 * x5
x21 = x18 * x20
x22 = S0 ^ 3
x23 = 1 / x22
x24 = 8 * S0_t * x22
x25 = S0 ^ 4
x26 = x23 * (x24 + 8 * x25 * xi)
x27 = G * x1
x28 = cosh(x27)
x29 = Fx * u
x30 = S0 ^ 2
x31 = 1 / x30
x32 = S0 * S0_tx
x33 = S0_t * S0_x
x34 = x32 - x33
x35 = x31 * x34
x36 = x29 + x35 + xi_x
x37 = Sp * x36 + St
x38 = S0 * xi_x + S0_tx + S0_x * x11 + S0_x * xi - x18 * x36 + x37 * x6
x39 = 4 * x6
x40 = x19 * x6
x41 = Fy * u
x42 = S0 * S0_ty
x43 = S0_t * S0_y
x44 = x42 - x43
x45 = x31 * x44
x46 = x41 + x45 + xi_y
x47 = Sh + Sp * x46
x48 = S0 * xi_y + S0_ty + S0_y * x11 + S0_y * xi - x18 * x46 + x47 * x6
x49 = sinh(x27)
x50 = x38 * x49
x51 = Fx * x6
x52 = -Fxp * u + x51
x53 = -x52
x54 = Fy * x6
x55 = -Fyp * u + x54
x56 = -x55
x57 = Fxh + Fxp * x46
x58 = Fyp * x36 + Fyt
x59 = 2 * x1
x60 = 2 * x6
x61 = Fy * x59 - Fyp * x60 + Fypp * u
x62 = Fx * x59 - Fxp * x60 + Fxpp * u
x63 = Gh + Gp * x46
x64 = x1 * x63
x65 = u ^ 4
x66 = 3 * x65
x67 = -G * x66 + Gp * x1
x68 = x46 * x67
x69 = x64 - x68
x70 = Gp * x36 + Gt
x71 = x1 * x70
x72 = 2 * x71
x73 = 6 * x65
x74 = -G * x73 + 2 * Gp * x1
x75 = x20 * x60
x76 = -12 * G * x65 + 4 * Gp * x1
x77 = 4 * S0_xx
x78 = x36 ^ 2
x79 = Spp * x36
x80 = Spt + x79
x81 = 2 * x32 - 2 * x33
x82 = 2 * x29 + 2 * xi_x
x83 = x31 * x81 + x82
x84 = S0_x - x37 * x59 + x6 * x80
x85 = 2 * S0_x * x23
x86 = Fxp * x36 + Fxt
x87 = u * x86 + x31 * (S0 * S0_txx - S0_t * S0_xx) - x34 * x85 - x36 * x53 + xi_xx
x88 = S * x73 - 4 * Sp * x1 + Spp * x6
x89 = -x36 * x88 + x84
x90 = Bp * x36 + Bt
x91 = -B * x66 + Bp * x1
x92 = x36 * x91
x93 = x1 * x90 - x92
x94 = x12 * x6
x95 = Gph + Gpp * x46
x96 = 2 * x1 * x95 - x63 * x73
x97 = 2 * x41 + 2 * xi_y
x98 = 2 * x45 + x97
x99 = Gpp * x36
x100 = Gpt + x99
x101 = 12 * u ^ 5
x102 = G * x101 - Gp * x73 + Gpp * x1
x103 = x1 * x100 - x102 * x36 - x66 * x70
x104 = x36 * x67
x105 = -x104 + x71
x106 = Bh + Bp * x46
x107 = x1 * x106
x108 = x107 - x46 * x91
x109 = S0_y * x23
x110 = S0_ty * S0_x
x111 = S0_tx * S0_y
x112 = S0 * S0_txy - S0_t * S0_xy
x113 = u * x57 + u * x58 - x109 * x81 + x31 * (-x110 + x111 + x112) + x31 * (x110 - x111 + x112) - x36 * x56 - x44 * x85 - x46 * x53 + 2 * xi_xy
x114 = 4 * x48
x115 = 8 * xi_y
x116 = 8 * S0_xy
x117 = 8 * S0_y
x118 = Sph + Spp * x46
x119 = x118 * x6
x120 = 2 * x35 + x82
x121 = -B * x73 + 2 * Bp * x1
x122 = 2 * u
x123 = Bpp * x36 + Bpt
x124 = B * x101 - Bp * x73 + Bpp * x1
x125 = 4 * S0_yy
x126 = x46 ^ 2
x127 = x31 * (2 * x42 - 2 * x43) + x97
x128 = S0_y + x119 - x47 * x59
x129 = 4 * x46
x130 = Fyh + Fyp * x46
x131 = u * x130 - 2 * x109 * x44 + x31 * (S0 * S0_tyy - S0_t * S0_yy) - x46 * x56 + xi_yy
x132 = Bph + Bpp * x46
ABCS[1] = 0
ABCS[2] = -x3 ^ 3 * x8
ABCS[3] = -x3 ^ 2 * x8 * (-Sp * u + x10 + 3 * x9)

if u==0.0
	ABCS[4] = 4*S0^3*Sp
else
	ABCS[4] = u * x13 * x14 * x5 + u * x21 * x26 - 12 * x12 ^ 4 * x15 + x13 * x15 * x26 + x16 * x21 + x18 * x23 * x40 * x5 * (4 * S0_t ^ 2 * x30 + 4 * S0_x ^ 2 - S0_xx * x16 + 4 * S0_y ^ 2 - S0_yy * x16 + x24 * xi + 4 * x25 * xi ^ 2) - x28 * x38 ^ 2 * x39 + x28 * x75 * (x103 * x98 + x105 * (-x108 - x55) + x113 * x67 + x36 * x96 - x59 * (Gc + x100 * x46 + x36 * x95 + x46 * x99) + x69 * (x1 * x90 - x52 - x92)) + x40 * (x28 * (2 * x105 ^ 2 + x120 * x62 + x120 * (x1 * x123 - x124 * x36 - x66 * x90) + x121 * x87 - x122 * (Fxpp * x36 + Fxpt) + x36 * (2 * x1 * x123 - x73 * x90) + x53 ^ 2 - x59 * (Bb - Bpp * x78 + x123 * x83) + x60 * x86 + 2 * x93 ^ 2 + x93 * (2 * Fxp * u - 2 * x51)) + x49 * (2 * x1 * (Gb - Gpp * x78 + x100 * x83) - x103 * x120 - x36 * (2 * x1 * x100 - x70 * x73) - x74 * x87 - (-2 * x104 + x72) * (2 * x1 * x90 - x121 * x36 - x52))) + x48 * x50 * x8 + x49 * x75 * (u * (Fxph + Fxpp * x46) + u * (Fypp * x36 + Fypt) - x36 * x61 - x46 * x62 - x53 * x56 - x57 * x6 - x58 * x6 - x69 * (-x36 * x74 + x72)) + x5 * x94 * (x28 * (-x105 * x114 - 4 * x38 * x69) + x49 * (-8 * S0_txy - S0_x * x115 - x11 * x116 + x113 * (x16 + 4 * x17 - 8 * x2) - x116 * xi - x117 * xi_x - x14 * xi_xy + x36 * (-16 * x1 * x47 + x117 + 8 * x119) - x7 * (Sc + x118 * x36 + x46 * x79 + x46 * x80) + x89 * (x115 + 8 * x41 + 8 * x45))) + x6 * (x12 * (x28 * (4 * S0_tyy + x108 * x114 + x11 * x125 + x117 * xi_y + x125 * xi - x128 * x129 - x129 * (x128 - x46 * x88) - 4 * x131 * x18 + x16 * xi_yy + x39 * (Spp * x126 + Ss + x118 * x127)) + x48 * x49 * (-x46 * x76 + 4 * x64)) + x19 * (x28 * (2 * x108 ^ 2 - x108 * (2 * Fyp * u - 2 * x54) - x121 * x131 - x122 * (Fyph + Fypp * x46) + x130 * x60 - x46 * (2 * x1 * x132 - x106 * x73) + x56 ^ 2 + x59 * (Bpp * x126 + Bs + x127 * x132) + x61 * x98 + 2 * x69 ^ 2 - x98 * (x1 * x132 - x106 * x66 - x124 * x46)) + x49 * (2 * x1 * (Gpp * x126 + Gs + x127 * x95) - x131 * x74 - x46 * x96 - x98 * (x1 * x95 - x102 * x46 - x63 * x66) + (2 * x64 - 2 * x68) * (2 * x107 - x121 * x46 + x55))) - 4 * x28 * x48 ^ 2) * exp(2 * x4) + x94 * (-x28 * (-4 * S0_txx - 8 * S0_x * xi_x - x11 * x77 - x16 * xi_xx + 4 * x18 * x87 + 4 * x36 * x84 + 4 * x36 * x89 + 4 * x38 * x93 - x39 * (Sb - Spp * x78 + x80 * x83) - x77 * xi) + x50 * (-x36 * x76 + 4 * x71))    
end
	

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

   

   x0 = u .^ 3
x1 = B .* x0
x2 = exp(x1)
x3 = 8 * x2
x4 = u .^ 2
x5 = S0_t .* u
x6 = u .* xi
x7 = S0 .* x6
x8 = S .* x0
x9 = S0 + x7 + x8
x10 = x5 + x9
x11 = x10 .^ 3
x12 = x11 .* x4
x13 = x10 .^ 2
x14 = G .* x0
x15 = tanh(x14)
x16 = 3 * u
x17 = G .* x16
x18 = x15 .* x4 .* (-Gp + x17)
x19 = S0 + 2 * x7
x20 = x11 .* (B .* x16 - Bp)
x21 = sech(x14)
x22 = S0 .* x21
x23 = -Fxp
x24 = Fx .* u
x25 = x23 + x24
x26 = -x25
x27 = S0 .^ 2
x28 = x6 + 1
x29 = 2 * x8
x30 = x27 .* x28 + x5 .* (S0 - x29)
x31 = S0_x .* x30
x32 = 2 * x4
x33 = S .* x32
x34 = S0_tx .* x33
x35 = 2 * S0
x36 = S .* x35 .* x4
x37 = x36 .* xi_x
x38 = St .* u
x39 = Fx .* S0
x40 = Fx .* x29
x41 = S0 .* x4
x42 = 4 * x26 .* (x31 + x41 .* (S0 .* (x38 - x39 + x40) + x34 + x37))
x43 = 2 * x24
x44 = x23 + x43
x45 = -x44
x46 = 2 * x5
x47 = S0_tx .* u
x48 = 5 * x4
x49 = x10 .* (S0 .* (-Fxpt .* x35 + 2 * S0 .* u .* x44 .* xi_x + S0 .* u .* (Fx .^ 2 .* x48 + Fxp .^ 2 - 4 * Fxp .* x24 + 2 * Fxt) - 2 * x45 .* x47) + S0_x .* x45 .* x46)
x50 = 6 * u
x51 = S0 .^ 4
x52 = S0 .^ 3
x53 = 2 * x0
x54 = x13 .* x2
x55 = x54 .* (S0_t .^ 2 .* x27 .* x4 + S0_x .^ 2 .* x4 - S0_xx .* x41 + S0_y .^ 2 .* x4 - S0_yy .* x41 + Sd .* x52 .* x53 + x28 .* x46 .* x52 + x4 .* x51 .* xi .^ 2 + 2 * x51 .* x6 + x51) ./ u
x56 = Fy .* u
x57 = S0_y .* x30
x58 = S0_ty .* x33
x59 = x36 .* xi_y
x60 = Sh .* u
x61 = Fy .* S0
x62 = Fy .* x29
x63 = -Fyp
x64 = 2 * x56 + x63
x65 = -x64
x66 = S0_ty .* u
x67 = (x10 .* (S0 .* (-Fyph .* x35 + 2 * S0 .* u .* x64 .* xi_y + S0 .* u .* (Fy .^ 2 .* x48 + 2 * Fyh + Fyp .^ 2 - 4 * Fyp .* x56) - 2 * x65 .* x66) + S0_y .* x46 .* x65) + (4 * Fyp - 4 * x56) .* (x41 .* (S0 .* (x60 - x61 + x62) + x58 + x59) + x57)) .* exp(2 * x1)
x68 = Fy .* Gp
x69 = -3 * Fyp .* G + x68
x70 = 3 * G
x71 = Fxp .* x70
x72 = Fx .* Gp - x71
x73 = Gp .* S0
x74 = Fyp .* S0
x75 = Gh .* S0
x76 = Fxh .* S0
x77 = Fyp .* S0_tx
x78 = Fyt .* S0
x79 = x39 .* x4
x80 = x4 .* x61
x81 = S0 .* u
x82 = x81 .* xi_y
x83 = x81 .* xi_x
x84 = x10 .* x2 .* x35
x85 = u ./ x52
x86 = 4 * x2
x87 = 4 * u
x88 = S0 .* sinh(x14)
x89 = x56 + x63
x90 = Fy .* x27
x91 = Fxp .* S0
x92 = u .^ 4
x93 = x61 .* x92
x94 = Fyp .* x39
x95 = Fxp .* x61
x96 = cosh(x14)
x97 = 3 * x1
x98 = Fxp .* x97
x99 = Bp .* x4
x100 = Fyp .* x97 + Fyp - x56 .* (x99 + 2)
x101 = Fy .* S0_tx
x102 = S0_ty .* x92
x103 = Bp .* x0
x104 = 3 * B
x105 = x104 .* x92
x106 = Fxp .* x104
x107 = u .^ 5
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x12 .* x3
BB[1,2] = 0
CC[1,1] = -u .* x13 .* x3 .* (-Sp .* x4 + x18 .* x9 + x19 + x5 .* (x18 + 2) + 4 * x8)
CC[1,2] = -x0 .* x15 .* x20 .* x3
SS[1] = 2 * x85 .* (x0 .* x21 .* x84 .* (S0 .* (Fxh .* x73 - Fxp .* S0_ty .* x17 - Fxp .* x75 + Fyp .* x70 .* x79 - Fyt .* x73 + Gp .* S0_ty .* x24 - Gt .* u .* x61 + Gt .* x74 - x17 .* x76 + x17 .* x77 + x17 .* x78 + x24 .* x75 - x47 .* x68 - x69 .* x83 - x71 .* x80 + x72 .* x82) + x5 .* (S0_x .* x69 - S0_y .* x72)) + x22 .* x42 + x22 .* x49 - x22 .* x67 - x55 .* (B .* x50 - 2 * Bp))
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x12 .* x86
CC[2,1] = x2 .* x20 .* x53 .* sinh(2 * x14)
CC[2,2] = -x54 .* x87 .* (x19 + x4 .* (S .* x87 - Sp) + x46)
SS[2] = -x85 .* (S0 .* x86 .* x96 .* (S0 .* x4 .* (Fx .* Fyp .* x27 + 4 * Fx .* S .* x93 - Fxp .* x58 + Fxp .* x90 - Fyp .* x34 + S0_tx .* x62 + S0_ty .* x40 + Sh .* x79 + St .* x80 + x25 .* x59 - x29 .* x94 - x29 .* x95 + x37 .* x89 - x38 .* x74 - x43 .* x90 - x60 .* x91) - x26 .* x57 + x31 .* x89) + x42 .* x88 + x49 .* x88 - x55 .* (-G .* x50 + 2 * Gp) + x67 .* x88 + x84 .* x96 .* (S0 .* (-Bh .* x0 .* x91 + Bh .* x39 .* x92 + Bp .* Fx .* x102 - Bp .* x101 .* x92 + Bt .* x0 .* x74 - Bt .* x93 - Fx .* S0_ty .* x32 - 5 * Fx .* x0 .* x61 - Fxh .* x81 - Fxp .* Fyp .* x81 + Fxp .* x66 + Fxph .* S0 + Fyp .* x47 + Fypt .* S0 - Fyt .* x81 + x100 .* x83 - x101 .* x32 - x102 .* x106 + x103 .* x76 - x103 .* x78 + x104 .* x107 .* x94 - x105 .* x76 + x105 .* x77 + x105 .* x78 - x106 .* x107 .* x61 + x32 .* x94 + x32 .* x95 + x82 .* (Bp .* Fx .* x0 - x44 - x98)) - x5 .* (S0_x .* x100 + S0_y .* (Fxp + x24 .* (x99 - 2) - x98))))

    
    nothing
end



function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
         S0, S0_x, S0_y, S0_t, S0_tx, S0_ty,u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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
x85 = G * x18
x86 = 12 * x85
x87 = x70 * x86
x88 = x24 * x85
x89 = Fy * xi_y
x90 = S * x70
x91 = 28 * x90
x92 = 20 * x34
x93 = x24 * xi
x94 = Gpp * x18
x95 = 8 * x94
x96 = x89 * x95
x97 = S * x58
x98 = x58 * xi
x99 = 16 * Spp
x100 = x86 * x90
x101 = Fyh * xi
x102 = 12 * x88
x103 = 8 * x97
x104 = 24 * x2
x105 = x104 * x65
x106 = G * Gd
x107 = x105 * x106
x108 = x22 * xi
x109 = x104 * x106
x110 = x17 * x65
x111 = x45 * x82
x112 = S * x65
x113 = x86 * xi_yy
x114 = x65 * x85
x115 = 4 * x33
x116 = S * x75
x117 = x33 * xi_y
x118 = 32 * x117
x119 = x65 * xi
x120 = S * x119
x121 = Sp * x54
x122 = S * x49
x123 = Sp * x104
x124 = S * x17
x125 = xi_y ^ 2
x126 = x116 * x20
x127 = x119 * x20
x128 = Bp * x60
x129 = xi ^ 2
x130 = x60 * x65
x131 = Gpp * x22
x132 = x18 * x60
x133 = 4 * x132
x134 = 42 * x124
x135 = x17 * x60
x136 = 8 * x63 * x71
x137 = x125 * x35
x138 = x0 * x125
x139 = 4 * x94
x140 = x63 * x75
x141 = 2 * x19
x142 = x129 * x22
x143 = xi ^ 3
x144 = xi ^ 4
x145 = x24 * x63
x146 = x12 * x23
x147 = 8 * Spp
x148 = u ^ 11
x149 = x124 * x148
x150 = 24 * B
x151 = Bh * Fy * x150
x152 = x75 * xi
x153 = x150 * x152
x154 = Bh * xi_y
x155 = x124 * x150 * x154
x156 = Bp * x89
x157 = x40 * x70
x158 = x119 * x40
x159 = 24 * x82
x160 = Fy * x149
x161 = x152 * x36
x162 = x84 * x86
x163 = S * x148
x164 = x85 * x89
x165 = x116 * xi
x166 = 32 * x165
x167 = 16 * x89
x168 = x167 * x94
x169 = x124 * x89
x170 = x36 * xi_y
x171 = x116 * x86
x172 = 72 * x106
x173 = x70 * xi
x174 = S * x173 * x2
x175 = x159 * xi_y
x176 = x71 * xi
x177 = S * x67
x178 = Sh * x45
x179 = S * x173
x180 = B ^ 2
x181 = 36 * x180
x182 = x181 * x75
x183 = u ^ 12
x184 = x67 * xi
x185 = x184 * x40
x186 = x63 * x67
x187 = 2 * x17
x188 = x187 * x28
x189 = x17 * x31
x190 = G ^ 2
x191 = 36 * x190
x192 = x191 * x75
x193 = x148 * x63
x194 = x129 * x65
x195 = 12 * x194
x196 = 72 * x85
x197 = x177 * x196
x198 = x60 * x85
x199 = 48 * x173
x200 = 8 * Gpp
x201 = x132 * x200
x202 = S * x71
x203 = x60 * xi
x204 = 30 * x202
x205 = Fyh * x38
x206 = 36 * x63
x207 = x2 * x67
x208 = x206 * x207
x209 = u ^ 13
x210 = x209 * x64
x211 = x2 * x42
x212 = u ^ 16
x213 = x212 * x66
x214 = x129 * x2
x215 = x106 * x214
x216 = 36 * x24
x217 = x144 * x70
x218 = x125 * x196
x219 = x125 * xi
x220 = x129 * x24
x221 = u ^ 15
x222 = x125 * x95
x223 = x47 * x51
x224 = Sp * x14
x225 = x119 * x51
x226 = x124 * xi
x227 = x183 * x226
x228 = x148 * xi
x229 = x177 * xi
x230 = x124 * x228
x231 = 18 * x180
x232 = x135 * x67
x233 = x125 * x71
x234 = 2 * x57
x235 = x183 * x63
x236 = x17 * x235
x237 = x129 * x71
x238 = 18 * x190
x239 = 15 * x63
x240 = 2 * x62
x241 = u ^ 14
x242 = x241 * x63
x243 = x36 * x69
x244 = x129 * x67
x245 = x209 * x63
x246 = x69 * x73
x247 = x129 * x75
x248 = 72 * x180
x249 = x169 * x183
x250 = x170 * x184
x251 = x36 * x83
x252 = x129 * x70
x253 = 72 * x190
x254 = 84 * S
x255 = x143 * x177
x256 = x2 * x63
x257 = x228 * x256
x258 = x241 * x64 * xi
x259 = x111 * x17
x260 = x11 * x18
x261 = x10 * x260
x262 = x124 * x181
x263 = x209 * x60
x264 = x135 * x228
x265 = x125 * x148
x266 = x17 * x219
x267 = 2 * x237
x268 = 54 * x245
x269 = 30 * x247
x270 = x124 * x191
x271 = x63 * x85
x272 = x125 * x129
x273 = x125 * x139
x274 = x169 * x209 * xi
x275 = exp(2 * x1)
x276 = Bb * x275
x277 = Gb * x275
x278 = Sb * x275
x279 = 2 * Gt
x280 = Bh * x2
x281 = x279 * x280
x282 = x18 * x41
x283 = Fx * Fy
x284 = x0 * x260
x285 = Fx * xi_y
x286 = 2 * Gp * x25
x287 = x146 * x18
x288 = x284 * xi
x289 = x18 * x2
x290 = Fxp * Fyp
x291 = 2 * x290
x292 = Gp * xi_xy
x293 = 8 * Sh
x294 = x170 * x206 * x221
x295 = x129 * x148 * x170
x296 = x203 * x241
x297 = x183 * x219
x298 = x183 * x206
x299 = 6 * B
x300 = x17 * x299
x301 = Fxt * x275
x302 = x22 * x301
x303 = x135 * x231
x304 = x212 * x63
x305 = x129 * x183
x306 = x17 * x231
x307 = x125 * x242
x308 = x272 * x67
x309 = Fx * x275
x310 = x309 * x32
x311 = Bt * x17
x312 = Bt * x275
x313 = Gt * x18
x314 = x313 * x44
x315 = Fxp * x309
x316 = St * x309
x317 = x135 * x238
x318 = x17 * x238
x319 = x18 * x277
x320 = St * x275
x321 = Gt * x2
x322 = x300 * x321
x323 = x11 * x17
x324 = Bh * Gt
x325 = Gp * x283
x326 = x18 * x283
x327 = x326 * x51
x328 = Fx * x2
x329 = 2 * xi_y
x330 = Gp * x329
x331 = x12 * x18 * x285
x332 = x124 * x65
x333 = Fxh * Gp
x334 = x11 * x333
x335 = x108 * x323
x336 = Fxh * x12 * x18
x337 = Fyt * Gp
x338 = x11 * x337
x339 = Fyt * x12 * x18
x340 = x292 * x50
x341 = Bt ^ 2
x342 = x275 * x59
x343 = Fx ^ 2
x344 = Fxp ^ 2
x345 = Gt ^ 2
x346 = x17 ^ 2
x347 = Bd * x346
x348 = x299 * x347
x349 = Bt * x309
x350 = x68 * x71
x351 = x309 * x313
x352 = x312 * x313
x353 = St * x312
x354 = 12 * G
x355 = Gt * x309
x356 = x301 * xi
x357 = 4 * x313 * x320
x358 = x321 * x68
x359 = 12 * Bh * G * x328
x360 = x152 * x17
x361 = Gp * x36
x362 = Fx * x361 * x50
x363 = Gp * x11
x364 = x202 * xi
x365 = B * x17
x366 = x275 * x343
x367 = x366 * x65
x368 = 18 * x367
x369 = x187 * x276
x370 = x126 * x275
x371 = x127 * x275
x372 = 2 * x319
x373 = B * x347
x374 = x105 * x373
x375 = x104 * x373
x376 = x187 * x256 * x67
x377 = x214 * x59
x378 = x292 * x323
x379 = x150 * x349
x380 = x152 * x349
x381 = 24 * x85
x382 = x349 * x381
x383 = x309 * x311
x384 = 24 * G * x355
x385 = x316 * x86
x386 = Fx * x11
x387 = x366 * x67
x388 = 36 * B
x389 = 72 * B
x390 = 48 * B
x391 = x366 * xi
x392 = x300 * x301
x393 = 20 * x193
x394 = 4 * x352
x395 = x275 * x341
x396 = x187 * x315
x397 = x366 * x85
x398 = x301 * x38
x399 = x275 * x345
x400 = x347 * x389
x401 = Gt * x299 * x36
x402 = x241 * x256
x403 = x214 * x67
x404 = Bh * x17 * x37
x405 = Fx * x330
x406 = S * x183 * xi
x407 = 2 * x236
x408 = x2 * x348
x409 = x214 * x373
x410 = x383 * x68
x411 = x351 * x68
x412 = S * x389
x413 = x228 * x397
x414 = x17 * x354 * x355
x415 = x268 * x366
x416 = x209 * x366
x417 = x17 * x228 * x366
x418 = x241 * x391
x419 = x306 * x366
x420 = x318 * x366
x421 = 3 * G * u
x422 = x18 * x6
x423 = u * x6
x424 = x1 * x6
x425 = 3 * x424
x426 = 2 * x4
x427 = x426 + 2
x428 = 9 * (G * x4 + G) ^ 2
x429 = x190 * x24
x430 = 9 * x429
x431 = x63 * (x430 + 4)
x432 = 6 * S * x5 * (3 * x429 + 2)
x433 = x6 ^ 2
x434 = 2 * x3
x435 = Gpp * x6
x436 = x10 * x37
x437 = 4 * Spp
x438 = 17 * x3
x439 = 9 * x190
x440 = x10 * x6
x441 = 5 * x3 + 3 * x4 + x425 + 3
x442 = Gp * x6
x443 = 9 * x4 + x438 + 7
x444 = 3 * x16
x445 = 36 * x116
x446 = 36 * x119
x447 = 36 * x229
ABCS[1] = x8 * x9
ABCS[2] = 8 * x10 * x8
ABCS[3] = u * x11 * x7
ABCS[4] = 12 * B * Fx * Gh * S * x148 * x17 * x2 + 12 * B * Fx * Gh * S * x17 * x183 * x2 * xi + 6 * B * Fx * Gh * x129 * x17 * x2 * x67 + 6 * B * Fx * Gh * x17 * x2 * x241 * x63 + 6 * B * Fx * Gh * x17 * x2 * x70 + 12 * B * Fx * Gh * x17 * x2 * x75 * xi + 144 * B * Fy * G * S * x18 * x183 * xi_y + 144 * B * Fy * G * S * x18 * x209 * xi * xi_y + 72 * B * Fy * G * x129 * x148 * x18 * xi_y + 72 * B * Fy * G * x18 * x221 * x63 * xi_y + 144 * B * Fy * G * x18 * x67 * xi * xi_y + 72 * B * Fy * G * x18 * x75 * xi_y + 24 * B * Fy * Gh * S * x148 * x18 + 24 * B * Fy * Gh * S * x18 * x183 * xi + 12 * B * Fy * Gh * x129 * x18 * x67 + 12 * B * Fy * Gh * x18 * x241 * x63 + 12 * B * Fy * Gh * x18 * x70 + 24 * B * Fy * Gh * x18 * x75 * xi + 12 * B * Fy * S * Sh * x148 * x17 + 156 * B * Fy * S * x17 * x67 * xi * xi_y + 144 * B * Fy * S * x17 * x75 * xi_y + 12 * B * Fy * Sh * x17 * x70 + 12 * B * Fy * Sh * x17 * x75 * xi + 54 * B * Fy * x129 * x17 * x70 * xi_y + 102 * B * Fy * x17 * x183 * x63 * xi_y + 42 * B * Fy * x17 * x24 * xi_y + 96 * B * Fy * x17 * x65 * xi * xi_y + 12 * B * Fyh * S * x17 * x70 + 12 * B * Fyh * S * x17 * x75 * xi + 6 * B * Fyh * x129 * x17 * x65 + 6 * B * Fyh * x148 * x17 * x63 + 6 * B * Fyh * x17 * x22 + 12 * B * Fyh * x17 * x24 * xi + 72 * B * G * S * x125 * x148 * x18 + 72 * B * G * S * x125 * x18 * x183 * xi + 72 * B * G * S * x18 * x209 * x60 + 72 * B * G * S * x18 * x241 * x60 * xi + 36 * B * G * x125 * x129 * x18 * x67 + 36 * B * G * x125 * x18 * x241 * x63 + 36 * B * G * x125 * x18 * x70 + 72 * B * G * x125 * x18 * x75 * xi + 36 * B * G * x129 * x18 * x183 * x60 + 72 * B * G * x148 * x18 * x60 * xi + 36 * B * G * x18 * x212 * x60 * x63 + 36 * B * G * x18 * x60 * x67 + 24 * B * Gh * S * x148 * x18 * xi * xi_y + 24 * B * Gh * S * x18 * x67 * xi_y + 12 * B * Gh * x129 * x18 * x75 * xi_y + 12 * B * Gh * x18 * x209 * x63 * xi_y + 12 * B * Gh * x18 * x65 * xi_y + 24 * B * Gh * x18 * x70 * xi * xi_y + 12 * B * S * Sh * x17 * x67 * xi_y + 72 * B * S * x125 * x17 * x70 + 72 * B * S * x125 * x17 * x75 * xi + 84 * B * S * x148 * x17 * x60 * xi + 72 * B * S * x17 * x60 * x67 + 12 * B * S * x17 * x65 * xi_yy + 12 * B * S * x17 * x70 * xi * xi_yy + 12 * B * Sh * x17 * x65 * xi_y + 12 * B * Sh * x17 * x70 * xi * xi_y + 24 * B * x125 * x129 * x17 * x65 + 48 * B * x125 * x148 * x17 * x63 + 24 * B * x125 * x17 * x22 + 48 * B * x125 * x17 * x24 * xi + 6 * B * x129 * x17 * x24 * xi_yy + 30 * B * x129 * x17 * x60 * x75 + 54 * B * x17 * x209 * x60 * x63 + 12 * B * x17 * x22 * xi * xi_yy + 6 * B * x17 * x23 * xi_yy + 18 * B * x17 * x60 * x65 + 48 * B * x17 * x60 * x70 * xi + 6 * B * x17 * x63 * x67 * xi_yy - B * x206 * x212 * x397 - 84 * B * x230 * x366 + 24 * Bd * Bp * S * x129 * x2 * x346 * x70 + 8 * Bd * Bp * S * x143 * x2 * x346 * x75 + 8 * Bd * Bp * S * x2 * x24 * x346 + 24 * Bd * Bp * S * x2 * x346 * x65 * xi + 2 * Bd * Bp * x0 * x2 * x346 + 12 * Bd * Bp * x129 * x148 * x2 * x346 * x63 + 12 * Bd * Bp * x129 * x2 * x22 * x346 + 8 * Bd * Bp * x143 * x2 * x24 * x346 + 2 * Bd * Bp * x144 * x2 * x346 * x65 + 8 * Bd * Bp * x183 * x2 * x346 * x64 + 8 * Bd * Bp * x2 * x209 * x346 * x64 * xi + 2 * Bd * Bp * x2 * x221 * x346 * x66 + 8 * Bd * Bp * x2 * x23 * x346 * xi + 24 * Bd * Bp * x2 * x346 * x63 * x67 * xi + 12 * Bd * Bp * x2 * x346 * x63 * x75 + 24 * Bh * Fy * G * S * x148 * x18 + 24 * Bh * Fy * G * S * x18 * x183 * xi + 12 * Bh * Fy * G * x129 * x18 * x67 + 12 * Bh * Fy * G * x18 * x241 * x63 + 12 * Bh * Fy * G * x18 * x70 + 24 * Bh * Fy * G * x18 * x75 * xi + 28 * Bh * Fy * S * x17 * x70 + 32 * Bh * Fy * S * x17 * x75 * xi + 12 * Bh * Fy * x129 * x17 * x65 + 20 * Bh * Fy * x148 * x17 * x63 + 8 * Bh * Fy * x17 * x22 + 20 * Bh * Fy * x17 * x24 * xi + 24 * Bh * G * S * x148 * x18 * xi * xi_y + 24 * Bh * G * S * x18 * x67 * xi_y + 12 * Bh * G * x129 * x18 * x75 * xi_y + 12 * Bh * G * x18 * x209 * x63 * xi_y + 12 * Bh * G * x18 * x65 * xi_y + 24 * Bh * G * x18 * x70 * xi * xi_y + 8 * Bh * Gh * S * x18 * x67 * xi + 8 * Bh * Gh * S * x18 * x75 + 4 * Bh * Gh * x129 * x18 * x70 + 4 * Bh * Gh * x18 * x183 * x63 + 4 * Bh * Gh * x18 * x24 + 8 * Bh * Gh * x18 * x65 * xi + 4 * Bh * S * Sh * x17 * x75 + 32 * Bh * S * x17 * x65 * xi_y + 32 * Bh * S * x17 * x70 * xi * xi_y + 4 * Bh * Sh * x17 * x24 + 4 * Bh * Sh * x17 * x65 * xi + 12 * Bh * x129 * x17 * x24 * xi_y - Bh * x153 * x36 + 24 * Bh * x17 * x22 * xi * xi_y + 12 * Bh * x17 * x23 * xi_y + 20 * Bh * x17 * x63 * x67 * xi_y + 4 * Bp * Fxt * S * x17 * x275 * x65 + 4 * Bp * Fxt * S * x17 * x275 * x70 * xi + 2 * Bp * Fxt * x129 * x17 * x24 * x275 + 4 * Bp * Fxt * x17 * x22 * x275 * xi + 2 * Bp * Fxt * x17 * x23 * x275 + 2 * Bp * Fxt * x17 * x275 * x63 * x67 + 4 * Bp * S * x17 * x275 * x343 * x67 * xi + 4 * Bp * S * x17 * x275 * x343 * x75 - Bp * x101 * x157 + 2 * Bp * x129 * x17 * x275 * x343 * x70 - Bp * x129 * x30 * x74 - 2 * Bp * x135 * x235 + 2 * Bp * x17 * x183 * x275 * x343 * x63 + 2 * Bp * x17 * x24 * x275 * x343 + 4 * Bp * x17 * x275 * x343 * x65 * xi - Bp * x20 * x39 * xi + 16 * Bpp * Fy * S * x17 * x65 * xi_y + 16 * Bpp * Fy * S * x17 * x70 * xi * xi_y + 8 * Bpp * Fy * x129 * x17 * x24 * xi_y + 16 * Bpp * Fy * x17 * x22 * xi * xi_y + 8 * Bpp * Fy * x17 * x23 * xi_y + 8 * Bpp * Fy * x17 * x63 * x67 * xi_y + 8 * Bpp * S * x125 * x17 * x24 + 8 * Bpp * S * x125 * x17 * x65 * xi + 8 * Bpp * S * x17 * x60 * x70 + 8 * Bpp * S * x17 * x60 * x75 * xi + 4 * Bpp * x0 * x125 * x17 + 4 * Bpp * x125 * x129 * x17 * x22 + 8 * Bpp * x125 * x17 * x23 * xi + 4 * Bpp * x125 * x17 * x63 * x75 + 4 * Bpp * x129 * x17 * x60 * x65 + 4 * Bpp * x148 * x17 * x60 * x63 + 4 * Bpp * x17 * x22 * x60 + 8 * Bpp * x17 * x24 * x60 * xi + 4 * Bs * S * x17 * x24 + 4 * Bs * S * x17 * x65 * xi + 2 * Bs * x0 * x17 + 2 * Bs * x129 * x17 * x22 + 4 * Bs * x17 * x23 * xi + 2 * Bs * x17 * x63 * x75 + 12 * Bt * Fy * G * S * x148 * x17 * x2 + 12 * Bt * Fy * G * S * x17 * x183 * x2 * xi + 6 * Bt * Fy * G * x129 * x17 * x2 * x67 + 6 * Bt * Fy * G * x17 * x2 * x241 * x63 + 6 * Bt * Fy * G * x17 * x2 * x70 + 12 * Bt * Fy * G * x17 * x2 * x75 * xi + 12 * Bt * G * S * x148 * x17 * x2 * xi * xi_y + 12 * Bt * G * S * x17 * x2 * x67 * xi_y + 6 * Bt * G * x129 * x17 * x2 * x75 * xi_y + 6 * Bt * G * x17 * x2 * x209 * x63 * xi_y + 6 * Bt * G * x17 * x2 * x65 * xi_y + 12 * Bt * G * x17 * x2 * x70 * xi * xi_y + 4 * Bt * Gh * S * x17 * x2 * x67 * xi + 4 * Bt * Gh * S * x17 * x2 * x75 + 2 * Bt * Gh * x129 * x17 * x2 * x70 + 2 * Bt * Gh * x17 * x183 * x2 * x63 + 2 * Bt * Gh * x17 * x2 * x24 + 4 * Bt * Gh * x17 * x2 * x65 * xi - Bt * St * x275 * x78 + 168 * Fx * Fy * G * S * x148 * x17 * x2 * xi + 144 * Fx * Fy * G * S * x17 * x2 * x67 + 60 * Fx * Fy * G * x129 * x17 * x2 * x75 + 108 * Fx * Fy * G * x17 * x2 * x209 * x63 + 36 * Fx * Fy * G * x17 * x2 * x65 + 96 * Fx * Fy * G * x17 * x2 * x70 * xi + 16 * Fx * Fy * Gpp * S * x17 * x2 * x70 + 16 * Fx * Fy * Gpp * S * x17 * x2 * x75 * xi + 8 * Fx * Fy * Gpp * x129 * x17 * x2 * x65 + 8 * Fx * Fy * Gpp * x148 * x17 * x2 * x63 + 8 * Fx * Fy * Gpp * x17 * x2 * x22 + 16 * Fx * Fy * Gpp * x17 * x2 * x24 * xi + 16 * Fx * Fy * S * Spp * x18 * x2 * x70 + 72 * Fx * Fy * S * x18 * x190 * x2 * x209 + 72 * Fx * Fy * S * x18 * x190 * x2 * x241 * xi + 84 * Fx * Fy * S * x18 * x2 * x65 + 60 * Fx * Fy * S * x18 * x2 * x70 * xi + 16 * Fx * Fy * Spp * x18 * x2 * x22 + 16 * Fx * Fy * Spp * x18 * x2 * x24 * xi + 36 * Fx * Fy * x129 * x18 * x183 * x190 * x2 + 72 * Fx * Fy * x148 * x18 * x190 * x2 * xi + 36 * Fx * Fy * x18 * x190 * x2 * x212 * x63 + 36 * Fx * Fy * x18 * x190 * x2 * x67 + 30 * Fx * Fy * x18 * x2 * x63 * x67 + 4 * Fx * Fyp * S * x18 * x2 * x24 + 4 * Fx * Fyp * S * x18 * x2 * x65 * xi + 2 * Fx * Fyp * x0 * x18 * x2 + 2 * Fx * Fyp * x129 * x18 * x2 * x22 + 4 * Fx * Fyp * x18 * x2 * x23 * xi + 2 * Fx * Fyp * x18 * x2 * x63 * x75 + 24 * Fx * G * Gh * S * x148 * x18 * x2 + 24 * Fx * G * Gh * S * x18 * x183 * x2 * xi + 12 * Fx * G * Gh * x129 * x18 * x2 * x67 + 12 * Fx * G * Gh * x18 * x2 * x241 * x63 + 12 * Fx * G * Gh * x18 * x2 * x70 + 24 * Fx * G * Gh * x18 * x2 * x75 * xi + 12 * Fx * G * S * Sh * x148 * x17 * x2 + 156 * Fx * G * S * x17 * x2 * x67 * xi * xi_y + 144 * Fx * G * S * x17 * x2 * x75 * xi_y + 12 * Fx * G * Sh * x17 * x2 * x70 + 12 * Fx * G * Sh * x17 * x2 * x75 * xi + 54 * Fx * G * x129 * x17 * x2 * x70 * xi_y + 102 * Fx * G * x17 * x183 * x2 * x63 * xi_y + 42 * Fx * G * x17 * x2 * x24 * xi_y + 96 * Fx * G * x17 * x2 * x65 * xi * xi_y + 28 * Fx * Gh * S * x17 * x2 * x70 + 32 * Fx * Gh * S * x17 * x2 * x75 * xi + 12 * Fx * Gh * x129 * x17 * x2 * x65 + 20 * Fx * Gh * x148 * x17 * x2 * x63 + 8 * Fx * Gh * x17 * x2 * x22 + 20 * Fx * Gh * x17 * x2 * x24 * xi + 16 * Fx * Gpp * S * x17 * x2 * x65 * xi_y + 16 * Fx * Gpp * S * x17 * x2 * x70 * xi * xi_y + 8 * Fx * Gpp * x129 * x17 * x2 * x24 * xi_y + 16 * Fx * Gpp * x17 * x2 * x22 * xi * xi_y + 8 * Fx * Gpp * x17 * x2 * x23 * xi_y + 8 * Fx * Gpp * x17 * x2 * x63 * x67 * xi_y + 16 * Fx * S * Spp * x18 * x2 * x65 * xi_y + 72 * Fx * S * x18 * x183 * x190 * x2 * xi_y + 72 * Fx * S * x18 * x190 * x2 * x209 * xi * xi_y + 68 * Fx * S * x18 * x2 * x24 * xi_y + 56 * Fx * S * x18 * x2 * x65 * xi * xi_y + 24 * Fx * Sh * x18 * x2 * x22 + 16 * Fx * Sh * x18 * x2 * x24 * xi + 16 * Fx * Spp * x18 * x2 * x22 * xi * xi_y + 16 * Fx * Spp * x18 * x2 * x23 * xi_y + 36 * Fx * x129 * x148 * x18 * x190 * x2 * xi_y + 36 * Fx * x18 * x190 * x2 * x221 * x63 * xi_y + 72 * Fx * x18 * x190 * x2 * x67 * xi * xi_y + 36 * Fx * x18 * x190 * x2 * x75 * xi_y + 24 * Fx * x18 * x2 * x63 * x75 * xi_y - Fx * x18 * x214 * x24 * x30 - Fx * x280 * x37 * x71 - Fx * x402 * x404 - Fx * x403 * x404 + 12 * Fxh * G * S * x17 * x2 * x70 + 12 * Fxh * G * S * x17 * x2 * x75 * xi + 6 * Fxh * G * x129 * x17 * x2 * x65 + 6 * Fxh * G * x148 * x17 * x2 * x63 + 6 * Fxh * G * x17 * x2 * x22 + 12 * Fxh * G * x17 * x2 * x24 * xi + 4 * Fxh * S * x18 * x2 * x22 + 8 * Fxh * S * x18 * x2 * x24 * xi + 8 * Fxh * x18 * x2 * x63 * x70 - Fxh * x2 * x286 - Fxh * x261 - Fxh * x287 - Fxh * x288 + 4 * Fxp * Fy * S * x18 * x2 * x24 + 4 * Fxp * Fy * S * x18 * x2 * x65 * xi + 2 * Fxp * Fy * x0 * x18 * x2 + 2 * Fxp * Fy * x129 * x18 * x2 * x22 + 4 * Fxp * Fy * x18 * x2 * x23 * xi + 2 * Fxp * Fy * x18 * x2 * x63 * x75 + 4 * Fxt * Gp * S * x18 * x275 * x65 + 4 * Fxt * Gp * S * x18 * x275 * x70 * xi + 2 * Fxt * Gp * x129 * x18 * x24 * x275 + 4 * Fxt * Gp * x18 * x22 * x275 * xi + 2 * Fxt * Gp * x18 * x23 * x275 + 2 * Fxt * Gp * x18 * x275 * x63 * x67 + 4 * Fxt * S * Sp * x17 * x275 * x65 + 4 * Fxt * Sp * x17 * x22 * x275 * xi + 4 * Fxt * Sp * x17 * x23 * x275 + 4 * Fxt * x0 * x17 * x275 * xi + 4 * Fxt * x10 * x17 * x275 + 24 * Fy * G * Gt * S * x148 * x18 * x2 + 24 * Fy * G * Gt * S * x18 * x183 * x2 * xi + 12 * Fy * G * Gt * x129 * x18 * x2 * x67 + 12 * Fy * G * Gt * x18 * x2 * x241 * x63 + 12 * Fy * G * Gt * x18 * x2 * x70 + 24 * Fy * G * Gt * x18 * x2 * x75 * xi + 12 * Fy * G * S * St * x148 * x17 * x2 + 12 * Fy * G * St * x17 * x2 * x70 + 12 * Fy * G * St * x17 * x2 * x75 * xi + 4 * Fy * Gp * S * x18 * x70 * xi_y + 4 * Fy * Gp * S * x18 * x75 * xi * xi_y + 2 * Fy * Gp * x129 * x18 * x65 * xi_y + 2 * Fy * Gp * x148 * x18 * x63 * xi_y + 2 * Fy * Gp * x18 * x22 * xi_y + 4 * Fy * Gp * x18 * x24 * xi * xi_y + 28 * Fy * Gt * S * x17 * x2 * x70 + 32 * Fy * Gt * S * x17 * x2 * x75 * xi + 12 * Fy * Gt * x129 * x17 * x2 * x65 + 20 * Fy * Gt * x148 * x17 * x2 * x63 + 8 * Fy * Gt * x17 * x2 * x22 + 20 * Fy * Gt * x17 * x2 * x24 * xi + 4 * Fy * S * Sp * x17 * x70 * xi_y + 4 * Fy * Sp * x17 * x22 * xi_y + 4 * Fy * Sp * x17 * x24 * xi * xi_y + 24 * Fy * St * x18 * x2 * x22 + 16 * Fy * St * x18 * x2 * x24 * xi + 4 * Fy * x0 * x17 * xi_y - Fy * x159 * x227 + 4 * Fy * x17 * x23 * xi * xi_y - Fy * x227 * x321 * x68 - Fy * x55 * x99 - Fy * x56 * x98 + 4 * Fyh * Gp * S * x18 * x65 + 4 * Fyh * Gp * S * x18 * x70 * xi + 2 * Fyh * Gp * x129 * x18 * x24 + 4 * Fyh * Gp * x18 * x22 * xi + 2 * Fyh * Gp * x18 * x23 + 2 * Fyh * Gp * x18 * x63 * x67 + 4 * Fyh * S * Sp * x17 * x65 + 4 * Fyh * Sp * x17 * x22 * xi + 4 * Fyh * Sp * x17 * x23 + 4 * Fyh * x0 * x17 * xi + 4 * Fyh * x10 * x17 - Fyh * x100 - Fyh * x136 + 12 * Fyt * G * S * x17 * x2 * x70 + 12 * Fyt * G * S * x17 * x2 * x75 * xi + 6 * Fyt * G * x129 * x17 * x2 * x65 + 6 * Fyt * G * x148 * x17 * x2 * x63 + 6 * Fyt * G * x17 * x2 * x22 + 12 * Fyt * G * x17 * x2 * x24 * xi + 4 * Fyt * S * x18 * x2 * x22 + 8 * Fyt * S * x18 * x2 * x24 * xi + 8 * Fyt * x18 * x2 * x63 * x70 - Fyt * x2 * x286 - Fyt * x261 - Fyt * x287 - Fyt * x288 + 24 * G * Gt * S * x148 * x18 * x2 * xi * xi_y + 24 * G * Gt * S * x18 * x2 * x67 * xi_y + 12 * G * Gt * x129 * x18 * x2 * x75 * xi_y + 12 * G * Gt * x18 * x2 * x209 * x63 * xi_y + 12 * G * Gt * x18 * x2 * x65 * xi_y + 24 * G * Gt * x18 * x2 * x70 * xi * xi_y + 12 * G * S * St * x17 * x2 * x67 * xi_y + 24 * G * S * x17 * x2 * x65 * xi_xy + 24 * G * S * x17 * x2 * x70 * xi * xi_xy + 12 * G * St * x17 * x2 * x65 * xi_y + 12 * G * St * x17 * x2 * x70 * xi * xi_y + 12 * G * x129 * x17 * x2 * x24 * xi_xy + 24 * G * x17 * x2 * x22 * xi * xi_xy + 12 * G * x17 * x2 * x23 * xi_xy + 12 * G * x17 * x2 * x63 * x67 * xi_xy + 8 * Gc * S * x17 * x2 * x24 + 8 * Gc * S * x17 * x2 * x65 * xi + 4 * Gc * x0 * x17 * x2 + 4 * Gc * x129 * x17 * x2 * x22 + 8 * Gc * x17 * x2 * x23 * xi + 4 * Gc * x17 * x2 * x63 * x75 + 24 * Gd * Gp * S * x129 * x2 * x70 + 8 * Gd * Gp * S * x143 * x2 * x75 + 8 * Gd * Gp * S * x2 * x24 + 24 * Gd * Gp * S * x2 * x65 * xi + 2 * Gd * Gp * x0 * x2 + 12 * Gd * Gp * x129 * x148 * x2 * x63 + 12 * Gd * Gp * x129 * x2 * x22 + 8 * Gd * Gp * x143 * x2 * x24 + 2 * Gd * Gp * x144 * x2 * x65 + 8 * Gd * Gp * x183 * x2 * x64 + 8 * Gd * Gp * x2 * x209 * x64 * xi + 2 * Gd * Gp * x2 * x221 * x66 + 8 * Gd * Gp * x2 * x23 * xi + 24 * Gd * Gp * x2 * x63 * x67 * xi + 12 * Gd * Gp * x2 * x63 * x75 + 8 * Gh * Gt * S * x18 * x2 * x67 * xi + 8 * Gh * Gt * S * x18 * x2 * x75 + 4 * Gh * Gt * x129 * x18 * x2 * x70 + 4 * Gh * Gt * x18 * x183 * x2 * x63 + 4 * Gh * Gt * x18 * x2 * x24 + 8 * Gh * Gt * x18 * x2 * x65 * xi + 4 * Gh * S * St * x17 * x2 * x75 + 4 * Gh * St * x17 * x2 * x24 + 4 * Gh * St * x17 * x2 * x65 * xi + 4 * Gp * S * x18 * x24 * xi_yy + 4 * Gp * S * x18 * x275 * x343 * x67 * xi + 4 * Gp * S * x18 * x275 * x343 * x75 + 4 * Gp * S * x18 * x60 * x67 * xi + 4 * Gp * S * x18 * x60 * x75 + 4 * Gp * S * x18 * x65 * xi * xi_yy + 2 * Gp * x0 * x18 * xi_yy + 2 * Gp * x129 * x18 * x22 * xi_yy + 2 * Gp * x129 * x18 * x275 * x343 * x70 + 2 * Gp * x129 * x18 * x60 * x70 - Gp * x165 * x386 * x73 + 2 * Gp * x18 * x183 * x275 * x343 * x63 + 2 * Gp * x18 * x183 * x60 * x63 + 4 * Gp * x18 * x23 * xi * xi_yy + 2 * Gp * x18 * x24 * x275 * x343 + 2 * Gp * x18 * x24 * x60 + 4 * Gp * x18 * x275 * x343 * x65 * xi + 4 * Gp * x18 * x60 * x65 * xi + 2 * Gp * x18 * x63 * x75 * xi_yy - Gpp * x133 * x193 + 4 * Gt * S * Sh * x17 * x2 * x75 + 32 * Gt * S * x17 * x2 * x65 * xi_y + 32 * Gt * S * x17 * x2 * x70 * xi * xi_y + 4 * Gt * Sh * x17 * x2 * x24 + 4 * Gt * Sh * x17 * x2 * x65 * xi - Gt * x124 * x207 * x68 * xi_y + 12 * Gt * x129 * x17 * x2 * x24 * xi_y + 24 * Gt * x17 * x2 * x22 * xi * xi_y + 12 * Gt * x17 * x2 * x23 * xi_y + 20 * Gt * x17 * x2 * x63 * x67 * xi_y - Gt * x209 * x256 * x300 * xi_y + 8 * S * Sc * x18 * x2 * x24 + 16 * S * Sd * x129 * x2 * x24 + 16 * S * Sd * x2 * x22 * xi + 4 * S * Sp * x17 * x24 * xi_yy + 4 * S * Sp * x17 * x275 * x343 * x75 + 4 * S * Sp * x17 * x60 * x75 + 16 * S * u * x2 + 72 * S * x0 * x129 * x2 + 56 * S * x10 * x2 * xi - S * x107 - S * x123 * x142 + 40 * S * x143 * x2 * x23 + 8 * S * x144 * x2 * x22 + 2 * S * x17 * x22 * x275 * x344 + 2 * S * x17 * x22 * x61 + 2 * S * x17 * x24 * x275 * x344 * xi + 2 * S * x17 * x24 * x61 * xi + 16 * S * x18 * x2 * x22 * xi * xi_xy + 16 * S * x18 * x2 * x23 * xi_xy - S * x22 * x260 * x290 - S * x25 * x52 - S * x374 + 8 * Sc * x0 * x18 * x2 + 8 * Sc * x18 * x2 * x23 * xi + 24 * Sd * x2 * x63 * x65 + 32 * Sd * x2 * x63 * x70 * xi + 16 * Sd * x2 * x64 * x67 + 4 * Sh ^ 2 * x17 * x24 - Sh * x114 * x45 - Sh * x115 * x116 - Sh * x115 * x119 - Sh * x33 * x44 - Sh * x35 * x36 + 4 * Sp * x0 * x17 * xi_yy + 4 * Sp * x17 * x23 * xi * xi_yy + 4 * Sp * x17 * x24 * x275 * x343 + 4 * Sp * x17 * x24 * x60 + 4 * Sp * x17 * x275 * x343 * x65 * xi + 4 * Sp * x17 * x60 * x65 * xi - 8 * Sp * x282 * xi * xi_xy - Sp * x53 - Spp * x135 * x32 + 4 * St ^ 2 * x17 * x24 * x275 + 16 * St * x18 * x2 * x22 * xi * xi_y + 16 * St * x18 * x2 * x23 * xi_y - St * x24 * x289 * x293 - u * x14 - x0 * x143 * x15 + 2 * x0 * x17 * x275 * x344 * xi + 2 * x0 * x17 * x61 * xi - x0 * x18 * x51 * xi_xy - x0 * x20 * x278 - x0 * x21 + 2 * x0 * x275 * x6 * xi_xx * (x17 * (-3 * B * x423 + Bp * x6 - 4 * S * u + 2 * Sp) + x422 * (Gp - x421)) - x10 * x123 * x129 + x10 * x17 * x275 * x344 + x10 * x17 * x61 - x10 * x289 * x291 - x10 * x54 * xi - x100 * x301 - x101 * x102 - x101 * x103 - x101 * x171 - x102 * x356 - x103 * x356 - 72 * x106 * x116 * x214 - x106 * x208 - x107 * x143 - x108 * x109 - x108 * x113 - x108 * x170 * x99 - 12 * x108 * x283 * x289 - x108 * x336 - x108 * x339 - x108 * x375 - x108 * x56 * x73 - x109 * x210 - x109 * x255 - x109 * x258 - x11 * x226 * x324 * x67 - x11 * x237 * x325 - x11 * x29 * x292 - x11 * x325 * x58 - x110 * x111 - x110 * x214 * x405 - x112 * x113 - x112 * x118 - x112 * x168 - x112 * x336 - x112 * x339 - x113 * x179 - 24 * x114 * x272 - 144 * x116 * x164 - x116 * x17 * x356 * x68 - x116 * x214 * x347 * x389 - x116 * x323 * x324 - x116 * x327 - 8 * x116 * x352 - x116 * x357 - x116 * x362 - 20 * x117 * x186 - x117 * x35 * xi - x118 * x179 - x119 * x124 * x340 - 96 * x119 * x164 - 56 * x119 * x169 - x119 * x323 * x324 - 8 * x119 * x352 - x119 * x357 - x119 * x362 - x12 * x145 - x12 * x252 * x63 - x12 - x120 * x121 - x120 * x222 - 4 * x120 * x319 - x120 * x48 - x121 * x47 - x121 * x49 - x122 * x123 - x124 * x137 - x124 * x175 * x67 - x124 * x184 * x325 * x50 - x124 * x32 * xi * xi_yy - x124 * x387 * x389 - x125 * x136 - x125 * x147 * x97 - x126 * x128 - x126 * x156 * xi - x126 * x353 - x126 * x57 - x126 * x62 - x127 * x128 - x127 * x353 - x127 * x57 - x127 * x62 - x128 * x185 - x128 * x267 - x128 * x59 - x129 * x130 * x139 + x129 * x17 * x23 * x275 * x344 + x129 * x17 * x23 * x61 + x129 * x17 * x24 * x275 * x343 + x129 * x17 * x24 * x60 + 68 * x129 * x2 * x24 * x63 + 8 * x129 * x2 * x64 * x75 - x129 * x282 * x291 - x129 * x53 - x129 * x59 * x77 - x130 * x134 - 18 * x130 * x85 - x131 * x133 - x131 * x167 * x18 * xi - x132 * x200 * x24 * xi - x134 * x367 - x136 * x301 - x137 * x85 - x138 * x139 - x138 * x147 * x17 - x140 * x141 - 24 * x140 * x170 - x140 * x188 - x140 * x189 - x140 * x224 - x140 * x273 - x140 * x369 - x140 * x372 - x140 * x378 - x140 * x396 - x141 * x142 - x142 * x188 - x142 * x189 - x142 * x224 - x142 * x273 - x142 * x369 - x142 * x372 - x142 * x378 - x142 * x396 + 16 * x143 * x2 * x63 * x65 - x143 * x223 - x143 * x374 - x144 * x146 - x147 * x202 * x60 - x147 * x203 * x58 - x147 * x219 * x25 - x148 * x17 * x256 * x405 - x149 * x151 - x149 * x316 * x68 - x149 * x359 - x149 * x379 - x149 * x384 - x15 * x4 - x150 * x154 * x71 * xi - x150 * x163 * x351 - x150 * x351 * x406 - x151 * x227 - x152 * x162 - x152 * x17 * x316 * x68 - x152 * x385 - x153 * x351 - x155 * x228 - x155 * x67 - x156 * x157 - x156 * x78 * xi - x158 * x276 - x158 * x28 - x158 * x315 - x158 * x81 - x159 * x160 - x159 * x161 - x160 * x358 - x161 * x358 - x162 * x163 - x163 * x382 - x163 * x385 - 156 * x164 * x229 - 102 * x164 * x235 - 54 * x164 * x252 - x165 * x201 - x165 * x218 - x166 * x34 - x166 * x351 - x166 * x383 - x168 * x179 - x169 * x65 * x99 - x17 * x186 * x26 + 6 * x17 * x22 * x275 * x343 * xi - x17 * x22 * x328 * x330 + 6 * x17 * x22 * x60 * xi + 9 * x17 * x23 * x275 * x343 + 9 * x17 * x23 * x60 - x17 * x239 * x387 + x17 * x275 * x344 * x63 * x70 - x17 * x316 * x35 - x17 * x52 * x63 * x65 + x17 * x61 * x63 * x70 - x170 * x182 - x170 * x192 - x171 * x356 - x172 * x174 - x172 * x257 - x173 * x178 * x85 - x174 * x400 - x175 * x176 - x175 * x230 - x176 * x358 * xi_y - x177 * x178 * x85 + 16 * x18 * x2 * x63 * x65 * xi_xy - x18 * x223 * xi_xy - x18 * x256 * x291 * x70 - x18 * x277 * x9 - x180 * x294 - x181 * x264 - x181 * x295 - x181 * x417 - x182 * x266 - x185 * x395 - x185 * x399 - x185 * x57 - x185 * x62 - x186 * x43 - x186 * x96 - x19 * x9 - x190 * x294 - x191 * x264 - x191 * x295 - x191 * x417 - x192 * x266 - x193 * x205 - x193 * x30 * x76 - x193 * x392 - x193 * x398 - x193 * x92 - x194 * x205 - x194 * x392 - x194 * x398 - x195 * x34 - x195 * x351 - x195 * x383 - x197 * x366 - x197 * x60 - x198 * x199 - x198 * x228 * x254 - x198 * x268 - x198 * x269 - x199 * x397 + 88 * x2 * x22 * x63 * xi + 36 * x2 * x23 * x63 + 24 * x2 * x64 * x65 + 32 * x2 * x64 * x70 * xi + 4 * x2 * x66 * x67 - x2 * x9 * xi_x * (x17 * (x2 * (Fx * (78 * B * x120 + 51 * B * x140 + 27 * B * x142 - Bp * x10 * x433 - 2 * Sp * x440 + 21 * x1 + 28 * x122 + 12 * x145 + x180 * x445 + x180 * x446 + x180 * x447 + x190 * x445 + x190 * x446 + x190 * x447 + x231 * x235 + x231 * x24 + x231 * x252 + x235 * x238 + x238 * x252 + 34 * x3 + x389 * x47 + x390 * x49 - x426 + 18 * x429 - 2) + 2 * u * (Bt * x441 * x6 + Gt * x433 * x444 + St * (4 * x4 + x425 + 4))) - x423 * (Fy * u * (x421 * x443 - x442) + 4 * Fy * x435 + Gh * (10 * x3 + 6 * x4 + x425 + 6) - x444 * (Bh * x3 + Bh * x4 + Bh - 2 * Sh))) + x18 * (u * (-6 * Gh * x16 * x433 + x2 * x6 * (u * (Fx * x421 * (12 * x424 + x443) - Fx * x442 + x436 * (Bt * x3 + Bt * x4 + Bt + St)) + x279 * x441) - x293 * x5) - x30 * (-Sp * x440 + x116 * x238 + x119 * x238 + 14 * x122 + 6 * x145 + x229 * x238 + x235 * x439 + x252 * x439 - x4 + x423 * x437 + x430 + x438 - 1)) - x329 * (x17 * (x427 + x434) * (x435 + x436 * (x434 + x5)) + x18 * (x10 * x432 + x22 * x428 + x22 * x431 + x437 * x6))) - x201 * x90 - x202 * x285 * x363 - x202 * x301 * x68 - 28 * x202 * x349 - x203 * x204 - x204 * x391 - x208 * x373 - x209 * x397 * x412 - x21 * x47 - x21 * x49 - x210 * x375 - x211 * x213 - x211 * x217 - x213 * x408 - x215 * x216 - x215 * x298 - x216 * x409 - x217 * x408 - x218 * x90 - 48 * x219 * x88 - 24 * x219 * x97 - 2 * x22 * x275 * xi_x ^ 2 * (x17 * (9 * x0 * x180 * x433 + x0 * x428 + x0 * x431 + x432 + x68 * (2 * x145 + 3 * x3 * x5 + x5 ^ 2)) + x37 * x422 * (4 * x3 + x425 + x427)) - x22 * x30 * x76 - x22 * x331 - x220 * x43 - x220 * x46 - x220 * x96 - x222 * x47 - x222 * x49 - x225 * x326 - x225 * x63 - x227 * x359 - x227 * x379 - x227 * x384 - 8 * x229 * x352 - x23 * x43 - x23 * x46 - x23 * x96 - x230 * x358 * xi_y - x231 * x232 - x231 * x233 - x232 * x238 - x232 * x239 - x233 * x238 - x234 * x236 - x234 * x237 - x235 * x361 * x386 - x235 * x394 - x236 * x240 - x236 * x281 - x237 * x240 - x237 * x281 - x24 * x327 - 12 * x241 * x271 * x349 - x241 * x391 * x412 * x85 - x242 * x243 - x242 * x251 - x242 * x410 - x242 * x411 - x242 * x414 - x243 * x244 - x244 * x251 - x244 * x349 * x86 - x244 * x410 - x244 * x411 - x244 * x414 - x245 * x246 - x245 * x259 - x246 * x247 - x247 * x259 - x247 * x322 * xi_y - x248 * x249 - x248 * x250 - x248 * x274 - x249 * x253 - x25 * x26 - x25 * x340 * xi - x250 * x253 - x252 * x394 - x253 * x274 - x254 * x413 - x255 * x375 - x257 * x400 - x258 * x375 - x260 * x285 * x49 - x260 * x290 * x47 * xi - x262 * x263 - x262 * x265 - x262 * x296 - x262 * x297 - x262 * x416 - x262 * x418 - x263 * x270 - x265 * x270 - 48 * x265 * x271 - x267 * x395 - x267 * x399 - x269 * x365 * x366 - x269 * x397 - x27 * x276 - x27 * x28 - x27 * x315 - x270 * x296 - x270 * x297 - x270 * x416 - x270 * x418 - x276 * x79 - x276 * x80 - x278 * x79 - x278 * x80 - x28 * x79 - x28 * x80 - x281 * x58 - 18 * x282 * x283 - x284 * x285 - x285 * x363 * x98 - x288 * x290 - x29 * x31 - x298 * x409 - x299 * x321 * x72 - x3 * x51 - x300 * x302 - x301 * x68 * x98 - x302 * x38 - x302 * x40 - x303 * x304 - x303 * x305 - x304 * x317 - x304 * x419 - x304 * x420 - x305 * x317 - x305 * x388 * x397 - x305 * x419 - x305 * x420 - x306 * x307 - x306 * x308 - x306 * x387 - x307 * x318 - x308 * x318 - x310 * x311 - x310 * x313 - x312 * x314 - x314 * x320 - x315 * x79 - x315 * x80 - x316 * x350 - x316 * x87 - 16 * x316 * x98 - x318 * x387 - 4 * x319 * x47 - 4 * x319 * x49 - x32 * x34 - x322 * x65 * xi_y - x331 * x90 - x331 * x93 - x332 * x334 - x332 * x338 - x333 * x335 - x333 * x376 - x333 * x377 - x334 * x364 - x335 * x337 - x337 * x376 - x337 * x377 - x338 * x364 - x34 * x91 - x340 * x97 - x341 * x342 - x341 * x370 - x341 * x371 - x342 * x345 - x345 * x370 - x345 * x371 - x348 * x41 - x349 * x350 - x349 * x87 - 20 * x349 * x98 - x351 * x393 - x351 * x68 * x70 - x351 * x91 - 20 * x351 * x93 - x354 * x355 * x71 - x359 * x360 - x360 * x384 - x365 * x368 - 24 * x365 * x380 - x365 * x415 - x368 * x85 - x38 * x39 - x380 * x381 - x382 * x406 - x383 * x393 - x387 * x388 * x85 - x389 * x413 - x39 * x40 - x390 * x391 * x71 - x395 * x407 - x399 * x407 - x40 * x65 * x77 - x401 * x402 - x401 * x403 - x41 * x42 - x415 * x85 - x47 * x48 - x48 * x49 - x55 * x56 - x57 * x59 - x59 * x62 - x69 * x72 - x69 * x74 - x72 * x83 - x79 * x81 - x80 * x81 - x84 * x87 - 42 * x88 * x89 - 68 * x89 * x97 - x92 * x93

    nothing
end
