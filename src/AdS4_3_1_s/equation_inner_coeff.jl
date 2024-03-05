
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

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")


  x0 = u ^ 3
x1 = 2 * x0
x2 = B * x0
x3 = exp(x2)
x4 = S0 * xi
x5 = S * x0
x6 = (S0 + S0_t * u + u * x4 + x5) ^ 4
x7 = x3 * x6
x8 = u ^ 2
x9 = x3 * x8
x10 = S * x8
x11 = 1 / u
x12 = S0 * x11
x13 = S0_t + x10 + x12 + x4
x14 = x13 ^ 4
x15 = 4 * x8
x16 = G * x0
x17 = cosh(x16)
x18 = S0 * xi_x
x19 = S0_x * xi
x20 = S0_x * x11
x21 = Fx * u
x22 = S0 ^ 2
x23 = 1 / x22
x24 = S0 * S0_tx
x25 = S0_t * S0_x
x26 = x24 - x25
x27 = x23 * x26
x28 = x21 + x27 + xi_x
x29 = Sp * x28 + St
x30 = x29 * x8
x31 = Sp * x8
x32 = S0 + x31 - 2 * x5
x33 = S0_tx + x18 + x19 + x20 - x28 * x32 + x30
x34 = sinh(x16)
x35 = Gp * x28 + Gt
x36 = x0 * x35
x37 = u ^ 4
x38 = G * x37
x39 = 4 * Gp * x0 - 12 * x38
x40 = 4 * S0
x41 = 4 * S0_xx
x42 = x28 ^ 2
x43 = Spp * x28
x44 = Spt + x43
x45 = 2 * x24 - 2 * x25
x46 = 2 * xi_x
x47 = 2 * x21 + x46
x48 = x23 * x45 + x47
x49 = S0_x - x1 * x29 + x44 * x8
x50 = S0 ^ 3
x51 = 1 / x50
x52 = 2 * S0_x * x51
x53 = Fx * x8
x54 = Fxp * u - x53
x55 = u * (Fxp * x28 + Fxt) + x23 * (S0 * S0_txx - S0_t * S0_xx) - x26 * x52 - x28 * x54 + xi_xx
x56 = 6 * x37
x57 = 4 * x0
x58 = S * x56 - Sp * x57 + Spp * x8
x59 = -x28 * x58 + x49
x60 = Bp * x28 + Bt
x61 = x0 * x60
x62 = 3 * x37
x63 = -B * x62 + Bp * x0
x64 = -x28 * x63 + x61
x65 = x13 ^ 2
x66 = Gpp * x28
x67 = Gpt + x66
x68 = 2 * Gp * x0 - 6 * x38
x69 = 2 * x27 + x47
x70 = 12 * u ^ 5
x71 = G * x70 - Gp * x56 + Gpp * x0
x72 = x0 * x67 - x28 * x71 - x35 * x62
x73 = Gp * x0 - 3 * x38
x74 = -x28 * x73 + x36
x75 = -B * x56 + 2 * Bp * x0
x76 = 2 * x75
x77 = Bpp * x28 + Bpt
x78 = B * x70 - Bp * x56 + Bpp * x0
x79 = Fy * u
x80 = S0 * S0_ty
x81 = S0_t * S0_y
x82 = x80 - x81
x83 = x23 * x82
x84 = x79 + x83 + xi_y
x85 = Sh + Sp * x84
x86 = S0 * xi_y + S0_ty + S0_y * x11 + S0_y * xi - x32 * x84 + x8 * x85
x87 = Gh + Gp * x84
x88 = x0 * x87
x89 = x34 * x86
x90 = 4 * S0_yy
x91 = x84 ^ 2
x92 = Sph + Spp * x84
x93 = 2 * xi_y
x94 = 2 * x79 + x93
x95 = x23 * (2 * x80 - 2 * x81) + x94
x96 = x8 * x92
x97 = S0_y - x1 * x85 + x96
x98 = 4 * x84
x99 = 2 * S0_y
x100 = -Fy * x8 + Fyp * u
x101 = u * (Fyh + Fyp * x84) - x100 * x84 + x23 * (S0 * S0_tyy - S0_t * S0_yy) - x51 * x82 * x99 + xi_yy
x102 = Bh + Bp * x84
x103 = x0 * x102
x104 = x103 - x63 * x84
x105 = Gph + Gpp * x84
x106 = 2 * x0 * x105 - x56 * x87
x107 = 2 * x83 + x94
x108 = x73 * x84
x109 = -x108 + x88
x110 = Bph + Bpp * x84
x111 = 2 * x73
x112 = 2 * x8
x113 = -16 * S * x0 + 8 * S0 + 8 * x31
x114 = 2 * xi_xy
x115 = 2 * S0_xy
x116 = S0_ty * S0_x
x117 = S0_tx * S0_y
x118 = S0 * S0_txy - S0_t * S0_xy
x119 = -S0_y * x45 * x51 + u * (Fxh + Fxp * x84) + u * (Fyp * x28 + Fyt) - x100 * x28 + x114 + x23 * (-x116 + x117 + x118) + x23 * (x116 - x117 + x118) - x52 * x82 - x54 * x84
x120 = 8 * S0_t * x50
x121 = S0 ^ 4
x122 = x51 / 8
ABCS[1] = x1 * x7
ABCS[2] = 8 * x6 * x9
ABCS[3] = 4 * u * x7
if u == 0.0 
	ABCS[4] = -4*S0^3*Sp
else
	ABCS[4] = x13 * x8 * (x17 * (-4 * S0_txx - 8 * S0_x * xi_x - x11 * x41 - x15 * (Sb - Spp * x42 + x44 * x48) + 4 * x28 * x49 + 4 * x28 * x59 + 4 * x32 * x55 + 4 * x33 * x64 - x40 * xi_xx - x41 * xi) + x33 * x34 * (x28 * x39 - 4 * x36)) + x14 * x15 * x3 + x15 * x17 * x33 ^ 2 + x65 * x8 * (x17 * (2 * x0 * (Bb - Bpp * x42 + x48 * x77) - x28 * (2 * x0 * x77 - x56 * x60) + x54 ^ 2 - x55 * x75 - 2 * x64 ^ 2 - x69 * (x0 * x77 - x28 * x78 - x60 * x62) - 2 * x74 ^ 2) + x34 * (-x1 * (Gb - Gpp * x42 + x48 * x67) + x28 * (2 * x0 * x67 - x35 * x56) + x55 * x68 + x69 * x72 + x74 * (-x28 * x76 + 4 * x61))) + x8 * (x13 * (-x17 * (4 * S0_tyy + 8 * S0_y * xi_y - 4 * x101 * x32 + 4 * x104 * x86 + x11 * x90 + x15 * (Spp * x91 + Ss + x92 * x95) + x40 * xi_yy + x90 * xi - x97 * x98 - x98 * (-x58 * x84 + x97)) + x89 * (x39 * x84 - 4 * x88)) + 4 * x17 * x86 ^ 2 + x65 * (x17 * (-x1 * (Bpp * x91 + Bs + x110 * x95) + x100 ^ 2 + x101 * x75 - 2 * x104 ^ 2 + x107 * (x0 * x110 - x102 * x62 - x78 * x84) - 2 * x109 ^ 2 + x84 * (2 * x0 * x110 - x102 * x56)) + x34 * (-x1 * (Gpp * x91 + Gs + x105 * x95) + x101 * x68 + x106 * x84 + x107 * (x0 * x105 - x62 * x87 - x71 * x84) - x109 * (4 * x103 - x76 * x84)))) * exp(2 * x2) + x9 * (x14 * (Bd * x112 * x17 ^ 2 * x63 + Gd * x111 * x8) + x65 * (-x113 * (S0 / (2 * x8) + Sd * u + x11 * x122 * (x120 + 8 * x121 * xi) + x122 * (4 * S0_t ^ 2 * x22 + 4 * S0_x ^ 2 - S0_xx * x40 + 4 * S0_y ^ 2 - S0_yy * x40 + x120 * xi + 4 * x121 * xi ^ 2)) + x17 * (4 * x0 * (Gc + x105 * x28 + x66 * x84 + x67 * x84) + 2 * x104 * x74 - 2 * x106 * x28 - 2 * x107 * x72 - 2 * x109 * x64 - x111 * x119) + x34 * (x100 * (-2 * Fxp * u + 2 * x53) + (-2 * x108 + 2 * x88) * (-x28 * x68 + 2 * x36))) - x89 * (8 * S0_tx - x113 * x28 + 8 * x18 + 8 * x19 + 8 * x20 + 8 * x30) + (x17 * (x109 * x33 + x74 * x86) + x34 * (S0 * x114 + 2 * S0_txy + S0_x * x93 + S0_y * x46 - x107 * x59 + x11 * x115 + x112 * (Sc + x28 * x92 + x43 * x84 + x44 * x84) + x115 * xi - x119 * x32 - x28 * (-x57 * x85 + 2 * x96 + x99))) * (4 * S0_t + 4 * x10 + 4 * x12 + 4 * x4))
end
    nothing
end
