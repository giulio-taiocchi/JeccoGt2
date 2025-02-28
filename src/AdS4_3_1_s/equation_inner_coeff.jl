
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


x0 = u ^ 2
x1 = Gp ^ 2
x2 = u ^ 4
x3 = G * Gp
x4 = G ^ 2
x5 = B * u
x6 = (Bp - 3 * x5) ^ 2
x7 = G * u ^ 3
ABCS[1] = 4 * x0
ABCS[2] = 24 * u
ABCS[3] = 9 * u ^ 6 * x4 - 6 * u ^ 5 * x3 + x1 * x2 + x2 * x6 * cosh(x7) ^ 2 + 24
ABCS[4] = u * (S0 * u * xi + S0 + S0_t * u) * (9 * B ^ 2 * x0 + Bp ^ 2 - 6 * Bp * x5 - 12 * u * x3 + 18 * x0 * x4 + 2 * x1 + x6 * cosh(2 * x7)) / 2

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
x4 = S0_t * u
x5 = u * xi
x6 = S0 * x5
x7 = S * x0
x8 = S0 + x6 + x7
x9 = x4 + x8
x10 = x9 ^ 2
x11 = 2 * x10
x12 = x11 * x3
x13 = 3 * x1
x14 = Bp * x3
x15 = -x14
x16 = -Bp
x17 = 3 * u
x18 = B * x17
x19 = x16 + x18
x20 = G * x0
x21 = 2 * x20
x22 = cosh(x21)
x23 = x19 * x22
x24 = x13 + x15 + x23 * x3
x25 = u * x2
x26 = x10 * x25
x27 = 2 * Gp
x28 = G * u
x29 = 6 * x28
x30 = -x19
x31 = sinh(x21)
x32 = x30 * x31
x33 = -x27 + x29 - x32
x34 = x0 * x10
x35 = 4 * S0
x36 = Spp * x35
x37 = S0 ^ 2
x38 = Bpp * x37
x39 = Sp * u
x40 = 24 * S0
x41 = x37 * xi
x42 = u * x37
x43 = S * x3
x44 = 4 * Spp
x45 = S ^ 2
x46 = u ^ 5
x47 = x45 * x46
x48 = Sp ^ 2 * x0
x49 = u ^ 6
x50 = 36 * S0
x51 = B * x46
x52 = S * x51
x53 = u ^ 7
x54 = S * x53
x55 = 6 * B
x56 = x54 * x55
x57 = u ^ 4
x58 = Sp * x57
x59 = x55 * x58
x60 = S * x57
x61 = Bp * S0
x62 = 2 * Sp
x63 = S * x49
x64 = Bp * x63
x65 = x0 * x61
x66 = 2 * S0
x67 = Bpp * x7
x68 = 2 * x5
x69 = S0 * xi
x70 = Sp * x3
x71 = u ^ 8
x72 = B * x71
x73 = 9 * x37
x74 = B * x3
x75 = Bp * x45
x76 = xi ^ 2
x77 = 4 * x42
x78 = B * x54
x79 = B * x63
x80 = x55 * x69
x81 = x61 * xi
x82 = x46 * x81
x83 = x57 * x61
x84 = x83 * xi
x85 = Bpp * x60
x86 = x66 * xi
x87 = G * Gp
x88 = Bp ^ 2
x89 = u ^ 9
x90 = x45 * x89
x91 = x0 * x37
x92 = u ^ 10
x93 = x75 * x92
x94 = Bp * x55
x95 = x37 * x57
x96 = 24 * x41
x97 = B ^ 2
x98 = 18 * S0
x99 = S * x71
x100 = 12 * x41
x101 = x63 * x66
x102 = 12 * G
x103 = Gp * x102
x104 = x45 * x92
x105 = G ^ 2
x106 = Gp ^ 2
x107 = S * x72
x108 = G * S
x109 = Gp * x71
x110 = x108 * x109
x111 = u ^ 11 * x45
x112 = 9 * x97
x113 = x46 * x73
x114 = 18 * x105
x115 = x37 * x46
x116 = x106 * x90
x117 = x106 * x91
x118 = Bp * x46
x119 = 18 * x97
x120 = S * x69 * x89
x121 = x54 * x88
x122 = Gp * x46
x123 = G * x122
x124 = 36 * x105
x125 = x106 * x54
x126 = x35 * xi
x127 = 15 * x57
x128 = x37 * x76
x129 = B * x128
x130 = x41 * x49
x131 = Bp * x91
x132 = x57 * x88
x133 = x106 * x57
x134 = 4 * x41
x135 = x46 * x88
x136 = x128 * x49
x137 = x53 * x73 * x76
x138 = x106 * x115 * x76
x139 = x108 * x72
x140 = 12 * S0
x141 = Gp * x31
x142 = x31 * x61
x143 = x102 * x142
x144 = 18 * G
x145 = x144 * x31
x146 = B * x111
x147 = x141 * x55
x148 = 6 * G
x149 = x31 * x93
x150 = Bp * x31
x151 = x150 * x95
x152 = x150 * x90
x153 = x27 * x31
x154 = B * x69
x155 = x108 * x154 * x31 * x89
x156 = x31 * xi
x157 = G * x31
x158 = B * x130
x159 = x31 * x41
x160 = x118 * x159
x161 = x41 * x57
x162 = x136 * x150
x163 = x118 * x128
x164 = Bpp * u
x165 = 15 * x1
x166 = 7 * x14
x167 = x102 * x122
x168 = 9 * x49
x169 = x168 * x97
x170 = x114 * x49
x171 = 2 * x133
x172 = -Gp
x173 = G * x17
x174 = x172 + x173
x175 = -x174
x176 = x175 * x32 * x57
x177 = 2 * x176
x178 = x3 * x88
x179 = 6 * x1
x180 = x13 + 5
x181 = Bpp + u * (-Bp * (x179 + 7) + x178 + x18 * x180)
x182 = u * x22
x183 = S0_t ^ 2 * u
x184 = Bpp * x8
x185 = x178 * x8
x186 = 7 * x5
x187 = x5 + 1
x188 = x179 * x187
x189 = x186 + x188
x190 = -x62
x191 = S * u
x192 = 5 * x5
x193 = x13 * x187 + x192
x194 = S * x17
x195 = x22 * x8
x196 = 8 * x70
x197 = Spp * u
x198 = 16 * x7
x199 = S0 * x132
x200 = -Sp
x201 = 9 * x191
x202 = x55 * x60
x203 = x13 + 7
x204 = 9 * S
x205 = 6 * S0
x206 = x204 * x89
x207 = S0 * x169
x208 = S * x114
x209 = x208 * x89
x210 = Bpp * S0
x211 = x3 * xi
x212 = x55 * x61
x213 = x212 * x46
x214 = x106 * x46
x215 = x53 * x69
x216 = x49 * xi
x217 = x140 * x87
x218 = x35 * x5 + x35
x219 = 2 * S0_t
x220 = 2 * Gpp
x221 = S0 * x220
x222 = x22 * x30
x223 = 2 * x0
x224 = x175 * x223
x225 = x222 * x224
x226 = -x13 + x14
x227 = 2 * u
x228 = Bpp - u * (Bp * (7 - x179) + x178 + x18 * (x13 - 5))
x229 = 22 * Gp * x60
x230 = x102 * x58
x231 = Gp * S0
x232 = 10 * u * x231
x233 = x220 * x7
x234 = 4 * x0
x235 = Gp * Sp * x234
x236 = G * x98
x237 = x236 * x3
x238 = 54 * x108 * x46
x239 = x31 * x88
x240 = x239 * x63
x241 = S0 * x0 * x239
x242 = 14 * x211 * x231
x243 = x31 * x71
x244 = x204 * x243 * x97
x245 = S0 * x31
x246 = x112 * x245 * x46
x247 = x220 * x6
x248 = 30 * x20 * x69
x249 = x156 * x199
x250 = x156 * x207
x251 = x31 * x54 * x94
x252 = x212 * x31 * x57
x253 = x156 * x213
x254 = x148 * x54
x255 = x31 * x69
x256 = Bp * x254 + Bpp * x31 * x6 + Gp * x56 - 5 * u * x142 + x0 * x150 * x62 + x122 * x80 - 18 * x139 - x144 * x154 * x49 + x148 * x82 + x148 * x83 - 11 * x150 * x60 + x165 * x255 - x166 * x255 + x195 * x224 * x30 + x210 * x31 + x231 * x55 * x57 - x236 * x51 + 9 * x245 * x74 - x27 * x64 - x27 * x65 - x27 * x84 + 27 * x31 * x52 - x31 * x59 + x31 * x67
x257 = u * x9
x258 = 6 * u
x259 = x10 * x258
x260 = 2 * x191 + x200
x261 = -x260
x262 = x261 * x3
x263 = S0 + x262
x264 = 4 * x263
x265 = 4 * x2
x266 = x265 * x9
x267 = x25 * x9
x268 = x10 * x31
x269 = x17 * x268
x270 = B_y * x34
x271 = 2 * x174
x272 = cosh(x20) ^ 2
x273 = x11 * x272
x274 = Bp_x * x2
x275 = x172 + 2 * x28
x276 = 12 * u
x277 = x10 * (x220 + x275 * x276)
x278 = 6 * x26
x279 = x266 * (Spp - 4 * x39 + 6 * x43)
x280 = 2 * x34
x281 = G_y * x280
x282 = x34 * xi_y
x283 = x30 ^ 2
x284 = x283 * x31
x285 = x271 * x282
x286 = B * x227 + x16
x287 = Bpp + x258 * x286
x288 = x268 * x287
x289 = x19 * x272
x290 = S0_x * x266 * x289
x291 = 4 * x289
x292 = S0_tx * x267
x293 = x33 * x9
x294 = x227 * x293
x295 = x267 * xi_x
x296 = S0_y * x33
x297 = x272 * x9
x298 = x19 * x297
x299 = x234 * x298
x300 = 1 / x37
x301 = x300 * (S0 * S0_ty - S0_t * S0_y)
x302 = x2 * x280
x303 = x2 * xi_x
x304 = x280 * x283
x305 = 2 * x7
x306 = x5 * x66 + x66
x307 = S0 * S0_tx
x308 = x300 * (-S0_t * S0_x + x307)
x309 = x273 * x287
x310 = x284 * x34
x311 = x271 * x34
x312 = x301 * x311
x313 = -Gp + x173 - x32
x314 = G_x * x302
x315 = x2 * x308
x316 = -S0 + x260 * x3
x317 = x295 * x316
x318 = x303 * x311
x319 = x267 * x308 * x316
x320 = x311 * x315
x321 = u * x263
x322 = x261 * x321
x323 = x261 * x263
x324 = x261 * x4 * (x262 + x66) + x37 * (Sp - x194 + x261 * x5)
x325 = x27 - x29 - x32
x326 = x2 * x325
x327 = Spp * x66
x328 = 24 * x43
x329 = Gp * x148
x330 = 9 * x105
x331 = x142 * x27
x332 = 3 * x141
x333 = B * x332
x334 = 3 * G
x335 = 6 * x108
x336 = x105 * x168
x337 = u * x272
x338 = x186 - x188
x339 = -3 * B * x0 * x187 + x192
x340 = x11 * x2
x341 = x9 * (x328 - 16 * x39 + x44)
x342 = B_x * x2
x343 = x272 * xi_y
x344 = x19 * x257
x345 = x343 * x344
x346 = 4 * S0_y
x347 = x340 * (Gpp + x258 * x275)
x348 = x10 * (2 * Bpp + x276 * x286)
x349 = x174 + x32
x350 = x326 * x9
x351 = S0_x * x350
x352 = 2 * x325
x353 = x272 * x301
AA[1,1] = x12 * x2
AA[1,2] = 0
BB[1,1] = x26 * (x24 + 8)
BB[1,2] = -x33 * x34
CC[1,1] = x25 * (-B * x100 * x118 + 4 * Bp * x141 * x161 - 5 * Bp * x42 + Bpp * x45 * x49 - Gp * x107 * x140 * x156 + Gp * x126 * x150 * x54 - 18 * S * x82 + 48 * S0 * x43 - S0 * x59 - Sp * x46 * x80 - Sp * x56 + x1 * x96 - x100 * x14 - x100 * x141 * x51 + x101 * x88 - x102 * x160 - x103 * x104 - x103 * x136 - x103 * x95 - x104 * x147 + x105 * x50 * x99 + x106 * x35 * x63 - 12 * x107 * x81 - 24 * x110 * x69 + x111 * x112 + x111 * x114 + x113 * x97 + x114 * x115 + x114 * x128 * x53 + 2 * x116 + 2 * x117 + x119 * x120 + x119 * x130 + x120 * x124 + x121 * x86 - x123 * x96 + x124 * x130 + x125 * x126 + x127 * x129 + x128 * x135 + x129 * x145 * x53 + x131 * x153 - 7 * x131 * x76 + 2 * x132 * x41 + x133 * x134 - x136 * x147 - x136 * x94 + x137 * x97 + 2 * x138 + x139 * x31 * x50 - x140 * x141 * x78 + x141 * x35 * x64 - x143 * x54 - x143 * x99 * xi + x145 * x146 + x145 * x37 * x51 - x147 * x95 - x148 * x149 - x148 * x151 - x148 * x162 + x152 * x27 + x153 * x163 + 36 * x155 + 36 * x157 * x158 + x183 * (-x118 * x55 + x132 + x164 + x165 - x166 - x167 + x169 + x170 + x171 + x177 + x181 * x182 + 4) + x195 * (u * (-Bp * (S0 * (x189 + 5) + x3 * (x190 + x191 * (x179 + 11))) + x18 * (S0 * (x193 + 3) + x3 * (x190 + x194 * (x1 + 3))) + x185) + x184) + x219 * (Bp * x58 + S0 * x164 + S0 * x170 - 3 * Sp * x51 + x1 * x140 - x103 * x99 + x112 * x215 + x114 * x215 - x118 * x204 + x121 - x123 * x140 + 2 * x125 + x127 * x154 + x133 * x66 + x135 * x69 - x14 * x205 + x177 * x8 + x182 * (u * (-Bp * (S0 * (x189 + 6) + x3 * (x200 + x201 + x202)) + x18 * (S0 * (x193 + 4) + x3 * (x191 * x203 + x200)) + x185) + x184) - x196 + 2 * x197 + x198 + x199 + x206 * x97 + x207 + x209 + x210 * x211 - x212 * x216 - x213 + x214 * x86 - x216 * x217 + x218 - 7 * x65 * xi + 21 * x79 + x85 - x94 * x99) + x3 * x38 * x76 + x36 * x5 + x36 + x38 * x68 + x38 - x39 * x40 - x40 * x54 * x87 + 8 * x41 + x44 * x7 + 27 * x45 * x72 + 12 * x47 - 4 * x48 + x50 * x52 - 11 * x53 * x75 - x55 * x93 - 16 * x60 * x61 - 12 * x61 * x78 + x62 * x64 + x62 * x65 + x62 * x84 + x66 * x67 + 32 * x69 * x7 - 16 * x69 * x70 + 42 * x69 * x79 + x73 * x74 + x76 * x77 + x85 * x86 + x88 * x90 + x88 * x91 - x94 * x95 + x97 * x98 * x99)
CC[1,2] = -x257 * (x221 - x229 - x230 - x232 + x233 + x235 + x237 + x238 - x240 - x241 - x242 - x244 - x246 + x247 + x248 - x249 - x250 + x251 + x252 + x253 + x256 + x4 * (x220 + x225 + x227 * (-Gp * (x226 + 7) + x173 * (x226 + 5)) + x228 * x31))
SS[1] = B_x * x272 * x278 + B_x * x289 * x302 - B_y * x269 + Bp_y * x268 - G_y * x259 + Gp_y * x11 - S0_ty * x294 + S_x * x2 * x299 + S_x * x25 * x264 + 8 * S_x * x267 - S_y * x223 * x293 - Sp_x * x266 + x19 * x285 + x19 * x312 + x222 * x281 + x222 * x285 - x23 * x312 - x257 * x33 * x66 * xi_y + x263 * x294 * x301 - x265 * x300 * (-S0_x * x324 + x307 * x322 + x323 * x42 * xi_x) + x270 * x271 - x270 * x32 + x272 * x283 * x302 * x308 + x272 * x303 * x304 - x273 * x274 - x277 * x301 - x277 * xi_y + x279 * x308 + x279 * xi_x + x282 * x284 - x288 * x301 - x288 * xi_y + x289 * x295 * x35 + x290 * x5 + x290 + x291 * x292 + x291 * x317 + x291 * x319 - x294 * x316 * xi_y - x296 * x68 * x9 - x296 * (x305 + x306 + 2 * x4) + x301 * x310 + x303 * x309 + x309 * x315 + x313 * x314 + x313 * x318 + x313 * x320
AA[2,1] = 0
AA[2,2] = x12
BB[2,1] = x326 * x34
BB[2,2] = -u * x10 * (x24 - 8)
CC[2,1] = x267 * (-x221 + x229 + x230 + x232 - x233 - x235 - x237 - x238 + x240 + x241 + x242 + x244 + x246 - x247 - x248 + x249 + x250 - x251 - x252 - x253 + x256 - x4 * (-x181 * x31 + x220 - x225 + x227 * (Gp * (-x15 - x203) + x173 * (x15 + x180))))
CC[2,2] = x227 * (-B * x137 * x157 - Gp * x152 - S0 * Sp * x276 + S0 * x208 * x71 + S0 * x328 + Spp * x305 + x101 * x106 - x104 * x329 + x104 * x333 + x105 * x113 + x105 * x137 - x107 * x236 * x31 - x110 * x140 * xi + x111 * x330 + x116 + x117 + x122 * x159 * x55 + x125 * x86 + x128 * x227 + x129 * x332 * x49 - x131 * x141 + x134 - x136 * x329 + x138 - x141 * x163 + x142 * x254 - x145 * x158 - 9 * x146 * x157 + x147 * x69 * x99 + x148 * x160 + x149 * x334 - x150 * x161 * x27 + x151 * x334 - 18 * x155 - x157 * x51 * x73 + x162 * x334 - x167 * x41 + x170 * x41 + x171 * x41 + x183 * (-6 * x123 + x133 - x176 - x228 * x337 + x336 + 2) - x196 * x69 + x198 * x69 + x209 * x69 - x217 * x54 + x219 * (S0 * x133 + S0 * x336 + x105 * x206 - x109 * x335 - x123 * x205 + x125 - x176 * x8 + x197 - x205 * x216 * x87 + x214 * x69 + x215 * x330 + x306 + x337 * (u * (-Bp * (S0 * (-x338 - 6) + x3 * (Sp - x201 + x202)) + x18 * (S0 * (-x339 - 4) + x3 * (Sp + x191 * (x13 - 7))) + x185) - x184) + 8 * x7 - 4 * x70) + x231 * x31 * x56 + x243 * x335 * x81 + x272 * x8 * (u * (S0 * x18 * (-x339 - 3) + x13 * (x194 * (x1 - 3) + x62) + x14 * (x190 + x191 * (11 - x179)) + x185 + x61 * (x338 + 5)) - x184) + x327 * x5 + x327 - x329 * x95 - x331 * x54 * xi - x331 * x63 + x333 * x95 + 6 * x47 - 2 * x48)
SS[2] = -B_y * x259 * x272 + Bp_y * x273 - G_x * x278 + Gp_x * x340 - S0_ty * x257 * x291 + S0_y * x297 * (-B * x276 + 4 * Bp) + S_x * x223 * x350 + 8 * S_y * x257 - S_y * x299 + 4 * S_y * x321 - Sp_y * (x218 + 4 * x4 + 4 * x7) + x19 * x31 * x34 * x342 - x19 * x318 - x19 * x320 + x23 * x314 + x23 * x318 + x23 * x320 + x264 * x344 * x353 - x268 * x274 + x269 * x342 + 2 * x270 * x289 + x281 * x349 + x285 * x349 + x288 * x303 + x288 * x315 + x292 * x352 + x295 * x325 * x66 - x298 * x346 * x5 - x300 * (S0_ty * x322 * x35 + x323 * x77 * xi_y - x324 * x346) + x301 * x341 + x303 * x310 + x304 * x343 + x304 * x353 - x308 * x347 + x310 * x315 - x311 * x342 + x312 * x349 - 4 * x316 * x345 + x317 * x352 + x319 * x352 + x341 * xi_y - x343 * x348 - x345 * x35 - x347 * xi_x - x348 * x353 + x351 * x68 + 2 * x351

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


if u==0.0
	#println("entered Sd u zero")
	ABCS[4] = 4*S0^3*Sp
else
	ABCS[4] = x6 * (S0 * x19 * x5 * x9 + S0_t * x17 + u * x37 * x50 * (-x228 * x61 * (B * x106 * x255 + Bh * Gh * x156 + Fy * x173 * x220 + Fyh * x121 - Fyh * x243 + Fyp * x161 + Gh * x124 - Gs * x32 - S0_t * S0_yy * x240 + S0_ty * (S0 * S0_y * (-x101 * x248 * x260 + x182) + x104 * (6 * Bh * x52 + Fyp * x83 + x131 * (x7 - 1) + x183 * x248 * x260 + x64 * (Gp + x129 * (x261 - 5)))) + S0_tyy * x167 + S0_yy * x242 + x103 * x239 + x106 * x120 * x2 - x106 * x247 + x120 * xi_yy + x126 * x229 * x260 + x131 * x134 - x131 * x148 * x220 - x139 * x48 * x71 + x140 * xi_y + x144 * x244 + x149 * x244 - 30 * x160 * xi_y - x164 * x181 + x165 * xi_y + x179 * x181 - x18 * x220 * x245 + x200 * x245 + x203 * x245 - x210 * x220 * x253 - x220 * x252 * x256 * xi_y - x231 * x250 + x231 * x258 - x235 * x249 + x235 * x257 - x236 * x238 - x237 * xi_yy - 7 * x238 * x86 - x241 * x27 + x251 * x253 * xi_y) + x53 * (x109 * x227 * x228 + x119 ^ 2 * x92 + x119 * x2 * x223 * x66 + x122 * x66 ^ 2 + x223 * (u * x224 * x39 + x225 * x85 - x30 * (Fyph - u * (Fyh + x224 * x64) + x183 * x225)) - x228 * (12 * B * x229 + x128 * x39 * (-x232 * x85 + x234 * x30) + x230 * x231 - 6 * x234 * x30 * x85 + x32 * (Bs + x230 * x235 + x233 * x236 + 6 * x86 * (Bh + x143 * x228))) + x87 ^ 2 * x92)) + x1 * x17 - x11 * x186 * x53 * x81 * (3 * B * Fx * Gh * x32 * x91 + 3 * B * Gh * x18 * x32 * xi_x + 3 * B * Gt * S0_t * S0_y * x18 * x30 + 3 * Bh * G * S0_t * S0_x * x18 * x30 - Bh * Gt * x2 * x32 + 3 * Bt * Fy * G * x32 * x91 + 3 * Bt * G * S0_ty * x18 * x33 + 3 * Bt * G * x18 * x32 * xi_y + Bt * Gh * x2 * x32 - Bt * x172 * x18 * x71 + 36 * Fx * Fy * G * x18 * x32 - Fx * Fyp * x158 + 30 * Fx * G * S0_ty * x2 * x33 + 30 * Fx * G * x2 * x32 * xi_y + 7 * Fx * Gh * x32 * x6 + Fx * Gp * S0_t * S0_y * x30 * x6 + 3 * Fxh * G * x32 * x6 - Fxh * x121 + 3 * Fxp * G * S0_t * S0_y * x30 * x6 - Fxp * x161 + 30 * Fy * G * x2 * x32 * xi_x + Fy * Gp * S0_t * S0_x * x30 * x6 + 7 * Fy * Gt * x32 * x6 - Fy * x135 * x156 + 3 * Fyp * G * S0_t * S0_x * x30 * x6 + 3 * Fyt * G * x32 * x6 - Fyt * x121 + 12 * G * S0 * S0_t * S0_x * S0_y * u + 24 * G * S0_ty * x33 * x6 * xi_x + 24 * G * S0_x * S0_y * x31 * x6 + 6 * G * u * x32 * xi_xy + 24 * G * x32 * x6 * xi_x * xi_y + 2 * Gc * x32 + 6 * Gh * u * x32 * xi_x - Gh * x123 + 2 * Gp * S0_t * S0_xy * x30 + 2 * Gp * S0_ty * S0_x * x30 - Gp * x125 * x20 + 6 * Gt * S0_ty * u * x33 + 6 * Gt * u * x32 * xi_y - Gt * x124 - 3 * Gt * x147 * x148 + S0 * S0_tx * (S0_y * (-x101 * x141 - x182) + u * x30 * (-u * (Bh * x185 + x139 + x164 - x184 * x64) + x141 * x183 + x169 * (x7 + 2)) + x142) - 2 * x120 * xi_xy - x126 * x127 * x30 - x129 * x130 * x30 - x131 * x132 - x133 * x134 - x135 * x136 - x136 * x162 - x138 * xi_y - x140 * xi_x - x142 * x76 - x144 * x146 - x146 * x149 - x151 * x152 - x152 * x154 - x163 * xi_y - x165 * xi_x - x167 * x168 - x169 * x171 - x173 * x175 - x173 * x177 - x178 * x179 - x180 * x181) + x14 * x25 * x9 / x2 + 4 * x15 * x16 * x186 * (-x2 * x53 * (x49 * x90 + x60 * x87) + x61 * (Fx * Fy * Sp * x155 * x210 - Fx * x211 * x212 + Fxh * x189 + Fxh * x191 - Fxh * x197 + Fyt * x189 + Fyt * x191 - Fyt * x197 - 4 * S * x153 * xi_xy - 8 * S0 * x125 * x41 + S0_ty * x190 * x32 - S0_ty * x198 * x219 - S0_ty * x206 * x216 - S0_x * S0_y * x206 * x31 + S0_xy * x214 * x215 - S0_xy * x219 * x34 - Sc * x156 - Sh * x153 * x202 + Sp * x151 * xi_y + Sp * x156 * xi_xy + Sp * x200 * xi_x + x107 * x168 * x23 - x107 * x174 - x107 * x176 + x108 * x32 * x89 + 4 * x125 * x2 * x207 - 2 * x127 * x33 + 4 * x130 * x215 - x147 * x201 + x147 * x204 - x147 * x209 * x211 - x151 * x195 - 16 * x159 * x208 * x209 - x159 * x211 * x213 + x170 * x206 * xi_y + x175 * x221 - x176 * x218 * x91 + x177 * x221 - x187 * x42 - x187 + x188 * x190 * xi_y + x192 * x193 + 2 * x193 * x194 + x195 * x217 + x198 * x32 * x86 - x199 * x200 - x199 * x203 + x2 * x56 * (S0_y * (S0_t * x222 + x20 * x22 - 2 * x207) - x222 * x39 + x30 * (Fy * (S0 + x6 * (-x21 - 14 * x22)) - x222 * xi_y - 4 * x45)) + x201 * x220 - x204 * x220 - x205 * x212 * xi_x + x205 * x213 * x220 + x206 * x39 * x76)) - 12 * x16 ^ 4 * x8 + x16 * (4 * x49 * x50 * x6 * x61 * x87 - x54 * (-S0 * xi_yy - S0_tyy - S0_y * x47 - S0_yy * x15 - S0_yy * xi + x109 * x24 * x35 + x119 * x49 * x50 * x6 + x50 * (x30 * x74 + x72) * (S0_y * x113 + x112 * x40 + x48 * (Sph + x114 * x47 - x81 * (Sh + x110 * x64))) - x6 * (-Spp * x75 ^ 2 + Ss + (Sph + Spp * x75) * (x47 + 2 * x73 + x94)) + x75 * (S0_y * x102 + x100 * x95 + x6 * (Sph + u * (Fy * Spp - Fy * x96 - 2 * Sh) + x98 * xi_y)))) + x19 * x20 * x25 * x8 + 4 * x36 * x38 * (x26 + x27 - x28 - x29 + x30 * x31 + x32 * xi ^ 2 + x33 * x34 * xi) + x38 * x61 * x81 * (Fxph + Fxpp * x75 + Fypp * x80 + Fypt + u * x66 * x69 - u * (Fxh + Fxp * x75) - u * (Fyp * x80 + Fyt) - x75 * (Fxpp + x69 * x81) - x80 * (Fypp + x65 * x81) - x87 * x90 * x93) - x49 ^ 2 * x55 + x49 * x51 * x60 * x62 * x9 + (x16 * (4 * x262 * x61 * x90 + x54 * (S0 * xi_xx + S0_txx + S0_x * x59 + S0_xx * x15 + S0_xx * xi + x262 * x266 - x36 * (Fxt * x104 - S0_t * x28 + S0_txx * x30 + x107 * x265 + x216 * x89 + x263 - x264 * x76 + x33 * xi_xx + x56 * (-x198 + x264)) - x50 * (x30 * x79 + x77) * (S0_x * x113 + x112 * x57 + x48 * (Spt + x114 * x59 - x81 * (St + x110 * x68))) + x6 * (Sb - Spp * x80 ^ 2 + (Spp * x80 + Spt) * (x59 + 2 * x68 + 2 * x78)) - x80 * (S0_tx * x100 * x6 + S0_x * x102 + x6 * (Spt + u * (Fx * Spp - Fx * x96 - 2 * St) + x98 * xi_x)))) + x37 * (u * x53 * (Bb * x228 - Bp * Fxt * x284 - Bp * x265 * x92 - Fxp * u * x59 + Fxp * x0 * x198 * x70 - 2 * Fxpt + 6 * Fxt * x148 + Fxt * x81 - G * x290 * x297 + Gt ^ 2 * x92 + 12 * S0_t * x299 * x7 + S0_tx * x35 * x81 * (u * x198 * (-S0_t * (x230 + 9 * x292 * x91) + x227) + x30 * (u * (Fx * (-x285 + x301 + 27 * x7 + 2) + x128 * (Bt + Gt * x52)) + x67 + xi_x * (x301 + 24 * x7))) + S0_txx * x226 * x228 * x99 + S0_xx * x285 * x34 * x70 + u ^ 9 * x265 * x295 + u * x69 ^ 2 + 6 * x18 * x281 * x70 * (3 * x2 * x292 + x232) + x192 * x286 + x194 * x286 - x2 * x286 * x297 + x202 * x89 + x208 * x276 * x289 - x208 * x293 * x297 * xi_x + 30 * x210 * x279 - x214 * x285 * x299 - x228 * x271 - x261 * x267 * x70 + x261 * xi_xx + x266 ^ 2 * x93 + x266 * x284 * x69 * x70 - 48 * x269 * x297 + x272 * x287 + x272 * x294 + 4 * x273 + x274 * x291 + x274 * x296 + x280 * x290 - 2 * x285 * xi_xx - x287 * x298 + 2 * x288 * x297 - x288 * x59 + x291 * x300 - x294 * x298 + x296 * x300 - 4 * x297 * x89) + x284 * x62 * (Bt * Gt * x156 - Fx * x170 * x173 + Fx * x256 * x268 - Fxt * x121 + Fxt * x243 + Gb * x32 - Gt * x123 + S0_tx * (S0 * x198 * (S0 * x166 - x101 * x283) + x104 * (12 * u * x280 * x282 + u * (-x135 + x148 * x277 - x162 + x184 * x68 + x270 * x6) + x133 * (x7 + 1))) - S0_txx * x167 - S0_xx * x242 - x120 * x273 - x120 * xi_xx + x128 * x145 * xi_x - x132 * x133 - x133 * x171 + x135 * x48 * x76 - x138 * xi_x + 7 * x145 * x89 + x151 * x270 + x154 * x270 + 30 * x157 * x272 - x158 * x271 + x162 * x178 - x163 * xi_x - x170 * x210 * x278 - 36 * x170 * x256 * x280 - x178 * x180 - x217 * x270 + x237 * xi_xx - x239 * x263 + x240 * x267 + x241 * x26 + x247 * x265 + x249 * x274 + x250 * x275 + x251 * x278 * xi_x + x255 * x279 + x257 * x274 + x258 * x275 + x268 * x269 + x281 * x283 * x48)) - x55 * x60 ^ 2) * exp(2 * x7))
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



x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = 8 * x2
x4 = u ^ 2
x5 = S0_t * u
x6 = S0 * xi
x7 = u * x6
x8 = S * x0
x9 = S0 + x7 + x8
x10 = x5 + x9
x11 = x10 ^ 3
x12 = x11 * x4
x13 = x10 ^ 2
x14 = G * x0
x15 = tanh(x14)
x16 = 3 * u
x17 = G * x16
x18 = -Gp + x17
x19 = x15 * x18 * x4
x20 = S0 + 2 * x7
x21 = B * x16 - Bp
x22 = x0 * x11 * x21
x23 = S0 ^ 3
x24 = 1 / x23
x25 = sech(x14)
x26 = S0 * u
x27 = x25 * x26
x28 = u ^ 4
x29 = S0 * x28
x30 = 2 * S * x29
x31 = S0 ^ 2
x32 = 2 * x8
x33 = x31 * (u * xi + 1) + x5 * (S0 - x32)
x34 = S * x4
x35 = 2 * xi_y
x36 = x31 * x4
x37 = S0_ty * x30 + S0_y * x33 + x36 * (-Fy * S0 + Fy * x32 + Sh * u + x34 * x35)
x38 = -Fyp
x39 = Fy * u
x40 = -x38 - x39
x41 = x37 * x40
x42 = x38 + 2 * x39
x43 = -x42
x44 = 2 * x5
x45 = S0_ty * x26
x46 = 5 * x4
x47 = x10 * (S0_y * x43 * x44 + x31 * (-2 * Fyph - u * x35 * x43 + u * (Fy ^ 2 * x46 + 2 * Fyh + Fyp ^ 2 - 4 * Fyp * x39)) + 2 * x42 * x45)
x48 = 2 * x27
x49 = 4 * x2
x50 = S0 ^ 4
x51 = 2 * x23
x52 = x13 * (Sd * x0 * x51 + u * x51 * (S0_t + x6) + x4 * (-S0 * S0_xx - S0 * S0_yy + S0_t ^ 2 * x31 + S0_t * x51 * xi + S0_x ^ 2 + S0_y ^ 2 + x50 * xi ^ 2) + x50)
x53 = Fx * u
x54 = 2 * xi_x
x55 = S0_tx * x30 + S0_x * x33 + x36 * (-Fx * S0 + Fx * x32 + St * u + x34 * x54)
x56 = -Fxp
x57 = 2 * x53 + x56
x58 = -x57
x59 = S0_tx * x26
x60 = (x10 * (S0_x * x44 * x58 + x31 * (-2 * Fxpt - u * x54 * x58 + u * (Fx ^ 2 * x46 + Fxp ^ 2 - 4 * Fxp * x53 + 2 * Fxt)) + 2 * x57 * x59) + x55 * (4 * Fxp - 4 * x53)) * exp(2 * x1)
x61 = Gh * x31
x62 = Gp * x31
x63 = Fx * Gp
x64 = S0_y * x5
x65 = Gt * x31
x66 = 3 * G
x67 = Fxp * x66
x68 = Fy * Gp
x69 = S0_x * x5
x70 = Fyp * x66
x71 = x17 * x31
x72 = u * x31
x73 = Fxp * xi_y
x74 = Fyp * xi_x
x75 = Fx * x36
x76 = Fy * x36
x77 = 2 * x2
x78 = 4 * u
x79 = x2 * x78
x80 = sinh(x14)
x81 = x26 * x80
x82 = cosh(x14)
x83 = Fxp * x72
x84 = 2 * x4
x85 = S0_t * S0_y
x86 = Fx * x85
x87 = S0_t * S0_x
x88 = Fy * x87
x89 = Bp * x28
x90 = x0 * x31
x91 = Bp * x90
x92 = x28 * x31
x93 = Fy * x92
x94 = 3 * B
x95 = Fxp * x94
x96 = Fyp * x94
x97 = x92 * x94
x98 = Fx * x92
x99 = u ^ 5 * x31
x100 = 3 * x1
x101 = Bp * x4
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x12 * x3
BB[1,2] = 0
CC[1,1] = -u * x13 * x3 * (-Sp * x4 + x19 * x9 + x20 + x5 * (x19 + 2) + 4 * x8)
CC[1,2] = -x15 * x22 * x3
SS[1] = x24 * (x10 * x25 * x29 * x49 * (-Fxh * x62 + Fxh * x71 + Fxp * x61 - Fyp * x65 + Fyt * x62 - Fyt * x71 + x39 * x62 * xi_x + x39 * x65 + x45 * (3 * Fxp * G - x63) - x53 * x61 + x59 * (x68 - x70) + x63 * x64 - x63 * x72 * xi_y - x64 * x67 + x67 * x76 - x68 * x69 + x69 * x70 - x70 * x75 + x71 * x73 - x71 * x74) - x21 * x49 * x52 + 8 * x27 * x41 + x47 * x48 - x48 * x60)
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x12 * x49
CC[2,1] = x22 * x77 * sinh(2 * x14)
CC[2,2] = -x13 * x79 * (x20 + x4 * (S * x78 - Sp) + x44)
SS[2] = -x24 * (S0 * x41 * x78 * x80 - S0 * x79 * x82 * (x37 * (-x53 - x56) + x40 * x55) + x10 * x26 * x77 * x82 * (Bh * Fxp * x90 - Bh * x98 + Bp * x93 * xi_x - Bp * x98 * xi_y - Bt * Fyp * x90 + Bt * x93 - 5 * Fx * Fy * x90 - Fx * x96 * x99 - Fxh * x72 - Fxh * x91 + Fxh * x97 - Fxp * x64 + 2 * Fxp * x76 + Fxph * x31 + Fy * x95 * x99 - Fyp * x69 + 2 * Fyp * x75 - Fyp * x83 + Fypt * x31 - Fyt * x72 + Fyt * x91 - Fyt * x97 - x28 * x85 * x95 + x28 * x87 * x96 - x35 * x75 + x45 * (Fxp * x100 + Fxp - x53 * (x101 + 2)) - x54 * x76 + x59 * (-Fyp * x100 + Fyp + x39 * (x101 - 2)) + x72 * x74 + x73 * x97 - x74 * x97 + x83 * xi_y + x84 * x86 + x84 * x88 + x86 * x89 - x88 * x89) + x18 * x52 * x77 + x47 * x81 + x60 * x81)
    

    
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
	#println("entered A u zero")
	ABCS[4] = -4*S0^3*Sp
else
	ABCS[4] = x6 * (S0 * x19 * x5 * x9 + S0_t * x17 + u * x37 * x50 * (-x228 * x61 * (B * x106 * x255 + Bh * Gh * x156 + Fy * x173 * x220 + Fyh * x121 - Fyh * x243 + Fyp * x161 + Gh * x124 - Gs * x32 - S0_t * S0_yy * x240 + S0_ty * (S0 * S0_y * (-x101 * x248 * x260 + x182) + x104 * (6 * Bh * x52 + Fyp * x83 + x131 * (x7 - 1) + x183 * x248 * x260 + x64 * (Gp + x129 * (x261 - 5)))) + S0_tyy * x167 + S0_yy * x242 + x103 * x239 + x106 * x120 * x2 - x106 * x247 + x120 * xi_yy + x126 * x229 * x260 + x131 * x134 - x131 * x148 * x220 - x139 * x48 * x71 + x140 * xi_y + x144 * x244 + x149 * x244 - 30 * x160 * xi_y - x164 * x181 + x165 * xi_y + x179 * x181 - x18 * x220 * x245 + x200 * x245 + x203 * x245 - x210 * x220 * x253 - x220 * x252 * x256 * xi_y - x231 * x250 + x231 * x258 - x235 * x249 + x235 * x257 - x236 * x238 - x237 * xi_yy - 7 * x238 * x86 - x241 * x27 + x251 * x253 * xi_y) + x53 * (x109 * x227 * x228 + x119 ^ 2 * x92 + x119 * x2 * x223 * x66 + x122 * x66 ^ 2 + x223 * (u * x224 * x39 + x225 * x85 - x30 * (Fyph - u * (Fyh + x224 * x64) + x183 * x225)) - x228 * (12 * B * x229 + x128 * x39 * (-x232 * x85 + x234 * x30) + x230 * x231 - 6 * x234 * x30 * x85 + x32 * (Bs + x230 * x235 + x233 * x236 + 6 * x86 * (Bh + x143 * x228))) + x87 ^ 2 * x92)) + x1 * x17 - x11 * x186 * x53 * x81 * (3 * B * Fx * Gh * x32 * x91 + 3 * B * Gh * x18 * x32 * xi_x + 3 * B * Gt * S0_t * S0_y * x18 * x30 + 3 * Bh * G * S0_t * S0_x * x18 * x30 - Bh * Gt * x2 * x32 + 3 * Bt * Fy * G * x32 * x91 + 3 * Bt * G * S0_ty * x18 * x33 + 3 * Bt * G * x18 * x32 * xi_y + Bt * Gh * x2 * x32 - Bt * x172 * x18 * x71 + 36 * Fx * Fy * G * x18 * x32 - Fx * Fyp * x158 + 30 * Fx * G * S0_ty * x2 * x33 + 30 * Fx * G * x2 * x32 * xi_y + 7 * Fx * Gh * x32 * x6 + Fx * Gp * S0_t * S0_y * x30 * x6 + 3 * Fxh * G * x32 * x6 - Fxh * x121 + 3 * Fxp * G * S0_t * S0_y * x30 * x6 - Fxp * x161 + 30 * Fy * G * x2 * x32 * xi_x + Fy * Gp * S0_t * S0_x * x30 * x6 + 7 * Fy * Gt * x32 * x6 - Fy * x135 * x156 + 3 * Fyp * G * S0_t * S0_x * x30 * x6 + 3 * Fyt * G * x32 * x6 - Fyt * x121 + 12 * G * S0 * S0_t * S0_x * S0_y * u + 24 * G * S0_ty * x33 * x6 * xi_x + 24 * G * S0_x * S0_y * x31 * x6 + 6 * G * u * x32 * xi_xy + 24 * G * x32 * x6 * xi_x * xi_y + 2 * Gc * x32 + 6 * Gh * u * x32 * xi_x - Gh * x123 + 2 * Gp * S0_t * S0_xy * x30 + 2 * Gp * S0_ty * S0_x * x30 - Gp * x125 * x20 + 6 * Gt * S0_ty * u * x33 + 6 * Gt * u * x32 * xi_y - Gt * x124 - 3 * Gt * x147 * x148 + S0 * S0_tx * (S0_y * (-x101 * x141 - x182) + u * x30 * (-u * (Bh * x185 + x139 + x164 - x184 * x64) + x141 * x183 + x169 * (x7 + 2)) + x142) - 2 * x120 * xi_xy - x126 * x127 * x30 - x129 * x130 * x30 - x131 * x132 - x133 * x134 - x135 * x136 - x136 * x162 - x138 * xi_y - x140 * xi_x - x142 * x76 - x144 * x146 - x146 * x149 - x151 * x152 - x152 * x154 - x163 * xi_y - x165 * xi_x - x167 * x168 - x169 * x171 - x173 * x175 - x173 * x177 - x178 * x179 - x180 * x181) + x14 * x25 * x9 / x2 + 4 * x15 * x16 * x186 * (-x2 * x53 * (x49 * x90 + x60 * x87) + x61 * (Fx * Fy * Sp * x155 * x210 - Fx * x211 * x212 + Fxh * x189 + Fxh * x191 - Fxh * x197 + Fyt * x189 + Fyt * x191 - Fyt * x197 - 4 * S * x153 * xi_xy - 8 * S0 * x125 * x41 + S0_ty * x190 * x32 - S0_ty * x198 * x219 - S0_ty * x206 * x216 - S0_x * S0_y * x206 * x31 + S0_xy * x214 * x215 - S0_xy * x219 * x34 - Sc * x156 - Sh * x153 * x202 + Sp * x151 * xi_y + Sp * x156 * xi_xy + Sp * x200 * xi_x + x107 * x168 * x23 - x107 * x174 - x107 * x176 + x108 * x32 * x89 + 4 * x125 * x2 * x207 - 2 * x127 * x33 + 4 * x130 * x215 - x147 * x201 + x147 * x204 - x147 * x209 * x211 - x151 * x195 - 16 * x159 * x208 * x209 - x159 * x211 * x213 + x170 * x206 * xi_y + x175 * x221 - x176 * x218 * x91 + x177 * x221 - x187 * x42 - x187 + x188 * x190 * xi_y + x192 * x193 + 2 * x193 * x194 + x195 * x217 + x198 * x32 * x86 - x199 * x200 - x199 * x203 + x2 * x56 * (S0_y * (S0_t * x222 + x20 * x22 - 2 * x207) - x222 * x39 + x30 * (Fy * (S0 + x6 * (-x21 - 14 * x22)) - x222 * xi_y - 4 * x45)) + x201 * x220 - x204 * x220 - x205 * x212 * xi_x + x205 * x213 * x220 + x206 * x39 * x76)) - 12 * x16 ^ 4 * x8 + x16 * (4 * x49 * x50 * x6 * x61 * x87 - x54 * (-S0 * xi_yy - S0_tyy - S0_y * x47 - S0_yy * x15 - S0_yy * xi + x109 * x24 * x35 + x119 * x49 * x50 * x6 + x50 * (x30 * x74 + x72) * (S0_y * x113 + x112 * x40 + x48 * (Sph + x114 * x47 - x81 * (Sh + x110 * x64))) - x6 * (-Spp * x75 ^ 2 + Ss + (Sph + Spp * x75) * (x47 + 2 * x73 + x94)) + x75 * (S0_y * x102 + x100 * x95 + x6 * (Sph + u * (Fy * Spp - Fy * x96 - 2 * Sh) + x98 * xi_y)))) + x19 * x20 * x25 * x8 + 4 * x36 * x38 * (x26 + x27 - x28 - x29 + x30 * x31 + x32 * xi ^ 2 + x33 * x34 * xi) + x38 * x61 * x81 * (Fxph + Fxpp * x75 + Fypp * x80 + Fypt + u * x66 * x69 - u * (Fxh + Fxp * x75) - u * (Fyp * x80 + Fyt) - x75 * (Fxpp + x69 * x81) - x80 * (Fypp + x65 * x81) - x87 * x90 * x93) - x49 ^ 2 * x55 + x49 * x51 * x60 * x62 * x9 + (x16 * (4 * x262 * x61 * x90 + x54 * (S0 * xi_xx + S0_txx + S0_x * x59 + S0_xx * x15 + S0_xx * xi + x262 * x266 - x36 * (Fxt * x104 - S0_t * x28 + S0_txx * x30 + x107 * x265 + x216 * x89 + x263 - x264 * x76 + x33 * xi_xx + x56 * (-x198 + x264)) - x50 * (x30 * x79 + x77) * (S0_x * x113 + x112 * x57 + x48 * (Spt + x114 * x59 - x81 * (St + x110 * x68))) + x6 * (Sb - Spp * x80 ^ 2 + (Spp * x80 + Spt) * (x59 + 2 * x68 + 2 * x78)) - x80 * (S0_tx * x100 * x6 + S0_x * x102 + x6 * (Spt + u * (Fx * Spp - Fx * x96 - 2 * St) + x98 * xi_x)))) + x37 * (u * x53 * (Bb * x228 - Bp * Fxt * x284 - Bp * x265 * x92 - Fxp * u * x59 + Fxp * x0 * x198 * x70 - 2 * Fxpt + 6 * Fxt * x148 + Fxt * x81 - G * x290 * x297 + Gt ^ 2 * x92 + 12 * S0_t * x299 * x7 + S0_tx * x35 * x81 * (u * x198 * (-S0_t * (x230 + 9 * x292 * x91) + x227) + x30 * (u * (Fx * (-x285 + x301 + 27 * x7 + 2) + x128 * (Bt + Gt * x52)) + x67 + xi_x * (x301 + 24 * x7))) + S0_txx * x226 * x228 * x99 + S0_xx * x285 * x34 * x70 + u ^ 9 * x265 * x295 + u * x69 ^ 2 + 6 * x18 * x281 * x70 * (3 * x2 * x292 + x232) + x192 * x286 + x194 * x286 - x2 * x286 * x297 + x202 * x89 + x208 * x276 * x289 - x208 * x293 * x297 * xi_x + 30 * x210 * x279 - x214 * x285 * x299 - x228 * x271 - x261 * x267 * x70 + x261 * xi_xx + x266 ^ 2 * x93 + x266 * x284 * x69 * x70 - 48 * x269 * x297 + x272 * x287 + x272 * x294 + 4 * x273 + x274 * x291 + x274 * x296 + x280 * x290 - 2 * x285 * xi_xx - x287 * x298 + 2 * x288 * x297 - x288 * x59 + x291 * x300 - x294 * x298 + x296 * x300 - 4 * x297 * x89) + x284 * x62 * (Bt * Gt * x156 - Fx * x170 * x173 + Fx * x256 * x268 - Fxt * x121 + Fxt * x243 + Gb * x32 - Gt * x123 + S0_tx * (S0 * x198 * (S0 * x166 - x101 * x283) + x104 * (12 * u * x280 * x282 + u * (-x135 + x148 * x277 - x162 + x184 * x68 + x270 * x6) + x133 * (x7 + 1))) - S0_txx * x167 - S0_xx * x242 - x120 * x273 - x120 * xi_xx + x128 * x145 * xi_x - x132 * x133 - x133 * x171 + x135 * x48 * x76 - x138 * xi_x + 7 * x145 * x89 + x151 * x270 + x154 * x270 + 30 * x157 * x272 - x158 * x271 + x162 * x178 - x163 * xi_x - x170 * x210 * x278 - 36 * x170 * x256 * x280 - x178 * x180 - x217 * x270 + x237 * xi_xx - x239 * x263 + x240 * x267 + x241 * x26 + x247 * x265 + x249 * x274 + x250 * x275 + x251 * x278 * xi_x + x255 * x279 + x257 * x274 + x258 * x275 + x268 * x269 + x281 * x283 * x48)) - x55 * x60 ^ 2) * exp(2 * x7))
	
end


    nothing
    
end
