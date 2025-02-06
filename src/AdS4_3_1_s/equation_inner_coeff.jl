
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
x11 = x10 * x3
x12 = 2 * x11
x13 = 3 * x1
x14 = Bp * x3
x15 = -x14
x16 = G * x0
x17 = 2 * x16
x18 = cosh(x17)
x19 = -Bp
x20 = 3 * u
x21 = B * x20
x22 = x19 + x21
x23 = x22 * x3
x24 = x13 + x15 + x18 * x23
x25 = u * x2
x26 = x10 * x25
x27 = 2 * Gp
x28 = G * u
x29 = 6 * x28
x30 = -x22
x31 = sinh(x17)
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
x77 = B * x54
x78 = B * x63
x79 = x55 * x69
x80 = x61 * xi
x81 = x46 * x80
x82 = x57 * x61
x83 = x82 * xi
x84 = Bpp * x60
x85 = x66 * xi
x86 = G * Gp
x87 = Bp ^ 2
x88 = u ^ 9
x89 = x45 * x88
x90 = x0 * x37
x91 = u ^ 10
x92 = x75 * x91
x93 = Bp * x55
x94 = x37 * x57
x95 = 24 * x41
x96 = B ^ 2
x97 = 18 * S0
x98 = S * x71
x99 = 12 * x41
x100 = x63 * x66
x101 = 12 * G
x102 = Gp * x101
x103 = x45 * x91
x104 = G ^ 2
x105 = Gp ^ 2
x106 = S * x72
x107 = G * S
x108 = Gp * x71
x109 = x107 * x108
x110 = u ^ 11 * x45
x111 = 9 * x96
x112 = x46 * x73
x113 = 18 * x104
x114 = x37 * x46
x115 = x105 * x89
x116 = x105 * x90
x117 = Bp * x46
x118 = 18 * x96
x119 = S * x69 * x88
x120 = x54 * x87
x121 = Gp * x46
x122 = G * x121
x123 = 36 * x104
x124 = x105 * x54
x125 = x35 * xi
x126 = 15 * x57
x127 = x37 * x76
x128 = B * x127
x129 = x41 * x49
x130 = Bp * x90
x131 = x57 * x87
x132 = x105 * x57
x133 = 4 * x41
x134 = x46 * x87
x135 = x127 * x49
x136 = x53 * x73 * x76
x137 = x105 * x114 * x76
x138 = x107 * x72
x139 = 12 * S0
x140 = Gp * x31
x141 = x31 * x61
x142 = x101 * x141
x143 = 18 * G
x144 = x143 * x31
x145 = B * x110
x146 = x140 * x55
x147 = 6 * G
x148 = x147 * x31
x149 = Bp * x31
x150 = x149 * x94
x151 = x149 * x89
x152 = x27 * x31
x153 = B * x69
x154 = x107 * x153 * x31 * x88
x155 = x31 * xi
x156 = G * x31
x157 = B * x129
x158 = x31 * x41
x159 = x41 * x57
x160 = x135 * x149
x161 = x117 * x127
x162 = Bpp * u
x163 = 15 * x1
x164 = 7 * x14
x165 = x101 * x121
x166 = 9 * x49
x167 = x166 * x96
x168 = x113 * x49
x169 = 2 * x132
x170 = -Gp
x171 = G * x20
x172 = x170 + x171
x173 = -x172
x174 = x32 * x57
x175 = x173 * x174
x176 = 2 * x175
x177 = x3 * x87
x178 = 6 * x1
x179 = x13 + 5
x180 = Bpp + u * (-Bp * (x178 + 7) + x177 + x179 * x21)
x181 = u * x18
x182 = S0_t ^ 2 * u
x183 = Bpp * x8
x184 = x177 * x8
x185 = 7 * x5
x186 = x5 + 1
x187 = x178 * x186
x188 = x185 + x187
x189 = -x62
x190 = S * u
x191 = 5 * x5
x192 = x13 * x186 + x191
x193 = S * x20
x194 = x18 * x8
x195 = 8 * x70
x196 = Spp * u
x197 = 16 * x7
x198 = S0 * x131
x199 = -Sp
x200 = 9 * x190
x201 = x55 * x60
x202 = x13 + 7
x203 = 9 * S
x204 = 6 * S0
x205 = x203 * x88
x206 = S0 * x167
x207 = S * x113
x208 = x207 * x88
x209 = Bpp * S0
x210 = x3 * xi
x211 = x55 * x61
x212 = x211 * x46
x213 = x105 * x46
x214 = x53 * x69
x215 = x49 * xi
x216 = x139 * x86
x217 = x35 * x5 + x35
x218 = 2 * S0_t
x219 = 2 * Gpp
x220 = S0 * x219
x221 = x18 * x30
x222 = 2 * x0
x223 = x173 * x222
x224 = x221 * x223
x225 = -x13 + x14
x226 = 2 * u
x227 = Bpp - u * (Bp * (7 - x178) + x177 + x21 * (x13 - 5))
x228 = 22 * Gp * x60
x229 = x101 * x58
x230 = Gp * S0
x231 = 10 * u * x230
x232 = x219 * x7
x233 = 4 * x0
x234 = Gp * Sp * x233
x235 = G * x97
x236 = x235 * x3
x237 = 54 * x107 * x46
x238 = x31 * x87
x239 = x238 * x63
x240 = S0 * x0 * x238
x241 = 14 * x210 * x230
x242 = x31 * x71
x243 = x203 * x242 * x96
x244 = S0 * x31
x245 = x111 * x244 * x46
x246 = x219 * x6
x247 = 30 * x16 * x69
x248 = x155 * x198
x249 = x155 * x206
x250 = x31 * x54 * x93
x251 = x211 * x31 * x57
x252 = x155 * x212
x253 = x147 * x54
x254 = x31 * x69
x255 = Bp * x253 + Bpp * x31 * x6 + Gp * x56 - 5 * u * x141 + x0 * x149 * x62 + x121 * x79 - 18 * x138 - x143 * x153 * x49 + x147 * x81 + x147 * x82 - 11 * x149 * x60 + x163 * x254 - x164 * x254 + x194 * x223 * x30 + x209 * x31 + x230 * x55 * x57 - x235 * x51 + 9 * x244 * x74 - x27 * x64 - x27 * x65 - x27 * x83 + 27 * x31 * x52 - x31 * x59 + x31 * x67
x256 = u * x9
x257 = x10 * x226
x258 = 2 * x190 + x199
x259 = -x258 * x3
x260 = S0 + x259
x261 = 4 * x2
x262 = x260 * x261
x263 = 6 * x11
x264 = x260 * x35
x265 = x2 * xi_x
x266 = S0_x * xi
x267 = x260 ^ 2
x268 = x261 * x267
x269 = 1 / u
x270 = S0_x * x269
x271 = x25 * x9
x272 = 4 * x271
x273 = u * x10
x274 = x273 * x31
x275 = x261 * x9
x276 = x3 * x9
x277 = 3 * x31
x278 = 2 * x10
x279 = x278 * x57
x280 = x172 * x279
x281 = 1 / x37
x282 = x281 * (S0 * S0_tx - S0_t * S0_x)
x283 = x170 + 2 * x28
x284 = 6 * u
x285 = Gpp + x283 * x284
x286 = x257 * x285
x287 = cosh(x16) ^ 2
x288 = 2 * x26 * x287
x289 = B_x * x2
x290 = x272 * (Spp - 4 * x39 + 6 * x43)
x291 = x30 ^ 2
x292 = x291 * xi_y
x293 = x10 * x31
x294 = x293 * x57
x295 = x280 * xi_y
x296 = B * x226 + x19
x297 = Bpp + x284 * x296
x298 = x274 * x297
x299 = x22 * x287
x300 = S0_y * x9
x301 = x23 * x287
x302 = x275 * x301
x303 = 2 * x33
x304 = x276 * x303
x305 = S_y * x9
x306 = x281 * (S0 * S0_ty - S0_t * S0_y)
x307 = x22 * x289
x308 = x278 * x287
x309 = x308 * x57
x310 = x265 * x291
x311 = x291 * x306
x312 = x288 * x297
x313 = x22 * x280 * x306
x314 = x2 * x278
x315 = x57 * (-Gp + x171 - x32)
x316 = x314 * x315
x317 = x2 * x282
x318 = x291 * x317
x319 = -S0 + x258 * x3
x320 = x319 * xi_x
x321 = x319 * xi_y
x322 = x314 * xi_x
x323 = x282 * x319
x324 = x27 - x29 - x32
x325 = x2 * x324
x326 = Spp * x66
x327 = 24 * x43
x328 = Gp * x147
x329 = 9 * x104
x330 = x141 * x27
x331 = Gp * x277
x332 = B * x331
x333 = 3 * G
x334 = 6 * x107
x335 = x104 * x166
x336 = u * x287
x337 = x185 - x187
x338 = -3 * B * x0 * x186 + x191
x339 = 4 * x260
x340 = x269 * xi_y
x341 = 4 * x259 + x35
x342 = x269 * x341
x343 = S0_y / x3
x344 = 4 * x267
x345 = x217 + 4 * x4 + 4 * x7
x346 = x10 * x287
x347 = x9 * (x327 - 16 * x39 + x44)
x348 = x2 * x293
x349 = 12 * u
x350 = x287 * x300
x351 = x256 * x299
x352 = 4 * x351
x353 = 2 * x34
x354 = x172 * x353
x355 = x287 * x353
x356 = x346 * (2 * Bpp + x296 * x349)
x357 = x31 * x34
x358 = x353 * (x172 + x32)
x359 = x297 * x348
x360 = x22 * x354
x361 = x265 * x360
x362 = x325 * x9
x363 = S0_x * x362
x364 = x271 * x324
x365 = 2 * x364
x366 = x172 * x358
x367 = x317 * x360
AA[1,1] = x12 * x2
AA[1,2] = 0
BB[1,1] = x26 * (x24 + 8)
BB[1,2] = -x33 * x34
CC[1,1] = x25 * (-B * x117 * x99 + 4 * Bp * x140 * x159 - 5 * Bp * x42 + Bpp * x45 * x49 - Gp * x106 * x139 * x155 + Gp * x125 * x149 * x54 - 18 * S * x81 + 48 * S0 * x43 - S0 * x59 - Sp * x46 * x79 - Sp * x56 + x1 * x95 + x100 * x87 - x101 * x117 * x158 - x102 * x103 - x102 * x135 - x102 * x94 - x103 * x146 + x104 * x50 * x98 + x105 * x35 * x63 - 12 * x106 * x80 - 24 * x109 * x69 + x110 * x111 + x110 * x113 + x112 * x96 + x113 * x114 + x113 * x127 * x53 + 2 * x115 + 2 * x116 + x118 * x119 + x118 * x129 + x119 * x123 + x120 * x85 - x122 * x95 + x123 * x129 + x124 * x125 + x126 * x128 + x127 * x134 + x128 * x144 * x53 + x130 * x152 - 7 * x130 * x76 + 2 * x131 * x41 + x132 * x133 - x135 * x146 - x135 * x93 + x136 * x96 + 2 * x137 + x138 * x31 * x50 - x139 * x140 * x77 - x14 * x99 + x140 * x35 * x64 - x140 * x51 * x99 - x142 * x54 - x142 * x98 * xi + x144 * x145 + x144 * x37 * x51 - x146 * x94 - x147 * x150 - x147 * x160 - x148 * x92 + x151 * x27 + x152 * x161 + 36 * x154 + 36 * x156 * x157 + x182 * (-x117 * x55 + x131 + x162 + x163 - x164 - x165 + x167 + x168 + x169 + x176 + x180 * x181 + 4) + x194 * (u * (-Bp * (S0 * (x188 + 5) + x3 * (x189 + x190 * (x178 + 11))) + x184 + x21 * (S0 * (x192 + 3) + x3 * (x189 + x193 * (x1 + 3)))) + x183) + x218 * (Bp * x58 + S0 * x162 + S0 * x168 - 3 * Sp * x51 + x1 * x139 - x102 * x98 + x111 * x214 + x113 * x214 - x117 * x203 + x120 - x122 * x139 + 2 * x124 + x126 * x153 + x132 * x66 + x134 * x69 - x14 * x204 + x176 * x8 + x181 * (u * (-Bp * (S0 * (x188 + 6) + x3 * (x199 + x200 + x201)) + x184 + x21 * (S0 * (x192 + 4) + x3 * (x190 * x202 + x199))) + x183) - x195 + 2 * x196 + x197 + x198 + x205 * x96 + x206 + x208 + x209 * x210 - x211 * x215 - x212 + x213 * x85 - x215 * x216 + x217 - 7 * x65 * xi + 21 * x78 + x84 - x93 * x98) + x3 * x38 * x76 + x36 * x5 + x36 + x38 * x68 + x38 - x39 * x40 - x40 * x54 * x86 + 8 * x41 + 4 * x42 * x76 + x44 * x7 + 27 * x45 * x72 + 12 * x47 - 4 * x48 + x50 * x52 - 11 * x53 * x75 - x55 * x92 - 16 * x60 * x61 - 12 * x61 * x77 + x62 * x64 + x62 * x65 + x62 * x83 + x66 * x67 + 32 * x69 * x7 - 16 * x69 * x70 + 42 * x69 * x78 + x73 * x74 + x84 * x85 + x87 * x89 + x87 * x90 - x93 * x94 + x96 * x97 * x98)
CC[1,2] = -x256 * (x220 - x228 - x229 - x231 + x232 + x234 + x236 + x237 - x239 - x240 - x241 - x243 - x245 + x246 + x247 - x248 - x249 + x250 + x251 + x252 + x255 + x4 * (x219 + x224 + x226 * (-Gp * (x225 + 7) + x171 * (x225 + 5)) + x227 * x31))
if u==0.0
	SS[1] = -2*S0^2*(Bp_x) + 2*S0^2*(Gp_y) + 2*Bpp*S0*(S0_tx) + 4*Spp*(S0_tx) - 2*Gpp*S0*(S0_ty) - 4*Bp*S0*(S0_x) + 4*Sp*(S0_x) - 2*Bpp*(S0_t)*(S0_x) - (4*Spp*(S0_t)*(S0_x))/S0 + 4*Gp*S0*(S0_y) + 2*Gpp*(S0_t)*(S0_y) - 4*S0*(Sp_x) + 2*Bpp*S0^2*(xi_x) + 4*S0*Spp*(xi_x) - 2*Gpp*S0^2*(xi_y)
else
	SS[1] = -B_y * x10 * x174 - B_y * x11 * x277 + B_y * x280 - Bp_x * x288 + Bp_y * x274 + G_x * x316 + G_y * x221 * x279 - G_y * x263 + Gp_y * x257 + S0_tx * x262 + S0_tx * x302 - S0_ty * x304 + S0_x * x272 * x299 + 8 * S_x * x2 * x276 + S_x * x262 * x3 + S_x * x275 * x299 * x57 - Sp_x * x272 + x172 * x282 * x316 + x172 * x315 * x322 - x18 * x313 - x210 * x300 * x303 + x22 * x295 + x221 * x295 - x226 * x300 * x33 + x260 * x304 * x306 + x262 * x266 + x262 * x270 + x263 * x287 * x289 + x264 * x265 + x265 * x301 * x35 * x9 + x266 * x302 - x268 * x282 - x268 * xi_x - x270 * x275 - x276 * x33 * x66 * xi_y + x282 * x290 + x282 * x312 - x286 * x306 - x286 * xi_y + x290 * xi_x + x292 * x294 + x294 * x311 - x298 * x306 - x298 * xi_y + x302 * x320 + x302 * x323 - x303 * x305 * x57 - x304 * x321 + x307 * x309 + x309 * x310 + x309 * x318 + x312 * xi_x + x313
end
AA[2,1] = 0
AA[2,2] = x12
BB[2,1] = x325 * x34
BB[2,2] = -x273 * (x24 - 8)
CC[2,1] = x271 * (-x220 + x228 + x229 + x231 - x232 - x234 - x236 - x237 + x239 + x240 + x241 + x243 + x245 - x246 - x247 + x248 + x249 - x250 - x251 - x252 + x255 - x4 * (-x180 * x31 + x219 - x224 + x226 * (Gp * (-x15 - x202) + x171 * (x15 + x179))))
CC[2,2] = x226 * (-B * x136 * x156 + G * x277 * x92 - Gp * x151 + S0 * x207 * x71 + S0 * x327 + 2 * Spp * x7 + x100 * x105 - x103 * x328 + x103 * x332 + x104 * x112 + x104 * x136 - x106 * x235 * x31 - x109 * x139 * xi + x110 * x329 + x115 + x116 + x117 * x148 * x41 + x121 * x158 * x55 + x124 * x85 + x127 * x226 + x128 * x331 * x49 - x130 * x140 + x133 - x135 * x328 + x137 - x139 * x39 - x140 * x161 + x141 * x253 - x144 * x157 - 9 * x145 * x156 + x146 * x69 * x98 - x149 * x159 * x27 + x150 * x333 - 18 * x154 - x156 * x51 * x73 + x160 * x333 - x165 * x41 + x168 * x41 + x169 * x41 + x182 * (-6 * x122 + x132 - x175 - x227 * x336 + x335 + 2) - x195 * x69 + x197 * x69 + x208 * x69 - x216 * x54 + x218 * (S0 * x132 + S0 * x335 + x104 * x205 - x108 * x334 - x122 * x204 + x124 - x175 * x8 + x196 - x204 * x215 * x86 + x213 * x69 + x214 * x329 + x336 * (u * (-Bp * (S0 * (-x337 - 6) + x3 * (Sp - x200 + x201)) + x184 + x21 * (S0 * (-x338 - 4) + x3 * (Sp + x190 * (x13 - 7)))) - x183) + x5 * x66 + x66 + 8 * x7 - 4 * x70) + x230 * x31 * x56 + x242 * x334 * x80 + x287 * x8 * (u * (S0 * x21 * (-x338 - 3) + x13 * (x193 * (x1 - 3) + x62) + x14 * (x189 + x190 * (11 - x178)) + x184 + x61 * (x337 + 5)) - x183) + x326 * x5 + x326 - x328 * x94 - x330 * x54 * xi - x330 * x63 + x332 * x94 + 6 * x47 - 2 * x48)
if u==0.0
	SS[2] = 2*S0^2*(Bp_y) + 2*S0^2*(Gp_x) - 2*Gpp*S0*(S0_tx) - 2*Bpp*S0*(S0_ty) + 4*Spp*(S0_ty) + 4*Gp*S0*(S0_x) + 2*Gpp*(S0_t)*(S0_x) + 4*Bp*S0*(S0_y) + 4*Sp*(S0_y) + 2*Bpp*(S0_t)*(S0_y) - (4*Spp*(S0_t)*(S0_y))/S0 - 4*S0*(Sp_y) - 2*Gpp*S0^2*(xi_x) - 2*Bpp*S0^2*(xi_y) + 4*S0*Spp*(xi_y)
else
SS[2] = B_x * x20 * x348 - B_y * x284 * x346 + B_y * x299 * x353 - Bp_x * x348 + Bp_y * x308 + G_x * x18 * x2 * x22 * x353 - 6 * G_x * x26 + G_y * x358 + Gp_x * x314 + S0_tx * x365 + S0_ty * x342 - S0_ty * x352 + S0_y * x342 * xi + S_x * x222 * x362 + S_y * u * x339 + 8 * S_y * x256 - Sp_y * x345 - x10 * x317 * (x219 + x283 * x349) + x18 * x361 + x18 * x367 - 4 * x22 * x350 * x5 - x233 * x299 * x305 + x264 * x340 - x269 * x306 * x344 + x282 * x359 - x285 * x322 - x289 * x354 + x292 * x355 + x306 * x339 * x351 + x306 * x347 - x306 * x356 + x306 * x366 + x307 * x357 + x310 * x357 + x311 * x355 + x318 * x357 + x320 * x365 - x321 * x352 + x323 * x365 - x340 * x344 + x341 * x343 - x343 * x345 + x347 * xi_y - x35 * x351 * xi_y + x350 * (-B * x349 + 4 * Bp) - x356 * xi_y + x359 * xi_x - x361 + x363 * x68 + 2 * x363 + x364 * x66 * xi_x + x366 * xi_y - x367
end
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

x0 = S0_t * u
x1 = S0 * xi
x2 = u ^ 3
x3 = S * x2
x4 = S0 + u * x1 + x0 + x3
x5 = x4 ^ 3
x6 = B * x2
x7 = exp(x6)
x8 = u ^ 2
x9 = 8 * x8
x10 = x7 * x9
x11 = x4 ^ 2
x12 = Sp * u
x13 = S * x8
x14 = S0_t + x1
x15 = 1 / u
x16 = S0 * (x15 + xi) + S0_t + x13
x17 = 8 * x16 ^ 3 * x7
x18 = 8 * S0
x19 = u ^ 4
x20 = x7 / x19
x21 = -Sp
x22 = S * u
x23 = x8 * (-x21 - 2 * x22)
x24 = S0 + x23
x25 = 4 * S0
x26 = S0 * S0_ty
x27 = 2 * x26
x28 = S * x19
x29 = S0 ^ 2
x30 = 2 * x3
x31 = x0 * (S0 - x30) + x29 * (u * xi + 1)
x32 = Fy * S0
x33 = 2 * xi_y
x34 = x29 * x8
x35 = S0_y * x31 + x27 * x28 + x34 * (Fy * x30 + Sh * u + x13 * x33 - x32)
x36 = G * x2
x37 = cosh(x36)
x38 = S0 ^ 4
x39 = 1 / x38
x40 = 1 / x8
x41 = 4 * x37 * x39 * x40
x42 = S0 ^ 3
x43 = 1 / x42
x44 = 4 * x23 + x25
x45 = S0_x ^ 2
x46 = S0_y ^ 2
x47 = S0 * S0_xx
x48 = S0 * S0_yy
x49 = S0_t ^ 2
x50 = 2 * S0_t
x51 = x16 ^ 2
x52 = x51 * x7
x53 = S0 * S0_tx
x54 = S0_x * x31
x55 = Fx * S0
x56 = 2 * xi_x
x57 = x34 * (Fx * x30 + St * u + x13 * x56 - x55)
x58 = sinh(x36)
x59 = x39 * x58
x60 = -Fyp
x61 = Fy * u
x62 = x60 + x61
x63 = -x62
x64 = -Fxp
x65 = Fx * u
x66 = x64 + x65
x67 = 1 / x29
x68 = S0_t * S0_y
x69 = x26 - x68
x70 = x61 + xi_y
x71 = x67 * x69 + x70
x72 = S0_t * S0_x
x73 = x53 - x72
x74 = x67 * x73
x75 = x65 + xi_x
x76 = x74 + x75
x77 = 2 * u
x78 = 3 * u
x79 = G * x78
x80 = 3 * G
x81 = S0_y * x0
x82 = Fy * x8
x83 = x26 * x79 + x29 * (Gh + x79 * xi_y + x80 * x82) - x80 * x81
x84 = S0_x * x0
x85 = Fx * x8
x86 = x29 * (Gt + x79 * xi_x + x80 * x85) + x53 * x79 - x80 * x84
x87 = u ^ 5
x88 = 2 * x87
x89 = x39 * x88
x90 = 4 * x8
x91 = x39 * x90
x92 = x35 * x91
x93 = 8 * xi_y
x94 = 4 * S0_yy
x95 = Spp * x71
x96 = Sph + x95
x97 = 2 * x61
x98 = 2 * x12
x99 = -Spp + x98
x100 = -x99
x101 = 1 / S0
x102 = x100 * x101 * x8
x103 = S0_t * x8
x104 = x67 * (x103 * x99 + x29)
x105 = S0_ty * x102 + S0_y * x104 + x8 * (Sph + u * (Fy * Spp - Fy * x98 - 2 * Sh) + x100 * xi_y)
x106 = x46 * x50
x107 = u * x42
x108 = x42 * xi_y
x109 = x32 * x8
x110 = Fy ^ 2
x111 = x2 * x42
x112 = -2 * S0_y + x109
x113 = Fyh * x107 - S0_t * x48 + S0_tyy * x29 + x106 + x108 * x82 - x109 * x68 + x110 * x111 + x112 * x26 + x42 * xi_yy
x114 = 4 * x43
x115 = x114 * x24
x116 = x21 + 3 * x22
x117 = -x116
x118 = x117 * x2
x119 = -x118 * x50 + x29
x120 = u * x117
x121 = x29 * x70
x122 = x121 + x69
x123 = B * x78
x124 = 3 * B
x125 = B * x8
x126 = 3 * x125
x127 = x123 * x26 - x124 * x81 + x29 * (Bh + Fy * x126 + x123 * xi_y)
x128 = 2 * Gpp
x129 = x128 * x71
x130 = x29 * x75
x131 = 2 * x53
x132 = 2 * x72
x133 = x131 - x132
x134 = -Gp * x78 + Gpp
x135 = -Gp
x136 = G * u
x137 = x135 + 4 * x136
x138 = -x137
x139 = 6 * u
x140 = 2 * x29
x141 = 3 * x6
x142 = u * x39
x143 = B * x19
x144 = 3 * Fy
x145 = -x135 - x79
x146 = 2 * S0_x
x147 = x55 * x8
x148 = 2 * Fx * Fy * x111 + Fx * S0_ty * x34 + Fxh * x107 + Fyt * x107 - S0 * S0_xy * x50 + S0_txy * x140 + 4 * S0_x * x68 + x108 * x85 - x109 * x72 + x112 * x53 - x146 * x26 - x147 * x68 + x42 * x82 * xi_x + 2 * x42 * xi_xy
x149 = 2 * x2
x150 = x131 * x28 + x54 + x57
x151 = 8 * S0_xy
x152 = 8 * xi_x
x153 = S0_x * x119 + x118 * x131 + x34 * (Spt + x120 * x56 - x77 * (St + x116 * x65))
x154 = Spp * x76 + Spt
x155 = u * x38
x156 = x60 + x97
x157 = -x156
x158 = u * xi_y
x159 = -Bp + x123
x160 = -S0 * x159
x161 = 2 * x8
x162 = S0_ty ^ 2 * x34
x163 = 12 * x125
x164 = x46 * x49
x165 = 4 * B
x166 = Bh + x165 * x82
x167 = x158 * x165 + x166
x168 = xi_y ^ 2
x169 = x139 * xi_y
x170 = B * Fy
x171 = Gp * x38
x172 = u * x171
x173 = x38 * x79
x174 = Gh * x38
x175 = Gp * S0
x176 = Gp * x29
x177 = 6 * G
x178 = S0 * x0 * x177
x179 = x0 * x29 * x80
x180 = 6 * Gh
x181 = x38 * x8
x182 = x181 * x80
x183 = x145 * x42
x184 = x38 * x87
x185 = x143 * x180
x186 = x38 * xi_y
x187 = Bh * x177
x188 = x19 * x38
x189 = x36 * x38
x190 = 30 * x189
x191 = 18 * G
x192 = x188 * x191
x193 = 12 * G
x194 = x181 * x193
x195 = x193 * x8
x196 = u ^ 6
x197 = 36 * G
x198 = x170 * x196 * x197
x199 = x29 * x68
x200 = 30 * x36
x201 = x34 * x68
x202 = u ^ 7
x203 = x191 * x38
x204 = x202 * x203
x205 = B * x87
x206 = x203 * x205
x207 = x191 * x205
x208 = x141 - 2
x209 = 6 * x6
x210 = x150 * x91
x211 = 4 * S0_xx
x212 = x45 * x50
x213 = Fx * xi_x
x214 = Fx ^ 2
x215 = x123 * x53 - x124 * x84 + x29 * (Bt + Fx * x126 + x123 * xi_x)
x216 = S0_t * S0_xx
x217 = Gt * x38
x218 = 6 * Gt
x219 = x218 * x29
x220 = x218 * x38
x221 = x143 * xi_x
x222 = Bt * x177
x223 = Fx * Fxp
x224 = Fx * Gp
x225 = x181 * xi_x
x226 = Fxp * x80
x227 = x34 * x72
x228 = x2 * x214
x229 = xi_x ^ 2
x230 = x45 * x49
x231 = B * x196
x232 = x29 * x72
x233 = G * xi_x
x234 = B * x214
x235 = Fx * G
x236 = 36 * x235
x237 = S0_tx ^ 2
x238 = x141 + 2
x239 = x177 * x238
x240 = Bp * x8
x241 = 12 * Bt
x242 = Fx * x19
x243 = x2 * x241
x244 = 54 * x205
x245 = Bp * x242
x246 = 12 * Gt
x247 = x196 * x246
x248 = 24 * x143
x249 = G ^ 2
x250 = 36 * x249
x251 = u ^ 8 * x250
x252 = 18 * x249
x253 = x202 * x252
x254 = x67 * x72
x255 = Fx * x254
x256 = S0_t * x45
x257 = x230 * x39
x258 = x196 * x252
ABCS[1] = 0
ABCS[2] = -x10 * x5
ABCS[3] = -x10 * x11 * (-x12 + 3 * x13 + x14)

if u==0.0
	ABCS[4] = 4*S0^3*Sp
else
	ABCS[4] = x8 * (S0_t * x17 + x1 * x17 + x11 * x14 * x7 * (x18 + 8 * x23) / x2 + x11 * x20 * x24 * x25 + x142 * x51 * (-x161 * x58 * (B * x110 * x204 + Bh * x149 * x174 - Fy * Gp * x201 + Fy * x184 * x187 - Fy * x190 * xi_y + Fy * x199 * x200 + Fyh * x172 - Fyh * x182 + Fyp * Gh * x155 + Fyp * x144 * x189 + Fyp * x182 * xi_y - Fyp * x201 * x80 + 24 * G * x201 * xi_y - Gs * x38 - S0_t * S0_yy * x176 + S0_ty * (S0 * S0_y * (-2 * S0 * x145 - x103 * x193 * x208) + x107 * (6 * Bh * x36 + Fyp * x79 + x158 * x193 * x208 + x180 * (x6 - 1) + x61 * (G * x139 * (x209 - 5) + Gp))) + S0_tyy * x183 + S0_yy * x179 + x106 * x175 + x110 * x171 * x2 - x110 * x192 + x162 * x177 * x208 - x164 * x195 + x164 * x207 - x168 * x194 + x168 * x206 - x169 * x174 + x170 * x180 * x184 + x171 * x82 * xi_y + x171 * xi_yy - x173 * xi_yy - 7 * x174 * x82 - x178 * x46 + x180 * x29 * x81 + x185 * x186 - x185 * x199 + x186 * x198 + x187 * x188 * xi_y - x187 * x19 * x199 - x197 * x199 * x205 * xi_y - x198 * x199) + x37 * (x113 * x160 * x161 + x127 ^ 2 * x88 + x127 * x140 * x2 * x63 + x140 * (u * x156 * x26 + x157 * x81 - x29 * (Fyph - u * (Fyh + x156 * x61) + x157 * x158)) + x155 * x63 ^ 2 - x161 * (12 * B * x162 + x139 * x26 * (-x165 * x81 + x167 * x29) + x163 * x164 - 6 * x167 * x29 * x81 + x38 * (Bs + x163 * x168 + x166 * x169 + 6 * x82 * (Bh + x161 * x170))) + x83 ^ 2 * x88)) + x149 * x37 * x52 * (-2 * Gc + u * x39 * x86 * (x141 * x26 - x141 * x68 + x29 * (Bh * x8 + Fyp + x141 * xi_y + x143 * x144 - x61)) + x122 * x39 * (x138 * x139 * x53 - 6 * x138 * x84 + x140 * (Gpt + x138 * x78 * xi_x - x78 * (Gt + x137 * x65))) - x129 * x76 - x142 * x83 * (x141 * x53 - x141 * x72 + x29 * (u * (Bt * u + Fx * x141 + Fx) + x141 * xi_x + x64)) + x145 * x148 * x43 - x39 * (2 * x130 + x133) * (S0_t * S0_y * x134 - x134 * x26 - x29 * (Gph + u * (Fy * Gpp - 3 * Gh - 3 * Gp * x61) + x134 * xi_y)) - x71 * (2 * Gpt + x128 * x76) - x76 * (2 * Gph + x129)) - 12 * x16 ^ 4 * x7 + x16 * x7 * (-x37 * x91 * (x150 * x83 + x35 * x86) + x58 * (-8 * S0_txy - S0_x * x93 - S0_y * x152 + x105 * (x152 + 8 * x65 + 8 * x74) + x148 * x43 * x44 - x15 * x151 - x151 * xi + x153 * x39 * (8 * x121 + 8 * x26 - 8 * x68) - x18 * xi_xy - x9 * (Sc + x154 * x71 + x76 * x95 + x76 * x96))) + x16 * (x37 * (4 * S0_tyy + S0_y * x93 - 4 * x105 * x71 - x113 * x115 - 4 * x122 * x39 * (S0_y * x119 + x118 * x27 + x34 * (Sph + x120 * x33 - x77 * (Sh + x116 * x61))) - x127 * x92 + x15 * x94 + x25 * xi_yy + x90 * (-Spp * x71 ^ 2 + Ss + x96 * (x33 + x67 * (x27 - 2 * x68) + x97)) + x94 * xi) + x58 * x83 * x92) + x18 * x20 * x5 - x35 ^ 2 * x41 + x35 * x40 * x59 * x7 * (16 * x28 * x53 + 8 * x54 + 8 * x57) + x43 * x44 * x52 * (x29 * x49 + x38 * xi ^ 2 + x42 * x50 * xi + x45 + x46 - x47 - x48) + x52 * x58 * x77 * (Fxph + Fxpp * x71 + Fypp * x76 + Fypt + u * x63 * x66 - u * (Fxh + Fxp * x71) - u * (Fyp * x76 + Fyt) - x71 * (Fxpp + x66 * x77) - x76 * (Fypp + x62 * x77) - x83 * x86 * x89) + (-x150 ^ 2 * x41 + x16 * (x210 * x58 * x86 + x37 * (4 * S0_txx + S0_x * x152 - x115 * (Fxt * x107 - S0_t * x47 + S0_txx * x29 + x111 * x214 - x147 * x72 + x212 + x213 * x42 * x8 + x42 * xi_xx + x53 * (-x146 + x147)) + x15 * x211 - 4 * x153 * x39 * (x130 + x73) + x210 * x215 + x211 * xi + x25 * xi_xx - 4 * x76 * (S0_tx * x102 + S0_x * x104 + x8 * (Spt + u * (Fx * Spp - Fx * x98 - 2 * St) + x100 * xi_x)) + x90 * (Sb - Spp * x76 ^ 2 + x154 * (x133 * x67 + x56 + 2 * x65)))) + x51 * (u * x37 * (Bb * x161 - Bp * Fxt * x149 - Bp * x214 * x88 - Fxp * u * x56 + Fxp * x0 * x146 * x67 - 2 * Fxpt + 6 * Fxt * x143 + Fxt * x77 - G * x247 * x254 + Gt ^ 2 * x88 + S0_tx * x43 * x77 * (u * x146 * (-S0_t * (x163 + 9 * x249 * x87) + x160) + x29 * (u * (Fx * (-x240 + x258 + 27 * x6 + 2) + x139 * (Bt + Gt * x36)) + x64 + xi_x * (x258 + 24 * x6))) + S0_txx * x101 * x159 * x161 + S0_xx * x240 * x50 * x67 + u ^ 9 * x214 * x252 + u * x66 ^ 2 - x114 * x240 * x256 + x132 * x245 * x67 + x149 * x215 * x66 * x67 - x161 * x223 + 6 * x19 * x237 * x67 * (x165 + 3 * x2 * x249) + 30 * x196 * x234 + x202 * x235 * x246 - x202 * x250 * x254 * xi_x - x209 * x216 * x67 + x209 * xi_xx + x213 * x244 + x213 * x251 + x213 * x90 + x215 ^ 2 * x89 - 48 * x221 * x254 + 4 * x228 + x229 * x248 + x229 * x253 + x233 * x247 - 2 * x240 * xi_xx + x241 * x242 - x243 * x254 + x243 * xi_x - x244 * x255 - x245 * x56 + x248 * x257 - x251 * x255 + x253 * x257 - x255 * x90 + 12 * x256 * x43 * x6) + x149 * x59 * (Bt * x149 * x217 + Fx * x184 * x222 - Fx * x200 * x232 + Fx * x205 * x220 - Fxp * Gt * x155 - Fxt * x172 + Fxt * x182 + Gb * x38 + S0_tx * (S0 * x146 * (S0 * x145 - x103 * x239) + x107 * (u * (30 * G * x65 + x143 * x236 + x222 * x8 - x224 - x226) + 12 * x136 * x238 * xi_x + x218 * (x6 + 1))) - S0_txx * x183 - S0_xx * x179 + x139 * x217 * xi_x - x143 * x219 * x72 - x171 * x228 - x171 * xi_xx + x173 * xi_xx - x175 * x212 + x176 * x216 + x178 * x45 + x188 * x222 * xi_x - 3 * x189 * x223 - x19 * x222 * x232 + x190 * x213 + x192 * x214 + x194 * x229 + x195 * x230 + x197 * x213 * x231 * x38 + x204 * x234 - 36 * x205 * x232 * x233 + x206 * x229 + x207 * x230 + 7 * x217 * x85 - x219 * x84 + x220 * x221 - x224 * x225 + x224 * x227 - x225 * x226 + x226 * x227 - 24 * x227 * x233 - x231 * x232 * x236 + x237 * x239 * x34))) * exp(2 * x6))
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
x18 = x15 * x4 * (-Gp + x17)
x19 = S0 + 2 * x7
x20 = x0 * x11 * (B * x16 - Bp)
x21 = S0 ^ 3
x22 = 1 / x21
x23 = sech(x14)
x24 = S0 * u
x25 = x23 * x24
x26 = u ^ 4
x27 = S0 * x26
x28 = 2 * S * x27
x29 = S0 ^ 2
x30 = 2 * x8
x31 = x29 * (u * xi + 1) + x5 * (S0 - x30)
x32 = S * x4
x33 = 2 * xi_y
x34 = x29 * x4
x35 = S0_ty * x28 + S0_y * x31 + x34 * (-Fy * S0 + Fy * x30 + Sh * u + x32 * x33)
x36 = -Fyp
x37 = Fy * u
x38 = -x36 - x37
x39 = x35 * x38
x40 = x36 + 2 * x37
x41 = -x40
x42 = 2 * x5
x43 = S0_ty * x24
x44 = 5 * x4
x45 = x10 * (S0_y * x41 * x42 + x29 * (-2 * Fyph - u * x33 * x41 + u * (Fy ^ 2 * x44 + 2 * Fyh + Fyp ^ 2 - 4 * Fyp * x37)) + 2 * x40 * x43)
x46 = 2 * x25
x47 = S0 ^ 4
x48 = 2 * x21
x49 = Sd * x0 * x48 + u * x48 * (S0_t + x6) + x4 * (-S0 * S0_xx - S0 * S0_yy + S0_t ^ 2 * x29 + S0_t * x48 * xi + S0_x ^ 2 + S0_y ^ 2 + x47 * xi ^ 2) + x47
x50 = x13 * x2
x51 = Fx * u
x52 = 2 * xi_x
x53 = S0_tx * x28 + S0_x * x31 + x34 * (-Fx * S0 + Fx * x30 + St * u + x32 * x52)
x54 = -Fxp
x55 = 2 * x51 + x54
x56 = -x55
x57 = S0_tx * x24
x58 = (x10 * (S0_x * x42 * x56 + x29 * (-2 * Fxpt - u * x52 * x56 + u * (Fx ^ 2 * x44 + Fxp ^ 2 - 4 * Fxp * x51 + 2 * Fxt)) + 2 * x55 * x57) + x53 * (4 * Fxp - 4 * x51)) * exp(2 * x1)
x59 = Gh * x29
x60 = Gp * x29
x61 = Fx * Gp
x62 = S0_y * x5
x63 = Gt * x29
x64 = 3 * G
x65 = Fxp * x64
x66 = Fy * Gp
x67 = S0_x * x5
x68 = Fyp * x64
x69 = x17 * x29
x70 = u * x29
x71 = Fxp * xi_y
x72 = Fyp * xi_x
x73 = Fx * x34
x74 = Fy * x34
x75 = 4 * x2
x76 = 2 * x2
x77 = 4 * u
x78 = sinh(x14)
x79 = x24 * x78
x80 = cosh(x14)
x81 = Fxp * x70
x82 = 2 * x4
x83 = S0_t * S0_y
x84 = Fx * x83
x85 = S0_t * S0_x
x86 = Fy * x85
x87 = Bp * x26
x88 = x0 * x29
x89 = Bp * x88
x90 = x26 * x29
x91 = Fy * x90
x92 = 3 * B
x93 = Fxp * x92
x94 = Fyp * x92
x95 = x90 * x92
x96 = Fx * x90
x97 = u ^ 5 * x29
x98 = 3 * x1
x99 = Bp * x4
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x12 * x3
BB[1,2] = 0
CC[1,1] = -u * x13 * x3 * (-Sp * x4 + x18 * x9 + x19 + x5 * (x18 + 2) + 4 * x8)
CC[1,2] = -x15 * x20 * x3
SS[1] = x22 * (x10 * x23 * x27 * x75 * (-Fxh * x60 + Fxh * x69 + Fxp * x59 - Fyp * x63 + Fyt * x60 - Fyt * x69 + x37 * x60 * xi_x + x37 * x63 + x43 * (3 * Fxp * G - x61) - x51 * x59 + x57 * (x66 - x68) + x61 * x62 - x61 * x70 * xi_y - x62 * x65 + x65 * x74 - x66 * x67 + x67 * x68 - x68 * x73 + x69 * x71 - x69 * x72) + 8 * x25 * x39 + x45 * x46 - x46 * x58 + x49 * x50 * (-12 * B * u + 4 * Bp))
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x12 * x75
CC[2,1] = x20 * x76 * sinh(2 * x14)
CC[2,2] = -x50 * x77 * (x19 + x4 * (S * x77 - Sp) + x42)
SS[2] = x22 * (4 * S0 * u * x2 * x80 * (x35 * (-x51 - x54) + x38 * x53) - S0 * x39 * x77 * x78 - x10 * x24 * x76 * x80 * (Bh * Fxp * x88 - Bh * x96 + Bp * x91 * xi_x - Bp * x96 * xi_y - Bt * Fyp * x88 + Bt * x91 - 5 * Fx * Fy * x88 - Fx * x94 * x97 - Fxh * x70 - Fxh * x89 + Fxh * x95 - Fxp * x62 + 2 * Fxp * x74 + Fxph * x29 + Fy * x93 * x97 - Fyp * x67 + 2 * Fyp * x73 - Fyp * x81 + Fypt * x29 - Fyt * x70 + Fyt * x89 - Fyt * x95 - x26 * x83 * x93 + x26 * x85 * x94 - x33 * x73 + x43 * (Fxp * x98 + Fxp - x51 * (x99 + 2)) - x52 * x74 + x57 * (-Fyp * x98 + Fyp + x37 * (x99 - 2)) + x70 * x72 + x71 * x95 - x72 * x95 + x81 * xi_y + x82 * x84 + x82 * x86 + x84 * x87 - x86 * x87) + x13 * x2 * x49 * (-6 * G * u + 2 * Gp) - x45 * x79 - x58 * x79)
    
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

x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = S0_t * u
x4 = u * xi
x5 = S * x0
x6 = S0 * x4 + S0 + x3 + x5
x7 = x6 ^ 4
x8 = x2 * x7
x9 = 2 * x0
x10 = u ^ 2
x11 = 4 * u
x12 = S * x10
x13 = 1 / u
x14 = S0 * (x13 + xi)
x15 = S0_t + x12 + x14
x16 = Fy * S0
x17 = 2 * x5
x18 = 2 * xi_y
x19 = Fy * x17 + Sh * u + x12 * x18 - x16
x20 = S0 ^ 2
x21 = x10 * x20
x22 = x20 * (x4 + 1) + x3 * (S0 - x17)
x23 = S0 * S0_ty
x24 = 2 * x23
x25 = u ^ 4
x26 = S * x25
x27 = S0_y * x22 + x24 * x26
x28 = x19 * x21 + x27
x29 = G * x0
x30 = cosh(x29)
x31 = 1 / x10
x32 = S0 ^ 4
x33 = 1 / x32
x34 = 4 * x33
x35 = x30 * x31 * x34
x36 = sinh(x29)
x37 = 3 * u
x38 = G * x37
x39 = 3 * x3
x40 = S0_y * x39
x41 = G * x10
x42 = 3 * Fy
x43 = -G * x40 + x20 * (Gh + x38 * xi_y + x41 * x42) + x23 * x38
x44 = 4 * x10
x45 = x33 * x44
x46 = 4 * S0
x47 = 4 * S0_yy
x48 = 1 / x20
x49 = S0_t * S0_y
x50 = x23 - x49
x51 = Fy * u
x52 = x51 + xi_y
x53 = x48 * x50 + x52
x54 = Spp * x53
x55 = Sph + x54
x56 = 2 * x51
x57 = x24 - 2 * x49
x58 = 2 * u
x59 = Sp * x58 - Spp
x60 = -x59
x61 = x10 * x60 / S0
x62 = S0_t * x10
x63 = x48 * (x20 + x59 * x62)
x64 = S0_ty * x61 + S0_y * x63 + x10 * (Sph + u * (Fy * Spp - 2 * Sh - Sp * x56) + x60 * xi_y)
x65 = S0 ^ 3
x66 = 1 / x65
x67 = -Sp
x68 = x10 * (-S * x58 - x67)
x69 = S0 + x68
x70 = S0 * S0_yy
x71 = S0_y ^ 2
x72 = 2 * S0_t
x73 = x71 * x72
x74 = u * x65
x75 = x10 * x65
x76 = Fy * xi_y
x77 = x10 * x16
x78 = Fy ^ 2
x79 = x0 * x65
x80 = 2 * S0_y
x81 = x77 - x80
x82 = Fyh * x74 - S0_t * x70 + S0_tyy * x20 + x23 * x81 - x49 * x77 + x65 * xi_yy + x73 + x75 * x76 + x78 * x79
x83 = x20 * x52
x84 = x50 + x83
x85 = S * x37 + x67
x86 = -x85
x87 = x0 * x86
x88 = x20 - x72 * x87
x89 = u * x86
x90 = B * x37
x91 = B * x10
x92 = -B * x40 + x20 * (Bh + x42 * x91 + x90 * xi_y) + x23 * x90
x93 = x6 ^ 2
x94 = Fyp - x51
x95 = x57 + 2 * x83
x96 = -Bp * x37 + Bpp
x97 = 3 * Bh
x98 = 3 * Bp
x99 = -Bp
x100 = B * x11 + x99
x101 = -x100
x102 = x101 * x37
x103 = -x90 - x99
x104 = S0 * x103
x105 = x33 * x58
x106 = Gp * x32
x107 = u * x106
x108 = x32 * x38
x109 = 6 * Gh
x110 = x109 * x32
x111 = u * x110
x112 = Gp * S0
x113 = Gp * x20
x114 = S0_t * x113
x115 = 6 * G
x116 = S0 * x115 * x3
x117 = G * x20
x118 = x117 * x39
x119 = x20 * x3
x120 = 2 * x32
x121 = x0 * x120
x122 = x10 * x110
x123 = 3 * x32 * x41
x124 = Gp - x38
x125 = x124 * x65
x126 = x0 * x106
x127 = u ^ 5
x128 = Fy * x127
x129 = x25 * x32
x130 = B * x129
x131 = 6 * Bh
x132 = G * x32
x133 = x128 * x132
x134 = G * x129
x135 = x134 * xi_y
x136 = 27 * x29
x137 = x32 * x76
x138 = 15 * x134
x139 = xi_y ^ 2
x140 = 12 * x41
x141 = x140 * x32
x142 = S0_t ^ 2
x143 = x140 * x142
x144 = u ^ 6
x145 = 36 * B * x144
x146 = G * x145
x147 = x20 * x49
x148 = B * x25
x149 = x117 * x25
x150 = 27 * Fy
x151 = x150 * x29
x152 = x21 * x49
x153 = Fy * Gp
x154 = 18 * B
x155 = x132 * x154
x156 = u ^ 7 * x155
x157 = x127 * x155
x158 = G * x127
x159 = x142 * x154 * x158
x160 = x117 * x49
x161 = B * x127
x162 = 36 * x161
x163 = 3 * x1
x164 = x163 - 2
x165 = 2 * S0 * x124
x166 = G * x62
x167 = S0 * S0_y
x168 = G * u
x169 = 12 * x168
x170 = 4 * x1
x171 = x105 * x36
x172 = S0 * S0_tx
x173 = 2 * x172
x174 = S0_x * x22
x175 = Fx * S0
x176 = 2 * xi_x
x177 = Fx * x17 + St * u + x12 * x176 - x175
x178 = x173 * x26 + x174 + x177 * x21
x179 = S0_x * x39
x180 = 3 * Fx
x181 = -G * x179 + x172 * x38 + x20 * (Gt + x180 * x41 + x38 * xi_x)
x182 = x178 * x45
x183 = 4 * S0_xx
x184 = S0_t * S0_x
x185 = x172 - x184
x186 = x185 * x48
x187 = Fx * u
x188 = x187 + xi_x
x189 = x186 + x188
x190 = Spp * x189 + Spt
x191 = x173 - 2 * x184
x192 = 2 * x187
x193 = x176 + x192
x194 = S0 * S0_xx
x195 = S0_x ^ 2
x196 = x195 * x72
x197 = Fx * x75
x198 = x10 * x175
x199 = Fx ^ 2
x200 = 2 * S0_x
x201 = Fxt * x74 - S0_t * x194 + S0_txx * x20 + x172 * (x198 - x200) - x184 * x198 + x196 + x197 * xi_x + x199 * x79 + x65 * xi_xx
x202 = x188 * x20
x203 = x185 + x202
x204 = S0_x * x88 + x173 * x87 + x21 * (Spt + x176 * x89 - x58 * (St + x187 * x85))
x205 = -B * x179 + x172 * x90 + x20 * (Bt + x180 * x91 + x90 * xi_x)
x206 = -Fxp + x187
x207 = 3 * Bt
x208 = 6 * Gt
x209 = x208 * x32
x210 = u * x209
x211 = x10 * x209
x212 = x130 * xi_x
x213 = 6 * x32
x214 = x32 * xi_x
x215 = Fx * x136
x216 = Fx * Gp
x217 = x184 * x21
x218 = xi_x ^ 2
x219 = x184 * x20
x220 = G * x21
x221 = x117 * x184
x222 = x163 + 2
x223 = 6 * x222
x224 = 6 * x41
x225 = x30 * x33
x226 = 2 * S0_txy
x227 = 2 * xi_xy
x228 = 2 * S0_xy
x229 = 2 * x10
x230 = S0_xy * x72
x231 = Fx * S0_ty
x232 = Fx * Fy
x233 = S0_ty * x65
x234 = 4 * Gpp
x235 = x214 * xi_y
x236 = S0_x * x142
x237 = x187 * x234
x238 = x32 * xi_y
x239 = Fy * x234
x240 = 24 * x41
ABCS[1] = x8 * x9
ABCS[2] = 8 * x10 * x8
ABCS[3] = x11 * x8

 
if u == 0.0 
	ABCS[4] = -4*S0^3*Sp
else
	ABCS[4] = x10 * (4 * x15 ^ 4 * x2 + x15 * (-x28 * x36 * x43 * x45 + x30 * (-4 * S0_tyy - 8 * S0_y * xi_y + 4 * x10 * x28 * x33 * x92 - x13 * x47 + 4 * x33 * x84 * (S0_y * x88 + x21 * (Sph + x18 * x89 - x58 * (Sh + x51 * x85)) + x24 * x87) - x44 * (-Spp * x53 ^ 2 + Ss + x55 * (x18 + x48 * x57 + x56)) - x46 * xi_yy - x47 * xi + 4 * x53 * x64 + 4 * x66 * x69 * x82)) + x2 * (x15 ^ 2 * (-x225 * x9 * (3 * B * Fy * Gt * x127 * x32 + 3 * B * Gh * S0_t * S0_x * x20 * x25 + 3 * B * Gt * S0_ty * x25 * x65 - 3 * B * Gt * x147 * x25 + 3 * B * Gt * x25 * x32 * xi_y + 3 * Bh * Fx * G * x127 * x32 + 3 * Bh * G * x25 * x32 * xi_x + Bh * Gt * x0 * x32 + 3 * Bt * G * S0_t * S0_y * x20 * x25 - Bt * Gh * x0 * x32 + 2 * Fx * Fy * Gp * x0 * x32 + 27 * Fx * G * S0_t * S0_y * x0 * x20 + Fx * Gp * S0_ty * x10 * x65 + Fx * Gp * x10 * x32 * xi_y + 4 * Fx * Gpp * S0_t * S0_y * u * x20 - Fx * x10 * x239 * x32 - Fx * x122 + Fxh * Gp * u * x32 - Fxh * x123 + 27 * Fy * G * S0_t * S0_x * x0 * x20 + Fy * Gp * x10 * x32 * xi_x + 4 * Fy * Gpp * S0_t * S0_x * u * x20 - Fy * x211 + Fyt * Gp * u * x32 - Fyt * x123 + 24 * G * S0 * S0_t * S0_ty * S0_x * x10 + 24 * G * S0_t * S0_x * x10 * x20 * xi_y + 6 * G * S0_t * S0_xy * u * x20 + 24 * G * S0_t * S0_y * x10 * x20 * xi_x + 6 * G * S0_ty * S0_x * u * x20 - 12 * G * S0_x * x167 * x3 - G * x207 * x233 * x25 - Gc * x120 + 6 * Gh * S0_t * S0_x * u * x20 - Gh * x161 * x180 * x32 - 3 * Gh * x212 + 4 * Gp * S0 * S0_t * S0_x * S0_y + 2 * Gp * x32 * xi_xy + 4 * Gpp * S0 * S0_t * S0_ty * S0_x + 4 * Gpp * S0_t * S0_x * x20 * xi_y + 4 * Gpp * S0_t * S0_y * x20 * xi_x + 6 * Gt * S0_t * S0_y * u * x20 + 2 * S0_txy * x124 * x65 - S0_ty * x113 * x200 - S0_ty * x208 * x74 - S0_y * x234 * x236 - S0_y * x236 * x240 - x111 * xi_x - x113 * x230 - x133 * x207 - 30 * x134 * x232 - x135 * x207 - x136 * x231 * x65 - x149 * x184 * x97 - x151 * x214 - x152 * x216 - x153 * x217 - x168 * x213 * xi_xy - x172 * (S0_y * (x165 - x72 * (2 * Gpp + x140)) + x20 * (u * (Gh * x163 - Gp * x51 + x109 + x150 * x41 + x239 - x29 * x97) + xi_y * (x234 + x240)) + 4 * x23 * (Gpp + x224)) - x210 * xi_y - x214 * x234 * x51 - x215 * x238 - x233 * x234 * xi_x - x233 * x237 - x233 * x240 * xi_x - x234 * x235 - x235 * x240 - x237 * x238) + x36 * (x144 * x181 * x34 * x43 + x206 * x229 * x94) + (-8 * S0 - 8 * x68) * (S0 * x31 / 2 + Sd * u + x13 * (S0 * xi + S0_t) + x66 * (x142 * x20 - x194 + x195 + x32 * xi ^ 2 + x65 * x72 * xi - x70 + x71) / 2)) - x31 * x33 * x36 * (-x10 * x19 * x20 - x27) * (-8 * x10 * x177 * x20 - 16 * x172 * x26 - 8 * x174) - x58 * x7 * (-Bd * x103 * x30 ^ 2 - Gd * x124) + (x10 * x225 * (x178 * x43 + x181 * x28) + x36 * (S0 * x227 + S0_x * x18 + x13 * x228 - x204 * x33 * x95 + x226 + x228 * xi + x229 * (Sc + x189 * x54 + x189 * x55 + x190 * x53) - x64 * (2 * x186 + x193) - x66 * x69 * (Fxh * x74 + Fy * x75 * xi_x + Fyt * x74 - S0 * x230 + 4 * S0_x * x49 + x172 * x81 - x184 * x77 + x197 * xi_y - x198 * x49 + x20 * x226 - x200 * x23 + x21 * x231 + x227 * x65 + 2 * x232 * x79) + x80 * xi_x)) * (4 * S0_t + 4 * x12 + 4 * x14)) + x28 ^ 2 * x35 + x93 * (x171 * (B * x110 * x128 + Bh * Gh * x121 - Fy * x122 - Fy * x145 * x160 + Fyh * x107 - Fyh * x123 + 24 * G * x152 * xi_y - Gs * x32 + S0_ty ^ 2 * x115 * x164 * x21 + S0_ty * (x167 * (-12 * x164 * x166 - x165) + x74 * (x109 * (x1 - 1) + x131 * x29 + x164 * x169 * xi_y + x51 * (Gp + 9 * x168 * (x170 - 3)))) + S0_tyy * x125 + S0_y * x109 * x119 - S0_yy * x114 + S0_yy * x118 + x10 * x106 * x76 + x106 * xi_yy - x108 * xi_yy + x109 * x130 * xi_y - x109 * x147 * x148 - x111 * xi_y + x112 * x73 - x116 * x71 + x126 * x78 + x131 * x133 + x131 * x135 - x131 * x149 * x49 - x136 * x137 + x137 * x146 - x138 * x78 - x139 * x141 + x139 * x157 - x143 * x71 + x147 * x151 - x152 * x153 + x156 * x78 + x159 * x71 - x160 * x162 * xi_y) + x30 * (-x105 * (Bpp * x84 ^ 2 - Bs * x32 + x0 * x43 ^ 2 + x0 * x92 ^ 2 + x104 * x82 + x84 * (-x101 * x40 + x102 * x23 + x20 * (Bph + x102 * xi_y - x37 * (Bh + x100 * x51))) + x84 * (x20 * (Bph + u * (Bpp * Fy - x51 * x98 - x97) + x96 * xi_y) + x23 * x96 - x49 * x96) - x95 * (Bpp * x23 - Bpp * x49 + x20 * (Bph + Bpp * x51 + Bpp * xi_y))) + x94 ^ 2)) + (x15 * (-x181 * x182 * x36 + x30 * (-4 * S0_txx - 8 * S0_x * xi_x - x13 * x183 - x182 * x205 - x183 * xi + 4 * x189 * (S0_tx * x61 + S0_x * x63 + x10 * (Spt + u * (Fx * Spp - Sp * x192 - 2 * St) + x60 * xi_x)) + 4 * x201 * x66 * x69 + 4 * x203 * x204 * x33 - x44 * (Sb - Spp * x189 ^ 2 + x190 * (x191 * x48 + x193)) - x46 * xi_xx)) + x178 ^ 2 * x35 + x93 * (-x171 * (Bt * Fx * x158 * x213 + Bt * Gt * x121 - Bt * x115 * x219 * x25 + 6 * Bt * x134 * xi_x - Fx * x145 * x221 + Fx * x146 * x214 + Fx * x161 * x209 + Fx * x211 - Fxt * x107 + Fxt * x123 + Gb * x32 + S0_tx ^ 2 * x220 * x223 + S0_tx * (S0 * x200 * (S0 * x124 - x166 * x223) + x74 * (u * (Bt * x224 + 9 * G * x187 * (x170 + 3) - x216) + x169 * x222 * xi_x + x208 * (x1 + 1))) - S0_txx * x125 - S0_x * x119 * x208 + S0_xx * x114 - S0_xx * x118 - x10 * x214 * x216 - x106 * xi_xx + x108 * xi_xx - x112 * x196 + x116 * x195 - x126 * x199 + x138 * x199 + x141 * x218 + x143 * x195 - x148 * x208 * x219 + x156 * x199 + x157 * x218 + x159 * x195 - x162 * x221 * xi_x - 24 * x184 * x220 * xi_x + x208 * x212 + x210 * xi_x + x214 * x215 - x215 * x219 + x216 * x217) + x30 * (-x105 * (Bb * x32 - Bpp * x203 ^ 2 + x0 * x181 ^ 2 + x0 * x205 ^ 2 - x104 * x201 - x203 * (-x101 * x179 + x102 * x172 + x20 * (Bpt + x102 * xi_x - x37 * (Bt + x100 * x187))) - x203 * (x172 * x96 - x184 * x96 + x20 * (Bpt + u * (Bpp * Fx - x187 * x98 - x207) + x96 * xi_x)) + (x191 + 2 * x202) * (Bpp * x172 - Bpp * x184 + x20 * (Bpp * x187 + Bpp * xi_x + Bpt))) + x206 ^ 2))) * exp(2 * x1))
	
end

    nothing
    
end
