
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
x24 = u * x2
x25 = x10 * x24
x26 = x25 * (x13 + x15 + x23 * x3 + 8)
x27 = 4 * S0
x28 = Spp * x27
x29 = S0 ^ 2
x30 = Bpp * x29
x31 = Sp * u
x32 = 24 * S0
x33 = x29 * xi
x34 = u * x29
x35 = 5 * Bp
x36 = S * x3
x37 = 4 * Spp
x38 = S ^ 2
x39 = u ^ 5
x40 = 12 * x39
x41 = 4 * x0
x42 = u ^ 6
x43 = 36 * S0
x44 = B * x39
x45 = S * x44
x46 = u ^ 7
x47 = S * x46
x48 = 6 * B
x49 = x47 * x48
x50 = u ^ 4
x51 = Sp * x50
x52 = x48 * x51
x53 = S * x50
x54 = Bp * x53
x55 = S * x42
x56 = 2 * Sp
x57 = Bp * x56
x58 = S0 * x0
x59 = 2 * S0
x60 = Bpp * x7
x61 = 2 * x5
x62 = S0 * xi
x63 = Sp * x3
x64 = u ^ 8
x65 = B * x64
x66 = B * x29
x67 = 9 * x3
x68 = Bp * x38
x69 = xi ^ 2
x70 = B * Bp
x71 = 12 * S0
x72 = x47 * x71
x73 = B * x55
x74 = x39 * x62
x75 = Bp * x39
x76 = S * x75
x77 = x50 * x62
x78 = Bpp * x53
x79 = x59 * xi
x80 = G * Gp
x81 = Bp ^ 2
x82 = u ^ 9
x83 = x38 * x82
x84 = x0 * x29
x85 = u ^ 10
x86 = x68 * x85
x87 = Bp * x48
x88 = x29 * x50
x89 = 24 * x33
x90 = B ^ 2
x91 = 18 * S0
x92 = S * x64
x93 = 12 * x80
x94 = x38 * x85
x95 = G ^ 2
x96 = Gp ^ 2
x97 = x27 * x96
x98 = x71 * xi
x99 = S * x65
x100 = 9 * x90
x101 = u ^ 11 * x38
x102 = x29 * x39
x103 = 18 * x95
x104 = 2 * x96
x105 = x33 * x40
x106 = 18 * x90
x107 = S * x82
x108 = x107 * x62
x109 = x47 * x81
x110 = 36 * x95
x111 = 15 * x50
x112 = x66 * x69
x113 = x33 * x42
x114 = 7 * Bp
x115 = x50 * x81
x116 = x33 * x50
x117 = x39 * x81
x118 = x29 * x69
x119 = x42 * x87
x120 = x118 * x42
x121 = x118 * x46
x122 = sinh(x21)
x123 = G * S
x124 = x123 * x65
x125 = x122 * x72
x126 = B * Gp
x127 = Bp * G
x128 = Gp * x122
x129 = Bp * x27
x130 = G * x122
x131 = 18 * x130
x132 = Gp * x48
x133 = x122 * x132
x134 = 6 * x127
x135 = x122 * x134
x136 = 2 * Gp
x137 = Bp * x136
x138 = x122 * x137
x139 = B * x62
x140 = 36 * x130
x141 = x122 * xi
x142 = x141 * x71
x143 = x105 * x122
x144 = Bpp * u
x145 = 15 * x1
x146 = 7 * x14
x147 = x48 * x75
x148 = x40 * x80
x149 = x100 * x42
x150 = x103 * x42
x151 = -x19
x152 = x122 * x151
x153 = -Gp
x154 = G * x17
x155 = x153 + x154
x156 = -2 * x155
x157 = x152 * x156 * x50
x158 = x3 * x81
x159 = 6 * x1
x160 = x13 + 5
x161 = Bpp + u * (-Bp * (x159 + 7) + x158 + x160 * x18)
x162 = u * x22
x163 = Bpp * x8
x164 = x158 * x8
x165 = x5 + 1
x166 = x159 * x165 + 7 * x5
x167 = -x56
x168 = S * u
x169 = x13 * x165 + 5 * x5
x170 = x22 * x8
x171 = 2 * u
x172 = S0 * x115
x173 = -Sp
x174 = x13 + 7
x175 = 6 * S0
x176 = x59 * x96
x177 = S0 * x149
x178 = Bpp * S0
x179 = x3 * xi
x180 = x46 * x62
x181 = x27 * x5 + x27
x182 = x24 * (B * x101 * x131 + B * x113 * x140 + 4 * Bp * x116 * x128 - Bp * x123 * x142 * x64 - Bp * x98 * x99 + Bpp * x38 * x42 + Gp * x129 * x141 * x47 - Gp * x142 * x99 + 48 * S0 * x36 - S0 * x52 - 16 * S0 * x54 + S0_t ^ 2 * u * (x104 * x50 + x115 + x144 + x145 - x146 - x147 - x148 + x149 + x150 + x157 + x161 * x162 + 4) + 2 * S0_t * (Bp * x51 + S0 * x144 - S0 * x147 - S0 * x148 + S0 * x150 - 3 * Sp * x44 + Spp * x171 - x0 * x114 * x62 + x1 * x71 + x100 * x107 + x100 * x180 + x103 * x107 + x103 * x180 + x104 * x47 + x109 + x111 * x139 + x117 * x62 - x119 * x62 - x14 * x175 + x157 * x8 + x162 * (u * (-Bp * (S0 * (x166 + 6) + x3 * (9 * x168 + x173 + x48 * x53)) + x164 + x18 * (S0 * (x169 + 4) + x3 * (x168 * x174 + x173))) + x163) + x172 + x176 * x39 * xi + x176 * x50 + x177 + x178 * x179 + x181 - x42 * x80 * x98 - 8 * x63 + 16 * x7 + 21 * x73 - 9 * x76 + x78 - x87 * x92 - x92 * x93) - Sp ^ 2 * x41 - Sp * x48 * x74 - Sp * x49 + x1 * x89 + x100 * x101 + x100 * x102 + x100 * x121 + x101 * x103 + x102 * x103 + x102 * x104 * x69 + x103 * x121 + x104 * x83 + x104 * x84 - x105 * x70 + x106 * x108 + x106 * x113 + x107 * x139 * x140 + x108 * x110 + x109 * x79 + x110 * x113 + x111 * x112 + x112 * x131 * x46 - x114 * x69 * x84 + 2 * x115 * x33 + 4 * x116 * x96 + x117 * x118 - x118 * x119 + x118 * x122 * x136 * x75 - x120 * x133 - x120 * x135 - x120 * x93 + x122 * x124 * x43 - x125 * x126 - x125 * x127 - x126 * x143 - x127 * x143 + x128 * x129 * x55 - 6 * x130 * x86 + x131 * x29 * x44 - x133 * x88 - x133 * x94 - x135 * x88 + x138 * x83 + x138 * x84 - 12 * x14 * x33 + x170 * (u * (-Bp * (S0 * (x166 + 5) + x3 * (x167 + x168 * (x159 + 11))) + x164 + x18 * (S0 * (x169 + 3) + x3 * (S * x17 * (x1 + 3) + x167))) + x163) + x28 * x5 + x28 + x3 * x30 * x69 + x30 * x61 + x30 - x31 * x32 - x32 * x47 * x80 + 8 * x33 - x34 * x35 + 4 * x34 * x69 + x37 * x7 + x38 * x40 + 27 * x38 * x65 - x39 * x80 * x89 + x43 * x45 + x43 * x92 * x95 - 11 * x46 * x68 + x47 * x97 * xi - x48 * x86 + x55 * x57 + x55 * x59 * x81 + x55 * x97 + x57 * x58 + x57 * x77 + x59 * x60 - 16 * x62 * x63 + 32 * x62 * x7 + 42 * x62 * x73 - 18 * x62 * x76 - 24 * x62 * x80 * x92 + x66 * x67 - x70 * x72 + x78 * x79 + x81 * x83 + x81 * x84 - x87 * x88 - x88 * x93 + x90 * x91 * x92 - x93 * x94)
x183 = 6 * u
x184 = x10 * x183
x185 = S * x171 + x173
x186 = -x185 * x3
x187 = S0 + x186
x188 = 4 * x187
x189 = 4 * x2
x190 = x189 * x9
x191 = 1 / u
x192 = x187 * x189
x193 = x191 * x192
x194 = 1 / x3
x195 = S0_x * x194
x196 = x24 * x9
x197 = x10 * x122
x198 = x187 * x191
x199 = x2 * xi_x
x200 = x187 ^ 2 * x191
x201 = x189 * x200
x202 = x17 * x197
x203 = x0 * x11
x204 = x155 * x203
x205 = cosh(x20) ^ 2
x206 = x11 * x205
x207 = Bp_x * x2
x208 = 2 * Gpp
x209 = G * x171 + x153
x210 = 12 * u
x211 = x10 * (x208 + x209 * x210)
x212 = 6 * x25
x213 = x190 * (Spp - 4 * x31 + 6 * x36)
x214 = x0 * x10
x215 = x151 * x22
x216 = G_y * x203
x217 = 1 / x29
x218 = x217 * (S0 * S0_tx - S0_t * S0_x)
x219 = x197 * xi_y
x220 = x0 * x151 ^ 2
x221 = x204 * xi_y
x222 = B * x171 + x16
x223 = Bpp + x183 * x222
x224 = x19 * x205
x225 = S0_x * x190 * x224
x226 = x196 * x224
x227 = 4 * x226
x228 = G * x183
x229 = -x136 - x152 + x228
x230 = x229 * x9
x231 = x171 * x230
x232 = S0_y * x229
x233 = x205 * x9
x234 = x19 * x233
x235 = x234 * x41
x236 = 2 * x0
x237 = x217 * (S0 * S0_ty - S0_t * S0_y)
x238 = x0 * x19
x239 = x206 * x238
x240 = B_x * x2
x241 = x199 * x206
x242 = x197 * x220
x243 = x204 * x237
x244 = x197 * x223
x245 = x11 * x2
x246 = x0 * x245
x247 = x246 * (-Gp - x152 + x154)
x248 = x206 * x220
x249 = x2 * x218
x250 = -S0 + x185 * x3
x251 = x250 * xi_x
x252 = x155 * x247
x253 = x218 * x250
x254 = x136 - x152 - x228
x255 = x2 * x254
x256 = x214 * x255
x257 = Gpp * x59
x258 = S0 * u
x259 = G * x91
x260 = Gp * S0
x261 = S0 * x122
x262 = x122 * x81
x263 = x122 * x62
x264 = x0 * x156
x265 = x196 * (B * x261 * x67 + Bpp * x122 * x6 - 18 * G * x139 * x42 + G * x175 * x75 * xi + 12 * G * x51 - Gp * Sp * x41 + 10 * Gp * x258 + Gp * x49 + 22 * Gp * x53 + x0 * x122 * x57 + x100 * x122 * x92 + x100 * x261 * x39 + x122 * x178 - x122 * x258 * x35 + 27 * x122 * x45 - x122 * x47 * x87 - x122 * x52 - 11 * x122 * x54 + x122 * x60 - 54 * x123 * x39 - 18 * x124 + x127 * x175 * x50 + x132 * x74 + x134 * x47 - x137 * x55 - x137 * x58 - x137 * x77 + x141 * x172 + x141 * x177 + x145 * x263 - x146 * x263 - x147 * x263 + x151 * x170 * x264 + 14 * x179 * x260 - 30 * x20 * x62 - x208 * x7 - x257 * x5 - x257 - x259 * x3 - x259 * x44 + x260 * x48 * x50 - x261 * x50 * x87 + x262 * x55 + x262 * x58 - x4 * (-x122 * x161 + x171 * (Gp * (-x15 - x174) + x154 * (x15 + x160)) + x208 - x215 * x264))
x266 = S_y * u
x267 = x27 * xi_y
x268 = 4 * x186 + x27
x269 = x191 * x268
x270 = S0_y * x194
x271 = 4 * xi_y
x272 = x181 + 4 * x4 + 4 * x7
x273 = x9 * (-16 * x31 + 24 * x36 + x37)
x274 = S0_y * x233
x275 = u * x234
x276 = x155 * x246
x277 = x245 * (Gpp + x183 * x209)
x278 = x10 * x205 * (2 * Bpp + x210 * x222)
x279 = x152 + x155
x280 = x19 * x276
x281 = x255 * x9
x282 = S0_x * x281
x283 = x171 * x281
x284 = x23 * x276
AA[1,1] = x12 * x2
AA[1,2] = 0
BB[1,1] = x26
BB[1,2] = x26
CC[1,1] = x182
CC[1,2] = x182

if u==0.0
	println("entered fx u zero")
	SS[1] = -2*S0^2*(Bp_x) + 2*S0^2*(Gp_y) + 2*Bpp*S0*(S0_tx) + 4*Spp*(S0_tx) - 2*Gpp*S0*(S0_ty) - 4*Bp*S0*(S0_x) + 4*Sp*(S0_x) - 2*Bpp*(S0_t)*(S0_x) - (4*Spp*(S0_t)*(S0_x))/S0 + 4*Gp*S0*(S0_y) + 2*Gpp*(S0_t)*(S0_y) - 4*S0*(Sp_x) + 2*Bpp*S0^2*(xi_x) + 4*S0*Spp*(xi_x) - 2*Gpp*S0^2*(xi_y)
else
	SS[1] = B_x * x205 * x212 - B_y * x152 * x214 - B_y * x202 + B_y * x204 + Bp_y * x197 + G_x * x247 - G_y * x184 + Gp_y * x11 + S0_tx * x193 + S0_tx * x227 - S0_ty * x231 + S0_x * x193 * xi + S_x * x188 * x24 + 8 * S_x * x196 + S_x * x2 * x235 - S_y * x230 * x236 - Sp_x * x190 - u * x230 * x59 * xi_y + x187 * x231 * x237 + x19 * x221 + x19 * x243 - x190 * x195 + x192 * x195 + x198 * x199 * x27 - x201 * x218 - x201 * xi_x - x206 * x207 + x206 * x223 * x249 - x211 * x237 - x211 * xi_y + x213 * x218 + x213 * xi_x + x215 * x216 + x215 * x221 + x218 * x252 + x219 * x220 - x219 * x223 + x220 * x241 + x223 * x241 + x225 * x5 + x225 + x226 * x27 * xi_x + x227 * x251 + x227 * x253 - x23 * x243 - x231 * x250 * xi_y - x232 * x61 * x9 - x232 * (2 * x4 + x5 * x59 + x59 + 2 * x7) + x237 * x242 - x237 * x244 + x239 * x240 + x248 * x249 + x252 * xi_x
	end
AA[2,1] = 0
AA[2,2] = x12
BB[2,1] = x256
BB[2,2] = x256
CC[2,1] = x265
CC[2,2] = x265
if u==0.0
	println("entered fy u zero")
	SS[2] =2*S0^2*(Bp_y) + 2*S0^2*(Gp_x) - 2*Gpp*S0*(S0_tx) - 2*Bpp*S0*(S0_ty) + 4*Spp*(S0_ty) + 4*Gp*S0*(S0_x) + 2*Gpp*(S0_t)*(S0_x) + 4*Bp*S0*(S0_y) + 4*Sp*(S0_y) + 2*Bpp*(S0_t)*(S0_y) - (4*Spp*(S0_t)*(S0_y))/S0 - 4*S0*(Sp_y) - 2*Gpp*S0^2*(xi_x) - 2*Bpp*S0^2*(xi_y) + 4*S0*Spp*(xi_y)
else
	SS[2] = -B_x * x276 - B_y * x184 * x205 + B_y * x239 + Bp_y * x206 - G_x * x212 + G_x * x23 * x246 + Gp_x * x245 + S0_tx * x283 + S0_ty * x269 - 4 * S0_ty * x275 + S0_y * x269 * xi + S_x * x236 * x281 - S_y * x235 - Sp_y * x272 + x188 * x237 * x275 + x188 * x266 - 4 * x19 * x274 * x5 + x196 * x254 * x59 * xi_x - x197 * x207 + x197 * x238 * x240 + x198 * x267 + x199 * x242 + x199 * x244 - 4 * x200 * x237 - x200 * x271 + x202 * x240 + x216 * x279 - x218 * x277 - x218 * x280 + x218 * x284 + x221 * x279 + x237 * x248 + x237 * x273 - x237 * x278 + x242 * x249 + x243 * x279 + x244 * x249 + x248 * xi_y - x250 * x271 * x275 + x251 * x283 + x253 * x283 + 8 * x266 * x9 - x267 * x275 + x268 * x270 - x270 * x272 + x273 * xi_y + x274 * (-B * x210 + 4 * Bp) - x277 * xi_x - x278 * xi_y - x280 * xi_x + x282 * x61 + 2 * x282 + x284 * xi_x
 end


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
	println("entered Sd u zero")
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

 x0 = u ^ 2
x1 = S0_t * u
x2 = S0 * xi
x3 = u * x2
x4 = u ^ 3
x5 = S * x4
x6 = S0 + x3 + x5
x7 = x1 + x6
x8 = x7 ^ 3
x9 = B * x4
x10 = exp(x9)
x11 = 8 * x10
x12 = -x0 * x11 * x8
x13 = x7 ^ 2
x14 = G * x4
x15 = 3 * u
x16 = G * x15
x17 = -Gp + x16
x18 = x0 * x17 * tanh(x14)
x19 = -u * x11 * x13 * (S0 - Sp * x0 + x1 * (x18 + 2) + x18 * x6 + 2 * x3 + 4 * x5)
x20 = S0 ^ 3
x21 = 1 / x20
x22 = sech(x14)
x23 = S0 * u
x24 = x22 * x23
x25 = u ^ 4
x26 = S0 * x25
x27 = 2 * S * x26
x28 = S0 ^ 2
x29 = 2 * x5
x30 = x1 * (S0 - x29) + x28 * (u * xi + 1)
x31 = S * x0
x32 = 2 * xi_y
x33 = x0 * x28
x34 = S0_ty * x27 + S0_y * x30 + x33 * (-Fy * S0 + Fy * x29 + Sh * u + x31 * x32)
x35 = -Fyp
x36 = Fy * u
x37 = -x35 - x36
x38 = x34 * x37
x39 = x35 + 2 * x36
x40 = -x39
x41 = S0_y * x1
x42 = S0_ty * x23
x43 = 5 * x0
x44 = x7 * (x28 * (-2 * Fyph - u * x32 * x40 + u * (Fy ^ 2 * x43 + 2 * Fyh + Fyp ^ 2 - 4 * Fyp * x36)) + 2 * x39 * x42 + 2 * x40 * x41)
x45 = 2 * x24
x46 = B * x15 - Bp
x47 = 4 * x10
x48 = S0 ^ 4
x49 = 2 * x20
x50 = x13 * (Sd * x4 * x49 + u * x49 * (S0_t + x2) + x0 * (-S0 * S0_xx - S0 * S0_yy + S0_t ^ 2 * x28 + S0_t * x49 * xi + S0_x ^ 2 + S0_y ^ 2 + x48 * xi ^ 2) + x48)
x51 = Fx * u
x52 = 2 * xi_x
x53 = S0_tx * x27 + S0_x * x30 + x33 * (-Fx * S0 + Fx * x29 + St * u + x31 * x52)
x54 = -Fxp
x55 = 2 * x51 + x54
x56 = -x55
x57 = S0_x * x1
x58 = S0_tx * x23
x59 = (x53 * (4 * Fxp - 4 * x51) + x7 * (x28 * (-2 * Fxpt - u * x52 * x56 + u * (Fx ^ 2 * x43 + Fxp ^ 2 - 4 * Fxp * x51 + 2 * Fxt)) + 2 * x55 * x58 + 2 * x56 * x57)) * exp(2 * x9)
x60 = Gh * x28
x61 = Gp * x28
x62 = Fx * Gp
x63 = Gt * x28
x64 = 3 * G
x65 = Fxp * x64
x66 = Fy * Gp
x67 = Fyp * x64
x68 = x16 * x28
x69 = u * x28
x70 = Fxp * xi_y
x71 = Fyp * xi_x
x72 = Fx * x33
x73 = Fy * x33
x74 = 2 * x10
x75 = x4 * x46 * x74 * x8 * sinh(2 * x14)
x76 = x23 * sinh(x14)
x77 = x23 * cosh(x14)
x78 = Fxp * x69
x79 = 2 * x0
x80 = S0_t * S0_y
x81 = Fx * x80
x82 = S0_t * S0_x
x83 = Fy * x82
x84 = Bp * x25
x85 = x28 * x4
x86 = Bp * x85
x87 = x25 * x28
x88 = Fy * x87
x89 = 3 * B
x90 = Fxp * x89
x91 = Fyp * x89
x92 = x87 * x89
x93 = Fx * x87
x94 = u ^ 5 * x28
x95 = 3 * x9
x96 = Bp * x0
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = x12
BB[1,2] = x12
CC[1,1] = x19
CC[1,2] = x19
SS[1] = x21 * (x22 * x26 * x47 * x7 * (-Fxh * x61 + Fxh * x68 + Fxp * x60 - Fyp * x63 + Fyt * x61 - Fyt * x68 + x36 * x61 * xi_x + x36 * x63 + x41 * x62 - x41 * x65 + x42 * (3 * Fxp * G - x62) - x51 * x60 - x57 * x66 + x57 * x67 + x58 * (x66 - x67) - x62 * x69 * xi_y + x65 * x73 - x67 * x72 + x68 * x70 - x68 * x71) + 8 * x24 * x38 + x44 * x45 - x45 * x59 - x46 * x47 * x50)
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = 0
CC[2,1] = x75
CC[2,2] = x75
SS[2] = -x21 * (x17 * x50 * x74 + 4 * x38 * x76 + x44 * x76 - x47 * x77 * (x34 * (-x51 - x54) + x37 * x53) + x59 * x76 + x7 * x74 * x77 * (Bh * Fxp * x85 - Bh * x93 + Bp * x88 * xi_x - Bp * x93 * xi_y - Bt * Fyp * x85 + Bt * x88 - 5 * Fx * Fy * x85 - Fx * x91 * x94 - Fxh * x69 - Fxh * x86 + Fxh * x92 - Fxp * x41 + 2 * Fxp * x73 + Fxph * x28 + Fy * x90 * x94 - Fyp * x57 + 2 * Fyp * x72 - Fyp * x78 + Fypt * x28 - Fyt * x69 + Fyt * x86 - Fyt * x92 - x25 * x80 * x90 + x25 * x82 * x91 - x32 * x72 + x42 * (Fxp * x95 + Fxp - x51 * (x96 + 2)) - x52 * x73 + x58 * (-Fyp * x95 + Fyp + x36 * (x96 - 2)) + x69 * x71 + x70 * x92 - x71 * x92 + x78 * xi_y + x79 * x81 + x79 * x83 + x81 * x84 - x83 * x84))
    
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
