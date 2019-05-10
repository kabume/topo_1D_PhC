# TODO:没有画右边十层的图。
using Plots
using LaTeXStrings
gr()
#plotly()
function E_file(da, db, lc, na, nb, nc, w, dx)
    yeta = 0;
    c = 3e8;
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;
    kc = w * 1e15 * nc / c;
    ω₀ = c * pi / (2 * nc * lc);
    Δω = w .- ω₀;

    z11 = 0:0.1e-9:da/2; z21 = da/2:0.1e-9:dx; z31 = dx:0.1e-9:da/2+db; z41 = da/2+db:0.1e-9:da+db;
    z12 = z11 .+ D; z22 = z21 .+ D; z32 = z31 .+ D; z42 = z41 .+ D;
    z13 = z11 .+ 2D; z23 = z21 .+ 2D; z33 = z31 .+ 2D; z43 = z41 .+ 2D;
    z14 = z11 .+ 3D; z24 = z21 .+ 3D; z34 = z31 .+ 3D; z44 = z41 .+ 3D;
    z15 = z11 .+ 4D; z25 = z21 .+ 4D; z35 = z31 .+ 4D; z45 = z41 .+ 4D;
    z16 = z11 .+ 5D; z26 = z21 .+ 5D; z36 = z31 .+ 5D; z46 = z41 .+ 5D;
    z17 = z11 .+ 6D; z27 = z21 .+ 6D; z37 = z31 .+ 6D; z47 = z41 .+ 6D;
    z18 = z11 .+ 7D; z28 = z21 .+ 7D; z38 = z31 .+ 7D; z48 = z41 .+ 7D;
    z19 = z11 .+ 8D; z29 = z21 .+ 8D; z39 = z31 .+ 8D; z49 = z41 .+ 8D;
    z10 = z11 .+ 9D; z20 = z21 .+ 9D; z30 = z31 .+ 9D; z40 = z41 .+ 9D;

    x11 = z11 .+ 10D; x21 = z21 .+ 10D; x31 = z31 .+ 10D; x41 = z41 .+ 10D;
    x12 = x11 .+ D; x22 = x21 .+ D; x32 = x31 .+ D; x42 = x41 .+ D;
    x13 = x11 .+ 2D; x23 = x21 .+ 2D; x33 = x31 .+ 2D; x43 = x41 .+ 2D;
    x14 = x11 .+ 3D; x24 = x21 .+ 3D; x34 = x31 .+ 3D; x44 = x41 .+ 3D;
    x15 = x11 .+ 4D; x25 = x21 .+ 4D; x35 = x31 .+ 4D; x45 = x41 .+ 4D;
    x16 = x11 .+ 5D; x26 = x21 .+ 5D; x36 = x31 .+ 5D; x46 = x41 .+ 5D;
    x17 = x11 .+ 6D; x27 = x21 .+ 6D; x37 = x31 .+ 6D; x47 = x41 .+ 6D;
    x18 = x11 .+ 7D; x28 = x21 .+ 7D; x38 = x31 .+ 7D; x48 = x41 .+ 7D;
    x19 = x11 .+ 8D; x29 = x21 .+ 8D; x39 = x31 .+ 8D; x49 = x41 .+ 8D;
    x10 = x11 .+ 9D; x20 = x21 .+ 9D; x30 = x31 .+ 9D; x40 = x41 .+ 9D;
    x1 = [x11, x21, x31, x41, x12, x22, x32, x42];
    x2 = [x13, x23, x33, x43, x14, x24, x34, x44];
    x3 = [x15, x25, x35, x45, x16, x26, x36, x46];
    x4 = [x17, x27, x37, x47, x18, x28, x38, x48];
    x5 = [x19, x29, x39, x49, x10, x20, x30, x40];

    z1 = [z11, z21, z31, z41, z12, z22, z32, z42];
    z2 = [z13, z23, z33, z43, z14, z24, z34, z44];
    z3 = [z15, z25, z35, z45, z16, z26, z36, z46];
    z4 = [z17, z27, z37, z47, z18, z28, z38, z48];
    z5 = [z19, z29, z39, z49, z10, z20, z30, z40];
    #z = [z1, z2, z3, z4, z5, x1, x2, x3, x4, x5];
    z = [z1, z2, z3, z4, z5];
    m = [exp(im*ka*da/2) exp(-im*ka*da/2); ka*exp(im*ka*da/2) -ka*exp(-im*ka*da/2)];
    n = [1 1; kb -kb];
    o = [exp(im*kb*(dx-da/2)) exp(-im*kb*(dx-da/2)); kb*exp(im*kb*(dx-da/2)) -kb*exp(-im*kb*(dx-da/2))];
    p = [1 1;kc -kc];
    #q = [1+im*tan(kc*lc)/2/(w + im * yeta) im*tan(kc*lc)/2/(w + im * yeta); -im*tan(kc*lc)/2/(w + im * yeta) 1-im*tan(kc*lc)/2/(w + im * yeta)];
    q = [1+im*tan(kc*lc)/2 im*tan(kc*lc)/2; -im*tan(kc*lc)/2 1-im*tan(kc*lc)/2];
    r = [1 1; kc -kc];
    s = [1 1; kb -kb];
    t = [exp(im*kb*(D-dx-da/2)) exp(-im*kb*(D-dx-da/2)); kb*exp(im*kb*(D-dx-da/2)) -kb*exp(-im*kb*(D-dx-da/2))];
    u = [1 1; ka -ka];
    v = [exp(im*ka*da/2) 0; 0 exp(-im*ka*da/2)];
    T1 = v * inv(u) * t * inv(s) * r * q * inv(p) * o * inv(n) * m;

    m1 = [1+im*tan(kc*lc)*Δω/2/(Δω+im*yeta*w*1e15) im*tan(kc*lc)*Δω/2/(Δω+im*yeta*w*1e15); -im*tan(kc*lc)*Δω/2/(Δω+im*yeta*w*1e15) 1-im*tan(kc*lc)*Δω/2/(Δω+im*yeta*w*1e15)];
    n1 = [1 1; kc -kc];
    o1 = [1 1; kb -kb];
    p1 = [exp(im*kb*(D-dx-db/2)) exp(-im*kb*(D-dx-db/2)); kb*exp(im*kb*(D-dx-db/2)) -kb*exp(-im*kb*(D-dx-db/2))];
    q1 = [1 1; ka -ka];
    r1 = [exp(im*ka*da) exp(-im*ka*da); ka*exp(im*ka*da) -ka*exp(-im*ka*da)];
    s1 = [1 1; kb -kb];
    t1 = [exp(im*kb*(dx-db/2)) exp(-im*kb*(dx-db/2)); kb*exp(im*kb*(dx-db/2)) -kb*exp(-im*kb*(dx-db/2))];
    u1 = [1 1;kc -kc];
    T2 = inv(u1) * t1 * inv(s1) * r1 * inv(q1) * p1 * inv(o1) * n1 * m1;


    T = T1^10 * T2^10;
    s11 =  [T[2, 2]; -T[2, 1]]
    E11 =  s11[1] * exp.(im * ka * z11) +  s11[2] * exp.(-im * ka *z11);
    s21 = inv(n) * m * s11;
    E21 = s21[1] * exp.(im * kb * (z21 .- da / 2)) + s21[2] * exp.(-im * kb * (z21 .- da / 2));
    s31 = inv(s) * r * q * inv(p) * o * s21;
    E31 = s31[1] * exp.(im * kb * (z31 .- dx)) + s31[2] * exp.(-im * kb * (z31 .- dx));
    s41 = inv(u) * t * s31;
    E41 = s41[1] * exp.(im * ka * (z41 .- da / 2 .- db)) + s41[2] * exp.(-im * ka * (z41 .- da / 2 .- db));

    s12 = v * s41;
    E12 =  s12[1] * exp.(im * ka * (z12.-D)) +  s12[2] * exp.(-im * ka *(z12.-D));
    s22 = inv(n) * m * s12;
    E22 = s22[1] * exp.(im * kb * (z22 .- da / 2 .- D)) + s22[2] * exp.(-im * kb * (z22 .- da / 2 .- D));
    s32 = inv(s) * r * q * inv(p) * o * s22;
    E32 = s32[1] * exp.(im * kb * (z32 .- dx .- D)) + s32[2] * exp.(-im * kb * (z32 .- dx .- D));
    s42 = inv(u) * t * s32;
    E42 = s42[1] * exp.(im * ka * (z42 .- da / 2 .- db .- D)) + s42[2] * exp.(-im * ka * (z42 .- da / 2 .- db .- D));

    s13 = v * s42;
    E13 = s13[1] * exp.(im * ka * (z13.-2D)) +  s13[2] * exp.(-im * ka *(z13.-2D));
    s23 = inv(n) * m * s13;
    E23 = s23[1] * exp.(im * kb * (z23 .- da / 2 .- 2D)) + s23[2] * exp.(-im * kb * (z23 .- da / 2 .- 2D));
    s33 = inv(s) * r * q * inv(p) * o * s23;
    E33 = s33[1] * exp.(im * kb * (z33 .- dx .- 2D)) + s33[2] * exp.(-im * kb * (z33 .- dx .- 2D));
    s43 = inv(u) * t * s33;
    E43 = s43[1] * exp.(im * ka * (z43 .- da / 2 .- db .- 2D)) + s43[2] * exp.(-im * ka * (z43 .- da / 2 .- db .- 2D));

    s14 = v * s43;
    E14 =  s14[1] * exp.(im * ka * (z14.-3D)) +  s14[2] * exp.(-im * ka *(z14.-3D));
    s24 = inv(n) * m * s14;
    E24 = s24[1] * exp.(im * kb * (z24 .- da / 2 .- 3D)) + s24[2] * exp.(-im * kb * (z24 .- da / 2 .- 3D));
    s34 = inv(s) * r * q * inv(p) * o * s24;
    E34 = s34[1] * exp.(im * kb * (z34 .- dx .- 3D)) + s34[2] * exp.(-im * kb * (z34 .- dx .- 3D));
    s44 = inv(u) * t * s34;
    E44 = s44[1] * exp.(im * ka * (z44 .- da / 2 .- db .- 3D)) + s44[2] * exp.(-im * ka * (z44 .- da / 2 .- db .- 3D));

    s15 = v * s44;
    E15 =  s15[1] * exp.(im * ka * (z15.-4D)) +  s15[2] * exp.(-im * ka *(z15.-4D));
    s25 = inv(n) * m * s15;
    E25 = s25[1] * exp.(im * kb * (z25 .- da / 2 .- 4D)) + s25[2] * exp.(-im * kb * (z25 .- da / 2 .- 4D));
    s35 = inv(s) * r * q * inv(p) * o * s25;
    E35 = s35[1] * exp.(im * kb * (z35 .- dx .- 4D)) + s35[2] * exp.(-im * kb * (z35 .- dx .- 4D));
    s45 = inv(u) * t * s35;
    E45 = s45[1] * exp.(im * ka * (z45 .- da / 2 .- db .- 4D)) + s45[2] * exp.(-im * ka * (z45 .- da / 2 .- db .- 4D));

    s16 = v * s45;
    E16 =  s16[1] * exp.(im * ka * (z16.-5D)) +  s16[2] * exp.(-im * ka *(z16.-5D));
    s26 = inv(n) * m * s16;
    E26 = s26[1] * exp.(im * kb * (z26 .- da / 2 .- 5D)) + s26[2] * exp.(-im * kb * (z26 .- da / 2 .- 5D));
    s36 = inv(s) * r * q * inv(p) * o * s26;
    E36 = s36[1] * exp.(im * kb * (z36 .- dx .- 5D)) + s36[2] * exp.(-im * kb * (z36 .- dx .- 5D));
    s46 = inv(u) * t * s36;
    E46 = s46[1] * exp.(im * ka * (z46 .- da / 2 .- db .- 5D)) + s46[2] * exp.(-im * ka * (z46 .- da / 2 .- db .- 5D));

    s17 = v * s46;
    E17 =  s17[1] * exp.(im * ka * (z17.-6D)) +  s17[2] * exp.(-im * ka *(z17.-6D));
    s27 = inv(n) * m * s17;
    E27 = s27[1] * exp.(im * kb * (z27 .- da / 2 .- 6D)) + s27[2] * exp.(-im * kb * (z27 .- da / 2 .- 6D));
    s37 = inv(s) * r * q * inv(p) * o * s27;
    E37 = s37[1] * exp.(im * kb * (z37 .- dx .- 6D)) + s37[2] * exp.(-im * kb * (z37 .- dx .- 6D));
    s47 = inv(u) * t * s37;
    E47 = s47[1] * exp.(im * ka * (z47 .- da / 2 .- db .- 6D)) + s47[2] * exp.(-im * ka * (z47 .- da / 2 .- db .- 6D));

    s18 = v * s47;
    E18 =  s18[1] * exp.(im * ka * (z18.-7D)) +  s18[2] * exp.(-im * ka *(z18.-7D));
    s28 = inv(n) * m * s18;
    E28 = s28[1] * exp.(im * kb * (z28 .- da / 2 .- 7D)) + s28[2] * exp.(-im * kb * (z28 .- da / 2 .- 7D));
    s38 = inv(s) * r * q * inv(p) * o * s28;
    E38 = s38[1] * exp.(im * kb * (z38 .- dx .- 7D)) + s38[2] * exp.(-im * kb * (z38 .- dx .- 7D));
    s48 = inv(u) * t * s38;
    E48 = s48[1] * exp.(im * ka * (z48 .- da / 2 .- db .- 7D)) + s48[2] * exp.(-im * ka * (z48 .- da / 2 .- db .- 7D));

    s19 = v * s48;
    E19 =  s19[1] * exp.(im * ka * (z19.-8D)) +  s19[2] * exp.(-im * ka *(z19.-8D));
    s29 = inv(n) * m * s19;
    E29 = s29[1] * exp.(im * kb * (z29 .- da / 2 .- 8D)) + s29[2] * exp.(-im * kb * (z29 .- da / 2 .- 8D));
    s39 = inv(s) * r * q * inv(p) * o * s29;
    E39 = s39[1] * exp.(im * kb * (z39 .- dx .- 8D)) + s39[2] * exp.(-im * kb * (z39 .- dx .- 8D));
    s49 = inv(u) * t * s39;
    E49 = s49[1] * exp.(im * ka * (z49 .- da / 2 .- db .- 8D)) + s49[2] * exp.(-im * ka * (z49 .- da / 2 .- db .- 8D));

    s10 = v * s49;
    E10 =  s10[1] * exp.(im * ka * (z10.-9D)) +  s10[2] * exp.(-im * ka *(z10.-9D));
    s20 = inv(n) * m * s10;
    E20 = s20[1] * exp.(im * kb * (z20 .- da / 2 .- 9D)) + s20[2] * exp.(-im * kb * (z20 .- da / 2 .- 9D));
    s30 = inv(s) * r * q * inv(p) * o * s20;
    E30 = s30[1] * exp.(im * kb * (z30 .- dx .- 9D)) + s30[2] * exp.(-im * kb * (z30 .- dx .- 9D));
    s40 = inv(u) * t * s30;
    E40 = s40[1] * exp.(im * ka * (z40 .- da / 2 .- db .- 9D)) + s40[2] * exp.(-im * ka * (z40 .- da / 2 .- db .- 9D));


    dataE1 = [abs.(E11), abs.(E21), abs.(E31), abs.(E41),abs.(E12), abs.(E22), abs.(E32), abs.(E42)];
    dataE2 = [abs.(E13), abs.(E23), abs.(E33), abs.(E43),abs.(E14), abs.(E24), abs.(E34), abs.(E44)];
    dataE3 = [abs.(E15), abs.(E25), abs.(E35), abs.(E45),abs.(E16), abs.(E26), abs.(E36), abs.(E46)];
    dataE4 = [abs.(E17), abs.(E27), abs.(E37), abs.(E47),abs.(E18), abs.(E28), abs.(E38), abs.(E48)];
    dataE5 = [abs.(E19), abs.(E29), abs.(E39), abs.(E49),abs.(E10), abs.(E20), abs.(E30), abs.(E40)];
    dataE = [dataE1, dataE2, dataE3, dataE4, dataE5];

    return z,dataE;
end

da = 100e-9; db = 100e-9; lc = da;  # 单位m
na = 3.2;  # slab A
nb = 1;  # slab B
nc = 2;

dx = 100e-9;
#w = 1.557087;
w = 4;
z, dataE = E_file(da, db, lc, na, nb, nc, w, dx)
plot(z .* 1e9, dataE, legend = :none, dpi = 600,
ylabel = L"E(x)",
xlabel = L"x (nm)")
#savefig("E:\\OneDrive - email.ncu.edu.cn\\FDU\\topological\\manuscript\\image\\interface_efile.png")
