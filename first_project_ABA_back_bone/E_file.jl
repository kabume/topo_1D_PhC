using Plots
using LaTeXStrings
gr()
function E_file(da, db, lc, na, nb, nc, w, dx)
    c = 3e8;
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;
    kc = w * 1e15 * nc / c;
    if dx <= da/2
        z1 = 0:0.1e-9:dx;
        z2 = dx:0.1e-9:da/2;
        z3 = da/2:0.1e-9:da/2+db;
        z4 = da/2+db:0.1e-9:da+db;
        z = [z1, z2, z3, z4];
        m = [exp(im*ka*dx) exp(-im*ka*dx); ka*exp(im*ka*dx) -ka*exp(-im*ka*dx)];
        n = [1 1; kc -kc];
        o = [1+im*tan(kc*lc)/2 im*tan(kc*lc)/2; -im*tan(kc*lc)/2 1-im*tan(kc*lc)/2];
        p = [1 1; kc -kc];
        q = [1 1; ka -ka];
        r = [exp(im*ka*(da/2-dx)) exp(-im*ka*(da/2-dx)); ka*exp(im*ka*(da/2-dx)) -ka*exp(-im*ka*(da/2-dx))];
        s = [1 1; kb -kb];
        t = [exp(im*kb*db) exp(-im*kb*db); kb*exp(im*kb*db) -kb*exp(-im*kb*db)];
        u = [1 1; ka -ka];
        v = [exp(im*ka*da/2) 0; 0 exp(-im*ka*da/2)];
        T = v * inv(u) * t * inv(s) * r * inv(q) * p * o * inv(n) * m;

        AAA = T[1, 1] + T[2, 2];
        Temp0 = AAA^2 - 4;
        Temp1 = (AAA + sqrt(Temp0)) / 2;  #+ or -
        s1 =  [T[1, 2]; (Temp1 - T[1, 1])]
        E1 =  s1[1] * exp.(im * ka * z1) +  s1[2] * exp.(-im * ka *z1);
        s2 = inv(q) * p * o * inv(n) * m * s1;
        E2 = s2[1] * exp.(im * ka * (z2 .- dx)) + s2[2] * exp.(-im * ka * (z2 .- dx));
        s3 = inv(s) * r * s2;
        E3 = s3[1] * exp.(im * kb * (z3 .- da/2)) + s3[2] * exp.(-im * kb * (z3 .- da/2));
        s4 = inv(u) * t * s3;
        E4 = s4[1] * exp.(im * ka * (z4 .- da / 2 .- db)) + s4[2] * exp.(-im * ka * (z4 .- da / 2 .- db));
    elseif dx > da/2 && dx <= da/2 + db
        z1 = 0:0.1e-9:da/2;
        z2 = da/2:0.1e-9:dx;
        z3 = dx:0.1e-9:da/2+db;
        z4 = da/2+db:0.1e-9:da+db;
        z = [z1, z2, z3, z4];
        m = [exp(im*ka*da/2) exp(-im*ka*da/2); ka*exp(im*ka*da/2) -ka*exp(-im*ka*da/2)];
        n = [1 1; kb -kb];
        o = [exp(im*kb*(dx-da/2)) exp(-im*kb*(dx-da/2)); kb*exp(im*kb*(dx-da/2)) -kb*exp(-im*kb*(dx-da/2))];
        p = [1 1;kc -kc];
        q = [1+im*tan(kc*lc)/2 im*tan(kc*lc)/2; -im*tan(kc*lc)/2 1-im*tan(kc*lc)/2];
        r = [1 1; kc -kc];
        s = [1 1; kb -kb];
        t = [exp(im*kb*(D-dx-da/2)) exp(-im*kb*(D-dx-da/2)); kb*exp(im*kb*(D-dx-da/2)) -kb*exp(-im*kb*(D-dx-da/2))];
        u = [1 1; ka -ka];
        v = [exp(im*ka*da/2) 0; 0 exp(-im*ka*da/2)];
        T = v * inv(u) * t * inv(s) * r * q * inv(p) * o * inv(n) * m;

        AAA = T[1, 1] + T[2, 2];
        Temp0 = AAA^2 - 4;
        Temp1 = (AAA + sqrt(Temp0)) / 2;  #+ or -
        s1 =  [T[1, 2]; (Temp1 - T[1, 1])]
        E1 =  s1[1] * exp.(im * ka * z1) +  s1[2] * exp.(-im * ka *z1);
        s2 = inv(n) * m * s1;
        E2 = s2[1] * exp.(im * kb * (z2 .- da / 2)) + s2[2] * exp.(-im * kb * (z2 .- da / 2));
        s3 = inv(s) * r * q * inv(p) * o * s2;
        E3 = s3[1] * exp.(im * kb * (z3 .- dx)) + s3[2] * exp.(-im * kb * (z3 .- dx));
        s4 = inv(u) * t * s3;
        E4 = s4[1] * exp.(im * ka * (z4 .- da / 2 .- db)) + s4[2] * exp.(-im * ka * (z4 .- da / 2 .- db));
    else
        z1 = 0:0.1e-9:da/2;
        z2 = da/2:0.1e-9:da/2+db;
        z3 = da/2+db:0.1e-9:dx;
        z4 = dx:0.1e-9:da+db;
        z = [z1, z2, z3, z4];
        m = [exp(im*ka*da/2) exp(-im*ka*da/2); ka*exp(im*ka*da/2) -ka*exp(-im*ka*da/2)];
        n = [1 1; kb -kb];
        o = [exp(im*kb*db) exp(-im*kb*db); kb*exp(im*kb*db) -kb*exp(-im*kb*db)];
        p = [1 1; ka -ka];
        q = [exp(im*ka*(dx-da/2-db)) exp(-im*ka*(dx-da/2-db)); ka*exp(im*ka*(dx-da/2-db)) -ka*exp(-im*ka*(dx-da/2-db))];
        r = [1 1; kc -kc];
        s = [1+im*tan(kc*lc)/2 im*tan(kc*lc)/2; -im*tan(kc*lc)/2 1-im*tan(kc*lc)/2];
        t = [1 1; kc -kc];
        u = [1 1; ka -ka];
        v = [exp(im*ka*(D-dx)) 0; 0 exp(-im*ka*(D-dx))];
        T = v * inv(u) * t * s * inv(r) * q * inv(p) * o * inv(n) * m;


        AAA = T[1, 1] + T[2, 2];
        Temp0 = AAA^2 - 4;
        Temp1 = (AAA + sqrt(Temp0)) / 2;  #+ or -
        s1 =  [T[1, 2]; (Temp1 - T[1, 1])]
        E1 =  s1[1] * exp.(im * ka * z1) +  s1[2] * exp.(-im * ka *z1);
        s2 = inv(n) * m * s1;
        E2 = s2[1] * exp.(im * kb * (z2 .- da / 2)) + s2[2] * exp.(-im * kb * (z2 .- da / 2));
        s3 = inv(p) * o * s2;
        E3 = s3[1] * exp.(im * ka * (z3 .- da / 2 .- db)) + s3[2] * exp.(-im * ka * (z3 .- da / 2 .- db));
        s4 = inv(u) * t * s * inv(r) * q * s3;
        E4 = s4[1] * exp.(im * ka * (z4 .- dx)) + s4[2] * exp.(-im * ka * (z4 .- dx));
    end

    dataE = [abs.(E1), abs.(E2), abs.(E3), abs.(E4)];
    return z,dataE;
end

da = 100e-9; db = 100e-9; lc = da;  # 单位m
na = 3.2;  # slab A
nb = 1;  # slab B
nc = 2;
w = 14;
dx = 100e-9;


z, dataE = E_file(da, db, lc, na, nb, nc, w, dx)
plot(z .* 1e9, dataE, legend = :none)
