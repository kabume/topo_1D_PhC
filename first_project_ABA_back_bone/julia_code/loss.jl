using Plots
using LaTeXStrings
using SharedArrays
using Distributed
#plotlyjs()
gr()
#addprocs(8)  # 在Atom中，除了第一次运行，下次运行的时候要注释掉
function ABA(da, db, lc, na, nb, nc, w, dx, yeta)
    c = 3e8;
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;
    kc = w * 1e15 * nc / c;
    ω₀ = c * pi / (2 * nc * lc);
    Δω = w .- ω₀;

    num1 = length(w);
    num2 = length(dx);
    r0 = SharedArray{Float64}(num2, num1);
    t0 = SharedArray{Float64}(num2, num1);
    rr = CartesianIndices(size(r0));
    @sync begin#同步
        @distributed for nn = 1:num1
            m = [exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2); ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
            n = [1 1; kb[nn] -kb[nn]];
            o = [exp(im*kb[nn]*(dx-da/2)) exp(-im*kb[nn]*(dx-da/2)); kb[nn]*exp(im*kb[nn]*(dx-da/2)) -kb[nn]*exp(-im*kb[nn]*(dx-da/2))];
            p = [1 1;kc[nn] -kc[nn]];
            q = [1+im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*w[nn]*1e15) im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*w[nn]*1e15); -im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*w[nn]*1e15) 1-im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*w[nn]*1e15)];
            #q = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
            r = [1 1; kc[nn] -kc[nn]];
            s = [1 1; kb[nn] -kb[nn]];
            t = [exp(im*kb[nn]*(D-dx-da/2)) exp(-im*kb[nn]*(D-dx-da/2)); kb[nn]*exp(im*kb[nn]*(D-dx-da/2)) -kb[nn]*exp(-im*kb[nn]*(D-dx-da/2))];
            u = [1 1; ka[nn] -ka[nn]];
            v = [exp(im*ka[nn]*da/2) 0; 0 exp(-im*ka[nn]*da/2)];
            T11 = v * inv(u) * t * inv(s) * r * q * inv(p) * o * inv(n) * m;
            T = T11^10;
            r0[nn] = abs(-T[2,1]./T[2,2]);
            t0[nn] = abs(1 ./T[2,2]);
        end
    end
    return r0, t0;
    end
Arg = nothing; dx = nothing; r0 = nothing; w = nothing; dataPCLs = nothing; f1 = nothing; f2 = nothing;
da = 100e-9; db = 100e-9; lc = da / 2;  # 单位m
lc1 = da / 2; lc2 = da / 1.2;
gamma = 0;
na = 2 + gamma * im;  # slab A
nb = 5 + gamma * im;  # slab B
nc = 2;
w = 0.0001:0.001:4;
#w = 1.4:0.001:1.6
#dx = 0:0.1e-9:da+db;
dx = 100e-9;
yeta = 0;
r0, t0 = ABA(da, db, lc, na, nb, nc, w, dx, yeta)

#gif(anim, "/tmp/anim_fps15.gif", fps = 15)
plot(w, (r0)',xlabel = L"f\, (10^{15}/ 2 \pi\ Hz)", ylabel = L"|r|,|t|", dpi = 600,label = L"|r|")
plot!(w, (t0)',label = L"|t|")
savefig("E:\\XLL\\OneDrive - email.ncu.edu.cn\\FDU\\topological\\manuscript\\image\\da_100_db_100_lc_50_na_2_nb_5_nc_2.png")
