using Plots
using LaTeXStrings
using SharedArrays
using Distributed
#plotlyjs()
gr()
#plotly()
#addprocs(72)
function ACBCA(da, db, lc, na, nb, nc, ww, Δ, yeta)
    # 矩阵中有w[nn]会莫名报错数组越界
    c = 3e8;
    D = da + db;
    ka = ww * 1e15 * na / c;
    kb = ww * 1e15 * nb / c;
    kc = ww * 1e15 * nc / c;
    ω₀ = c * pi / (2 * nc * lc);
    Δω = ww .- ω₀;

    num1 = length(ww);
    num2 = length(Δ);
    r0 = SharedArray{Float64}(num2, num1);
    Arg = SharedArray{Float64}(num2, num1);
    rr = CartesianIndices(size(r0));
    @sync begin#同步
        @distributed for k in rr
            mm, nn=k.I;

            m = [exp(im*ka[nn]*(1-Δ[mm])*da/2) exp(-im*ka[nn]*(1-Δ[mm])*da/2); ka[nn]*exp(im*ka[nn]*(1-Δ[mm])*da/2) -ka[nn]*exp(-im*ka[nn]*(1-Δ[mm])*da/2)];
            n = [1 1; kc[nn] -kc[nn]];
            o = [1+im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15) im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15); -im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15) 1-im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15)];
            #o = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
            p = [1 1; kc[nn] -kc[nn]];
            q = [1 1; ka[nn] -ka[nn]];
            r = [exp(im*ka[nn]*(1+Δ[mm])*da/2) exp(-im*ka[nn]*(1+Δ[mm])*da/2); ka[nn]*exp(im*ka[nn]*(1+Δ[mm])*da/2) -ka[nn]*exp(-im*ka[nn]*(1+Δ[mm])*da/2)];
            s = [1 1; kb[nn] -kb[nn]];
            t = [exp(im*kb[nn]*(1+Δ[mm])*db/2) exp(-im*kb[nn]*(1+Δ[mm])*db/2); kb[nn]*exp(im*kb[nn]*(1+Δ[mm])*db/2) -kb[nn]*exp(-im*kb[nn]*(1+Δ[mm])*db/2)];
            u = [1 1; kc[nn] -kc[nn]];
            v = [1+im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15) im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15); -im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15) 1-im*tan(kc[nn]*lc)*Δω[nn]/2/(Δω[nn]+im*yeta*ww[nn]*1e15)];
            #v = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
            w = [1 1; kc[nn] -kc[nn]];
            x = [1 1; kb[nn] -kb[nn]];
            y = [exp(im*kb[nn]*(1-Δ[mm])*db/2) exp(-im*kb[nn]*(1-Δ[mm])*db/2); kb[nn]*exp(im*kb[nn]*(1-Δ[mm])*db/2) -kb[nn]*exp(-im*kb[nn]*(1-Δ[mm])*db/2)];
            z = [1 1; ka[nn] -ka[nn]];
            T =  inv(z) * y * inv(x) * w * v * inv(u) * t * inv(s) * r * inv(q) * p * o * inv(n) * m;
            AAA = T[1, 1] + T[2, 2];
            Temp0 = AAA^2 - 4;
            Temp1 = (AAA - sqrt(Temp0)) / 2;
            r0[mm, nn] = abs((Temp1 - T[1, 1]) / T[1, 2]);
            if r0[mm, nn] <= 1#导带取一种情况
                Temp1 = (AAA - sqrt(Temp0)) / 2;
                r0[mm, nn] = abs((Temp1 - T[1, 1]) / T[1, 2]);
                Arg[mm, nn] = angle((Temp1-T[1, 1]) / T[1, 2]);
            else#此时并不都是禁带，r0仍可能为1，所以取出r0等于1的值与上面导带的情况简并
                Temp1 = (AAA + sqrt(Temp0)) / 2;
                r0[mm, nn] = abs((Temp1 - T[1, 1]) / T[1, 2]);
                if r0[mm, nn] <= 1.0001 && r0[mm, nn] >= 0.9999#1为int所以加一些小数位
                    Temp1 = (AAA - sqrt(Temp0)) / 2;
                    r0[mm, nn] = abs((Temp1 - T[1, 1]) / T[1, 2]);
                    Arg[mm, nn] = angle((Temp1 - T[1, 1]) / T[1, 2]);
                else
                    Temp1 = (AAA + sqrt(Temp0)) / 2;
                    r0[mm, nn] = abs((Temp1 - T[1, 1]) / T[1, 2]);
                    Arg[mm, nn] = angle((Temp1 - T[1, 1])/T[1, 2]);
                end
            end
            if Arg[mm, nn]<-pi+0.001&&Arg[mm, nn]>-pi-0.001
                Arg[mm, nn] = pi;
            end
            #r0[mm,nn]=abs(-T[2,1]/T[2,2]);
            #Arg[mm,nn]=angle(-T[2,1]/T[2,2]);
        end
    end
    return r0, Arg;
end
Arg = nothing; dx = nothing; r0 = nothing; w = nothing; dataPCLs = nothing; f1 = nothing; f2 = nothing;
da = 100e-9; db = 100e-9; lc = da / 2;  # 单位m
na = 2;  # slab A
nb = 2;  # slab B
nc = 2;
w = 0.00001:0.001:5;
Δ = -1:0.001:1;
yeta = 0;

r0, Arg = ACBCA(da, db, lc, na, nb, nc, w, Δ, yeta)
f1 = heatmap(Δ, w, Arg', color = :bluesreds,
title = L"Arg(r)",
ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
xlabel = L"\Delta (nm)")

f2 = heatmap(Δ, w, r0', color = :bluesreds,
title = L"|r|",
xlabel = L"\Delta (nm)")

plot(f1, f2, layout = (1, 2))
