#并行版本 ABA结构
using Plots
using LaTeXStrings
using SharedArrays
using Distributed
plotlyjs()
#gr
addprocs(72)  # 在Atom中，除了第一次运行，下次运行的时候要注释掉
function ABA(da, db, lc, na, nb, nc, w, dx)
    c = 3e8;
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;
    kc = w * 1e15 * nc / c;
    num1 = length(w);
    num2 = length(dx);
    r0 = SharedArray{Float64}(num2, num1);
    Arg = SharedArray{Float64}(num2, num1);
    rr = CartesianIndices(size(r0));
    @sync begin#同步
        @distributed for k in rr
            mm, nn=k.I;
            if mm < da / 2D * num2
                m = [exp(im*ka[nn]*dx[mm]) exp(-im*ka[nn]*dx[mm]); ka[nn]*exp(im*ka[nn]*dx[mm]) -ka[nn]*exp(-im*ka[nn]*dx[mm])];
                n = [1 1; kc[nn] -kc[nn]];
                o = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
                p = [1 1; kc[nn] -kc[nn]];
                q = [1 1; ka[nn] -ka[nn]];
                r = [exp(im*ka[nn]*(da/2-dx[mm])) exp(-im*ka[nn]*(da/2-dx[mm])); ka[nn]*exp(im*ka[nn]*(da/2-dx[mm])) -ka[nn]*exp(-im*ka[nn]*(da/2-dx[mm]))];
                s = [1 1; kb[nn] -kb[nn]];
                t = [exp(im*kb[nn]*db) exp(-im*kb[nn]*db); kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
                u = [1 1; ka[nn] -ka[nn]];
                v = [exp(im*ka[nn]*da/2) 0; 0 exp(-im*ka[nn]*da/2)];
                T = v * inv(u) * t * inv(s) * r * inv(q) * p * o * inv(n) * m;
            elseif mm >= da/2D*num2 && mm < (da/2+db)/D*num2
                m = [exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2); ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
                n = [1 1; kb[nn] -kb[nn]];
                o = [exp(im*kb[nn]*(dx[mm]-da/2)) exp(-im*kb[nn]*(dx[mm]-da/2)); kb[nn]*exp(im*kb[nn]*(dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(dx[mm]-da/2))];
                p = [1 1;kc[nn] -kc[nn]];
                q = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
                r = [1 1; kc[nn] -kc[nn]];
                s = [1 1; kb[nn] -kb[nn]];
                t = [exp(im*kb[nn]*(D-dx[mm]-da/2)) exp(-im*kb[nn]*(D-dx[mm]-da/2)); kb[nn]*exp(im*kb[nn]*(D-dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(D-dx[mm]-da/2))];
                u = [1 1; ka[nn] -ka[nn]];
                v = [exp(im*ka[nn]*da/2) 0; 0 exp(-im*ka[nn]*da/2)];
                T = v * inv(u) * t * inv(s) * r * q * inv(p) * o * inv(n) * m;
            else
                m = [exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2); ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
                n = [1 1; kb[nn] -kb[nn]];
                o = [exp(im*kb[nn]*db) exp(-im*kb[nn]*db); kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
                p = [1 1; ka[nn] -ka[nn]];
                q = [exp(im*ka[nn]*(dx[mm]-da/2-db)) exp(-im*ka[nn]*(dx[mm]-da/2-db)); ka[nn]*exp(im*ka[nn]*(dx[mm]-da/2-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-da/2-db))];
                r = [1 1; kc[nn] -kc[nn]];
                s = [1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2; -im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2];
                t = [1 1; kc[nn] -kc[nn]];
                u = [1 1; ka[nn] -ka[nn]];
                v = [exp(im*ka[nn]*(D-dx[mm])) 0; 0 exp(-im*ka[nn]*(D-dx[mm]))];
                T = v * inv(u) * t * s * inv(r) * q * inv(p) * o * inv(n) * m;
            end
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
        end
    end
    return r0, Arg;
end

function PCLs(da, db, na, nb, w, dx)
    c = 3e8;
    num1 = length(w);
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;

    bj = [-2*pi, -pi, pi, 2*pi];
    Bn = SharedArray{Float64}(length(w), length(bj));
    rr0 = CartesianIndices(size(Bn));

    c1 = atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
    c2 = atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
    c3 = atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
    c4 = atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));

    aj1 = [0, pi, 2*pi, 3*pi, 4*pi];
    An11 = SharedArray{Float64}(length(w), length(aj1));
    An12 = SharedArray{Float64}(length(w), length(aj1));
    rr1 = CartesianIndices(size(An11));
    aj2 = [0, pi, 2*pi, 3*pi, 4*pi];
    An21 = SharedArray{Float64}(length(w), length(aj2));
    An22 = SharedArray{Float64}(length(w), length(aj2));
    rr2 = CartesianIndices(size(An21));

    @sync begin#同步
        @distributed for k in rr0
            nn, mm=k.I;
            Bn[nn, mm] = (((bj[mm] / kb[nn] + db) / 2) + da / 2);
            if Bn[nn, mm] <= da/2+db && Bn[nn, mm]>=da/2
                Bn[nn, mm] = Bn[nn, mm];
            else
                Bn[nn, mm] = NaN;
            end
        end
        @distributed for k in rr1
            nn, mm=k.I;
            An11[nn, mm] = (ka[nn] * da - c3[nn] - aj1[mm]) / 2 / ka[nn];
            An12[nn, mm] = (ka[nn] * da - c4[nn] - aj1[mm]) / 2 / ka[nn];
            if An11[nn, mm] <= da/2 && An11[nn, mm] >= 0
                An11[nn, mm] = An11[nn, mm];
            else
                An11[nn, mm] = NaN;
            end
            if An12[nn, mm] <= da/2 && An12[nn, mm] >= 0
                An12[nn, mm] = An12[nn, mm];
            else
                An12[nn, mm] = NaN;
            end
        end
        @distributed for k in rr2
            nn, mm=k.I;
            An21[nn, mm] = (ka[nn] * da + 2 * ka[nn] * db - c1[nn] + aj2[mm]) / ka[nn] / 2;
            An22[nn, mm] = (ka[nn] * da + 2 * ka[nn] * db - c2[nn] + aj2[mm]) / ka[nn] / 2;
            if An21[nn, mm] > da/2+db && An21[nn, mm] < da+db
                An21[nn, mm] = An21[nn, mm];
            else
                An21[nn, mm] = NaN;
            end
            if An22[nn, mm] > da/2+db && An22[nn, mm] < da+db
                An22[nn, mm] = An22[nn, mm];
            else
                An22[nn, mm] = NaN;
            end
        end
    end
    k_2 = Bn[1:num1, 1]; k_1 = Bn[1:num1, 2]; k1 = Bn[1:num1, 3]; k2 = Bn[1:num1, 4]; k0 = ones(num1) * (D / 2);

    p110 = An11[1:num1, 1]; p111 = An11[1:num1, 2]; p112 = An11[1:num1, 3]; p113 = An11[1:num1, 4]; p114 = An11[1:num1, 5];
    p120 = An12[1:num1, 1]; p121 = An12[1:num1, 2]; p122 = An12[1:num1, 3]; p123 = An12[1:num1, 4]; p124 = An12[1:num1, 5];

    p210 = An21[1:num1, 1]; p211 = An21[1:num1, 2]; p212 = An21[1:num1, 3]; p213 = An21[1:num1, 4]; p214 = An21[1:num1, 5];
    p220 = An22[1:num1, 1]; p221 = An22[1:num1, 2]; p222 = An22[1:num1, 3]; p223 = An22[1:num1, 4]; p224 = An22[1:num1, 5];

    dataA1 = [p110, p112, p114, p120, p122, p124];
    dataA2 = [p210, p212, p214, p220, p222, p224];
    dataB = [k0, k1, k_1, k2, k_2];
    dataPCLs = [dataA1, dataA2, dataB] * 1e9;
end

# 输入参数
Arg = nothing; dx = nothing; r0 = nothing; w = nothing; dataPCLs = nothing; f1 = nothing; f2 = nothing;
da = 100e-9; db = 100e-9; lc = da;  # 单位m
na = 3.2;  # slab A
nb = 1;  # slab B
nc = 2;
w = 0.001:0.01:15;
dx = 0:0.1*1e-9:da + db;

# 函数运算
@time r0, Arg = ABA(da, db, lc, na, nb, nc, w, dx)  # 计算ABA结构+共振柱的反射系数和反射相位
@time dataPCLs = PCLs(da, db, na, nb, w, dx)　　#　计算ABA结构+共振柱的相位切割线

# 画图

heatmap(dx * 1e9, w, Arg', color = :bluesreds, dpi = 600,
title = L"Arg(r)",
ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
xlabel = L"\Delta (nm)")
f1 = plot!(dataPCLs, w, legend = :none,  line = (:dot, 1.5));

heatmap(dx * 1e9, w, r0', color = :bluesreds, dpi=600,
title = L"|r|",
xlabel = L"\Delta (nm)")
f2 = plot!(dataPCLs, w, legend = :none,  line=(:dot, 1.5))

plot(f1, f2, layout = (1, 2))
savefig("test.html")
