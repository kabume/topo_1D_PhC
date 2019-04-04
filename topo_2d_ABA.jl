#并行版本
#ABA结构

using Plots
using Unicode
using LaTeXStrings
using SharedArrays
using Distributed
gr()
#pyplot()
Arg=nothing;dx=nothing;dx0=nothing;ka=nothing;kb=nothing;kc=nothing;p11=nothing;p12=nothing;r0=nothing;rr=nothing;w=nothing;R=0;

#addprocs(72)#除了第一次运行，下次运行的时候要注释掉

da=100e-9;db=200e-9;D=da+db;l0=da;#单位m
Dt=D*1e9;
dx0=0:0.1:Dt;
dx=dx0*1e-9
na=2;#slab A
nb=1;#slab B
nc=2;
c=3e8;
w=0.1:0.01:15;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;
num1=length(w);
num2=length(dx);
r0 = zeros(num2,num1);
Arg = zeros(num2,num1);
r0 = SharedArray{Float64}(num2, num1);
Arg = SharedArray{Float64}(num2, num1);
R = SharedArray{Float64}(num2, num1);
rr = CartesianIndices(size(r0));
@sync begin#同步
    @distributed for k in rr
        mm,nn=k.I;
        if mm<da/2D*num2
            m=[exp(im*ka[nn]*dx[mm]) exp(-im*ka[nn]*dx[mm]);ka[nn]*exp(im*ka[nn]*dx[mm]) -ka[nn]*exp(-im*ka[nn]*dx[mm])];
            n=[1 1;kc[nn] -kc[nn]];
            o=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
            p=[1 1;kc[nn] -kc[nn]];
            q=[1 1;ka[nn] -ka[nn]];
            r=[exp(im*ka[nn]*(da/2-dx[mm])) exp(-im*ka[nn]*(da/2-dx[mm]));ka[nn]*exp(im*ka[nn]*(da/2-dx[mm])) -ka[nn]*exp(-im*ka[nn]*(da/2-dx[mm]))];
            s=[1 1;kb[nn] -kb[nn]];
            t=[exp(im*kb[nn]*db) exp(-im*kb[nn]*db);kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
            u=[1 1;ka[nn] -ka[nn]];
            v=[exp(im*ka[nn]*da/2) 0;0 exp(-im*ka[nn]*da/2)];
            T=v*inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
        elseif mm>=da/2D*num2 && mm<(da/2+db)/D*num2
            m=[exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2);ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
            n=[1 1;kb[nn] -kb[nn]];
            o=[exp(im*kb[nn]*(dx[mm]-da/2)) exp(-im*kb[nn]*(dx[mm]-da/2));kb[nn]*exp(im*kb[nn]*(dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(dx[mm]-da/2))];
            p=[1 1;kc[nn] -kc[nn]];
            q=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
            r=[1 1;kc[nn] -kc[nn]];
            s=[1 1;kb[nn] -kb[nn]];
            t=[exp(im*kb[nn]*(D-dx[mm]-da/2)) exp(-im*kb[nn]*(D-dx[mm]-da/2));kb[nn]*exp(im*kb[nn]*(D-dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(D-dx[mm]-da/2))];
            u=[1 1;ka[nn] -ka[nn]];
            v=[exp(im*ka[nn]*da/2) 0;0 exp(-im*ka[nn]*da/2)];
            T=v*inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
        else
            m=[exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2);ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
            n=[1 1;kb[nn] -kb[nn]];
            o=[exp(im*kb[nn]*db) exp(-im*kb[nn]*db);kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
            p=[1 1;ka[nn] -ka[nn]];
            q=[exp(im*ka[nn]*(dx[mm]-da/2-db)) exp(-im*ka[nn]*(dx[mm]-da/2-db));ka[nn]*exp(im*ka[nn]*(dx[mm]-da/2-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-da/2-db))];
            r=[1 1;kc[nn] -kc[nn]];
            s=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
            t=[1 1;kc[nn] -kc[nn]];
            u=[1 1;ka[nn] -ka[nn]];
            v=[exp(im*ka[nn]*(D-dx[mm])) 0;0 exp(-im*ka[nn]*(D-dx[mm]))];
            T=v*inv(u)*t*s*inv(r)*q*inv(p)*o*inv(n)*m;
        end
        AAA=T[1,1]+T[2,2];
        Temp0=AAA^2-4;
        Temp1=(AAA-sqrt(Temp0))/2;
        r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
        if r0[mm,nn]<=1#导带取一种情况
            Temp1=(AAA-sqrt(Temp0))/2;
            r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
            Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
        else#此时并不都是禁带，r0仍可能为1，所以取出r0等于1的值与上面导带的情况简并
            Temp1=(AAA+sqrt(Temp0))/2;
            r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
            if r0[mm,nn]<=1.0001&&r0[mm,nn]>=0.9999#1为int所以加一些小数位
                Temp1=(AAA-sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
            else
                Temp1=(AAA+sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
            end
        end
        R[mm,nn]=real(T[1,2]);
        if R[mm,nn]<=0.01&&R[mm,nn]>=0.00
            R[mm,nn]=real(T[1,2]);
        else
            R[mm,nn]=1;
        end
    end
end

function draw_figure(dx0,w,Arg,r0)
    f11 = heatmap(dx0,w,Arg',color=:bluesreds,dpi=600,
    title=L"Arg(r)",
    ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
    xlabel = L"\Delta (nm)")

    f12 = heatmap(dx0,w,r0',color=:bluesreds,dpi=600,
    title=L"|r|",
    xlabel = L"\Delta (nm)")
    plot(f11,f12, layout=(1,2))
end

function shaixuank(k1,da,db)
    k=zeros(length(k1));
    for ii=1:length(k1)
        if k1[ii]<=da/2+db&&k1[ii]>=da/2
            k[ii]=k1[ii];
        else
            k[ii]=NaN;
        end
    end
    return k;
end
function shaixuanp1(p1,da,db)
    p=zeros(length(p1));
     for ii=1:length(p1)
        if p1[ii]<=da/2&&p1[ii]>=0
            p[ii]=p1[ii];
        else
            p[ii]=NaN;
        end
    end
    return p;
end
function shaixuanp2(p1,da,db)
    p=zeros(length(p1));
     for ii=1:length(p1)
        if p1[ii]>da/2+db&&p1[ii]<da+db
            p[ii]=p1[ii];
        else
            p[ii]=NaN;
        end
    end
    return p;
end

#draw_figure(dx0,w,Arg,r0)

#heatmap(dx0,w,R',dpi=600)
#savefig("pcls.png")
#plot(w,R[50,:])
#ylims!(-10,10)
#savefig("3p2lc=da=db=100_na=3p2_nb=2_nc=1.png")
#图片大小：pdf>>ps>svg>png
#latex支持pdf,png,eps
#gr支持输出pdf,ps,svg,png
#gr🐮🍺
k1=(((pi./kb.+db)/2).+da/2);
k1=shaixuank(k1,da,db);
k_1=(((-pi./kb.+db)/2).+da/2);
k_1=shaixuank(k_1,da,db);
k0=ones(num1)*(D/2);

k_2=(((-2*pi./kb.+db)/2).+da/2);
k_2=shaixuank(k_2,da,db);
k2=(((2*pi./kb.+db)/2).+da/2);
k2=shaixuank(k2,da,db);

c1=atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c2=atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c3=atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c4=atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
p110=(ka.*da.-(c3))./(2)./ka;
p111=(ka.*da.-(c3).-pi)./(2)./ka;
p112=(ka.*da.-(c3).-(2).*pi)./(2)./ka;
p113=(ka.*da.-(c3).-(3).*pi)./(2)./ka;
p114=(ka.*da.-(c3).-(4).*pi)./(2)./ka;

p120=(ka.*da.-(c4))./(2)./ka;
p121=(ka.*da.-(c4).-pi)./(2)./ka;
p122=(ka.*da.-(c4).-(2).*pi)./(2)./ka;
p123=(ka.*da.-(c4).-(3).*pi)./(2)./ka;
p124=(ka.*da.-(c4).-(4).*pi)./(2)./ka;

p220=(ka.*da.+(2).*ka.*db.-(c2))./ka./(2);
p221=(ka.*da.+(2).*ka.*db.-(c2).+pi)./ka./(2);
p222=(ka.*da.+(2).*ka.*db.-(c2).+(2).*pi)./ka./(2);
p223=(ka.*da.+(2).*ka.*db.-(c2).+(3).*pi)./ka./(2);
p224=(ka.*da.+(2).*ka.*db.-(c2).+(4).*pi)./ka./(2);

p210=(ka.*da.+(2).*ka.*db.-(c1))./ka./(2);
p211=(ka.*da.+(2).*ka.*db.-(c1).+pi)./ka./(2);
p212=(ka.*da.+(2).*ka.*db.-(c1).+(2).*pi)./ka./(2);
p213=(ka.*da.+(2).*ka.*db.-(c1).+(3).*pi)./ka./(2);
p214=(ka.*da.+(2).*ka.*db.-(c1).+(4).*pi)./ka./(2);

p110=shaixuanp1(p110,da,db);p111=shaixuanp1(p111,da,db);p112=shaixuanp1(p112,da,db);p113=shaixuanp1(p113,da,db);p114=shaixuanp1(p114,da,db);#p115=shaixuanp1(p115);
p120=shaixuanp1(p120,da,db);p121=shaixuanp1(p121,da,db);p122=shaixuanp1(p122,da,db);p123=shaixuanp1(p123,da,db);p124=shaixuanp1(p124,da,db);#p125=shaixuanp1(p125);
p210=shaixuanp2(p210,da,db);p211=shaixuanp2(p211,da,db);p212=shaixuanp2(p212,da,db);p213=shaixuanp2(p213,da,db);p214=shaixuanp2(p214,da,db);#p215=shaixuanp2(p215);
p220=shaixuanp2(p220,da,db);p221=shaixuanp2(p221,da,db);p222=shaixuanp2(p222,da,db);p223=shaixuanp2(p223,da,db);p224=shaixuanp2(p224,da,db);#p225=shaixuanp2(p225);

dataA1=[p110,p112,p114,p120,p122,p124];
dataA2=[p210,p212,p214,p220,p222,p224];
dataB=[k0,k1,k_1,k2,k_2]
data=[dataA1,dataA2,dataB]*1e9
#plot(data,w,legend=:none)
#ylim!([0,D*1e9])
heatmap(dx0,w,Arg',color=:bluesreds,dpi=600,
title=L"Arg(r)",
ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
xlabel = L"\Delta (nm)")
plot!(data,w,legend=:none,color=:black,line=(:dot,1.5))
#ylim!([0,D*1e9])

heatmap(dx0,w,r0',color=:bluesreds,dpi=600,
title=L"|r|",
ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
xlabel = L"\Delta (nm)")
plot!(data,w,legend=:none,color=:black,line=(:dot,1.5))
#ylim!([0,D*1e9])
#f1=plot(f11,f12,layout=(1,1))
#f21 = heatmap(dx0,w,r0',color=:bluesreds,dpi=600,
#title=L"|r|",
#xlabel = L"\Delta (nm)")
#f2=plot(f21,f12)
#plot(f1,f2, layout=(1,2))
