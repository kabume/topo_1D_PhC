using Plots
using Unicode
using LaTeXStrings
using SharedArrays
using Distributed
gr()
#pyplot()
Arg=nothing;dx=nothing;dx0=nothing;ka=nothing;kb=nothing;kc=nothing;p11=nothing;p12=nothing;r0=nothing;rr=nothing;w=nothing;

#addprocs(72)#除了第一次运行，下次运行的时候要注释掉

da=100e-9;db=100e-9;D=da+db;lc=da;#单位m
Dt=D*1e9;
dx0=0:0.1:Dt;
dx=dx0*1e-9
na=3.2;#slab A
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
R=SharedArray{Float64}(num2, num1);
rr = CartesianIndices(size(r0));
@sync begin#同步
    @distributed for k in rr
        mm,nn=k.I;
        if mm<db/D*num2
            #=
            m = [[exp(im*ka[nn]*da) exp(-im*ka[nn]*da)];[ka[nn]*exp(im*ka[nn]*da) -ka[nn]*exp(-im*ka[nn]*da)]];
            n = [[1 1];[kb[nn] -kb[nn]]];#[exp(im*kb[nn]*da),exp(-im*kb[nn]*da);kb[nn]*exp(im*kb[nn]*da),-kb[nn]*exp(-im*kb[nn]*da)];
            o = [[1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2];[-im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2]] * [[exp(im*kb[nn]*(db-dx[mm])) 0];[0 exp(-im*kb[nn]*(db-dx[mm]))]];
            p = [1 0; 0 1];
            q = [[exp(im*kb[nn]*dx[mm]) exp(-im*kb[nn]*dx[mm])];[kb[nn]*exp(im*kb[nn]*dx[mm]) -kb[nn]*exp(-im*kb[nn]*dx[mm])]];
            r = [[1 1];[ka[nn] -ka[nn]]];
            T = inv(r) * q * inv(p) * o * inv(n) * m;
            =#
            m=[[exp(im*ka[nn]*da) exp(-im*ka[nn]*da)];[ka[nn]*exp(im*ka[nn]*da) -ka[nn]*exp(-im*ka[nn]*da)]];
            n=[[1 1];[kb[nn] -kb[nn]]];#[exp(im*kb[nn]*da),exp(-im*kb[nn]*da);kb[nn]*exp(im*kb[nn]*da),-kb[nn]*exp(-im*kb[nn]*da)];
            o=[[exp(im*kb[nn]*(db-dx[mm])) exp(-im*kb[nn]*(db-dx[mm]))];[kb[nn]*exp(im*kb[nn]*(db-dx[mm])) -kb[nn]*exp(-im*kb[nn]*(db-dx[mm]))]];
            p=[[1 1];[kc[nn] -kc[nn]]];
            q=[[1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2];[-im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2]];
            r=[[1 1];[kc[nn] -kc[nn]]];
            s=[[1 1];[kb[nn] -kb[nn]]];
            t=[[exp(im*kb[nn]*dx[mm]) exp(-im*kb[nn]*dx[mm])];[kb[nn]*exp(im*kb[nn]*dx[mm]) -kb[nn]*exp(-im*kb[nn]*dx[mm])]];
            u=[[1 1];[ka[nn] -ka[nn]]];
            T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
        #=
        elseif mm>db/D*num2+5 && mm<num2

            m = [[1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2];[-im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2]] * [[exp(im*ka[nn]*(da+db-dx[mm])) 0];[0 exp(-im*ka[nn]*(da+db-dx[mm]))]];
            n = [1 0; 0 1];
            o = [[exp(im*ka[nn]*(dx[mm]-db)) exp(-im*ka[nn]*(dx[mm]-db))];[ka[nn]*exp(im*ka[nn]*(dx[mm]-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-db))]];
            p = [[1 1];[kb[nn] -kb[nn]]];
            q = [[exp(im*kb[nn]*db) exp(-im*kb[nn]*db)];[kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)]];
            r = [[1 1];[ka[nn] -ka[nn]]];
            T = inv(r) * q *inv(p) * o * inv(n) * m;
        =#
        else
            m=[[exp(im*ka[nn]*(da+db-dx[mm])) exp(-im*ka[nn]*(da+db-dx[mm]))];[ka[nn]*exp(im*ka[nn]*(da+db-dx[mm])) -ka[nn]*exp(-im*ka[nn]*(da+db-dx[mm]))]];
            n=[[1 1];[kc[nn] -kc[nn]]];
            o=[[1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2];[-im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2]];
            p=[[1 1];[kc[nn] -kc[nn]]];
            q=[[1 1];[ka[nn] -ka[nn]]];
            r=[[exp(im*ka[nn]*(dx[mm]-db)) exp(-im*ka[nn]*(dx[mm]-db))];[ka[nn]*exp(im*ka[nn]*(dx[mm]-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-db))]];
            s=[[1 1];[kb[nn] -kb[nn]]];
            t=[[exp(im*kb[nn]*db) exp(-im*kb[nn]*db)];[kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)]];
            u=[[1 1];[ka[nn] -ka[nn]]];
            T=inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
            #=
            m = [[exp(im*ka[nn]*da) exp(-im*ka[nn]*da)];[ka[nn]*exp(im*ka[nn]*da) -ka[nn]*exp(-im*ka[nn]*da)]];
            n = [[1 1];[kc[nn] -kc[nn]]];
            o = [[1+im*tan(kc[nn]*lc)/2 im*tan(kc[nn]*lc)/2];[-im*tan(kc[nn]*lc)/2 1-im*tan(kc[nn]*lc)/2]];
            p = [[1 1];[kc[nn] -kc[nn]]];
            r = [[1 1];[kb[nn] -kb[nn]]];
            s = [[exp(im*kb[nn]*db) exp(-im*kb[nn]*db)];[kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)]];
            t = [[1 1];[ka[nn] -ka[nn]]];
            T = inv(t) * s * inv(r) * p *o *inv(n) * m;
            =#
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
    end
end
function draw_figure(dx0,w,Arg,r0)
    p11 = heatmap(dx0,w,Arg',color=:bluesreds,dpi=600,
    title=L"Arg(r)",
    ylabel = L"f\, (10^{15}/ 2 \pi\ Hz)",
    xlabel = L"\Delta (nm)")

    p12 = heatmap(dx0,w,r0',color=:bluesreds,dpi=600,
    title=L"|r|",
    xlabel = L"\Delta (nm)")
    plot(p11,p12, layout=(1,2))
end
draw_figure(dx0,w,Arg,r0)
#savefig("第二种传输矩阵.png")
#savefig("3p2lc=da=db=100_na=3p2_nb=2_nc=1.png")
#图片大小：pdf>>ps>svg>png
#latex支持pdf,png,eps
#gr支持输出pdf,ps,svg,png
#gr🐮🍺


#savefig("lc=da相位.pdf")
#heatmap(dx,w,Em')
