const c = 3e8;

Arg=nothing;dx=nothing;r0=nothing;w=nothing;data=nothing;f1=nothing;f2=nothing;
da=100e-9;db=200e-9;lc=da;#单位m
na=2;#slab A
nb=1;#slab B
nc=2;
w=0.1:0.01:15;
dx=0:0.001*1e-9:da+db;

num1=length(w);
D=da+db;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;

bj=[-2*pi,-pi,pi,2*pi];
Bn=SharedArray{Float64}(length(bj), length(w));
rr0 = CartesianIndices(size(Bn));
@sync begin#同步
    @distributed for k in rr0
        mm,nn=k.I;
        Bn[mm,nn]=(((bj[mm]/kb[nn]+db)/2)+da/2);
        if Bn[mm,nn]<=da/2+db&&Bn[mm,nn]>=da/2
            Bn[mm,nn]=Bn[mm,nn];
        else
            Bn[mm,nn]=NaN;
        end
    end
end
k_2=Bn[1,1:num1];k_1=Bn[2,1:num1];k1=Bn[3,1:num1];k2=Bn[4,1:num1];k0=ones(num1)*(D/2);

c1=atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c2=atan.(.-(tan.(kb .* db) .* ka .^ (2) .+ tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (tan.(kb .* db) .* ka .^ (2) .- tan.(kb .* db) .* kb .^ (2)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c3=atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .- (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));
c4=atan.(.-(.-ka .^ (2) .* tan.(kb .* db) .- kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb .* (ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)) ./ (2) .+ (.-ka .^ (2) .* tan.(kb .* db) .+ kb .^ (2) .* tan.(kb .* db)) ./ ka ./ kb ./ (2), .-(ka .^ (4) .* tan.(kb .* db) .^ (2) .- kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* sqrt.(tan.(kb .* db) .^ (2) .* ka .^ (4) .* kb .^ (4) .+ ka .^ (4) .* kb .^ (4))) ./ (ka .^ (4) .* tan.(kb .* db) .^ (2) .+ (2) .* ka .^ (2) .* kb .^ (2) .* tan.(kb .* db) .^ (2) .+ kb .^ (4) .* tan.(kb .* db) .^ (2) .+ (4) .* ka .^ (2) .* kb .^ (2)));


aj1=[0,pi,2*pi,3*pi,4*pi];
An11=SharedArray{Float64}(length(aj1), length(w));
An12=SharedArray{Float64}(length(aj1), length(w));
rr1 = CartesianIndices(size(An11));
@sync begin#同步
    @distributed for k in rr1
        mm,nn=k.I;
        An11[mm,nn]=(ka[nn]*da-c3[nn]-aj1[mm])/2/ka[nn];
        An12[mm,nn]=(ka[nn]*da-c4[nn]-aj1[mm])/2/ka[nn];
        if An11[mm,nn]<=da/2&&An11[mm,nn]>=0
            An11[mm,nn]=An11[mm,nn];
        else
            An11[mm,nn]=NaN;
        end
        if An12[mm,nn]<=da/2&&An12[mm,nn]>=0
            An12[mm,nn]=An12[mm,nn];
        else
            An12[mm,nn]=NaN;
        end
    end
end
p110=An11[1,1:num1];p111=An11[2,1:num1];p112=An11[3,1:num1];p113=An11[4,1:num1];p114=An11[5,1:num1];
p120=An12[1,1:num1];p121=An12[2,1:num1];p122=An12[3,1:num1];p123=An12[4,1:num1];p124=An12[5,1:num1];

aj2=[0,pi,2*pi,3*pi,4*pi];
An21=SharedArray{Float64}(length(aj2), length(w));
An22=SharedArray{Float64}(length(aj2), length(w));
rr2 = CartesianIndices(size(An21));
@sync begin#同步
    @distributed for k in rr2
        mm,nn=k.I;
        An21[mm,nn]=(ka[nn]*da+2*ka[nn]*db-c1[nn]+aj2[mm])/ka[nn]/2;
        An22[mm,nn]=(ka[nn]*da+2*ka[nn]*db-c2[nn]+aj2[mm])/ka[nn]/2;
        if An21[mm,nn]>da/2+db&&An21[mm,nn]<da+db
            An21[mm,nn]=An21[mm,nn];
        else
            An21[mm,nn]=NaN;
        end
        if An22[mm,nn]>da/2+db&&An22[mm,nn]<da+db
            An22[mm,nn]=An22[mm,nn];
        else
            An22[mm,nn]=NaN;
        end
    end
end

p210=An21[1,1:num1];p211=An21[2,1:num1];p212=An21[3,1:num1];p213=An21[4,1:num1];p214=An21[5,1:num1];
p220=An22[1,1:num1];p221=An22[2,1:num1];p222=An22[3,1:num1];p223=An22[4,1:num1];p224=An22[5,1:num1];

dataA1=[p110,p112,p114,p120,p122,p124];
dataA2=[p210,p212,p214,p220,p222,p224];
dataB=[k0,k1,k_1,k2,k_2];
data=[dataA1,dataA2,dataB]*1e9;
