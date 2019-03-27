function onecell(k,h0,w,Delta,h1,h2)
    Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
    num1=length(w);
    num2=length(Delta);
    r=zeros(num2,num1);Arg=zeros(num2,num1);
    for m=1:num2
        for n=1:num1
            Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
            Th1=[exp(1j*k(n)*h1(m)),0;0,exp(-1j*k(n)*h1(m))];
            Th2=[exp(1j*k(n)*h2(m)),0;0,exp(-1j*k(n)*h2(m))];
            T=Th1*Tc*Th2*Tc*Th1;
            r(m,n)=abs(-T(2,1)./T(2,2));
            Arg(m,n)=angle(-T(2,1)./T(2,2));
        end
    end
    figure
    imagesc(Delta,w,r);
    figure
    imagesc(Delta,w,Arg);
end