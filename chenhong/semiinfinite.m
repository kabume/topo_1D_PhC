function semiinfinite(k,h0,w,Delta,h1,h2)
    %transfer matrix of Tc
    Z0=377;
    Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
    P=[1,1;1/Z0,-1/Z0];
    num1=length(w);
    num2=length(Delta);
    r=zeros(num2,num1);Arg=zeros(num2,num1);
    r1=zeros(num2,num1);Arg1=zeros(num2,num1);
    parfor n=1:num1
        for m=1:num2
            Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
            Th1=[exp(1j*k(n)*h1(m)),0;0,exp(-1j*k(n)*h1(m))];
            Th2=[exp(1j*k(n)*h2(m)),0;0,exp(-1j*k(n)*h2(m))];
            T=Th1*Tc*Th2*Tc*Th1;
            M=P*T*P^-1;
            Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            if real(Z)>0
                 Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            else
                Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            end
            r1(m,n)=abs((Z-Z0)/(Z+Z0));
            Arg1(m,n)=angle((Z-Z0)/(Z+Z0));
            if Arg1(m,n)<-pi+0.001&&Arg1(m,n)>-pi-0.001%·ÀÒç³ö
                Arg1(m,n)=pi;
            end
            r(m,n)=abs(-T(2,1)./T(2,2));
            Arg(m,n)=angle(-T(2,1)./T(2,2));
        end
    end
    figure
    subplot(2,2,1)
    imagesc(w,Delta,r);
    title("|r| of one cell")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,2)
    imagesc(w,Delta,Arg);
    title("Arg of one cell")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,3)
    imagesc(w,Delta,r1);
    title("|r| of semi-infinite systems")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,4)
    imagesc(w,Delta,Arg1);
    title("Arg of semi-infinite systems")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
%     figure
%     imagesc(w,Delta,Arg1);
%     title("Arg of semi-infinite systems")
%     
end
