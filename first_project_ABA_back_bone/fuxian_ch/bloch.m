%Calculate the band diagram 
function bloch(w,k,D,h0,Delta)
    A=cos(k*D)-(tan(k*h0).^2).*cos(k*D)/4-tan(k*h0).*sin(k*D)+(tan(k*h0).^2).*cos(k*Delta*D)/4;%cosqD
    B=real(acos(A));
    plot(B,w,'k');
    hold on
    plot(-B,w,'k');
    hold off
end