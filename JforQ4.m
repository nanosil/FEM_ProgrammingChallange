
function J=JforQ4(eta, exi,x1,y1,x2,y2,x3,y3,x4,y4)
    
    j11=-(1-eta)*x1+(1-eta)*x2+(1+eta)*x3-(1+eta)*x4;
    j21=-(1-exi)*x1-(1+exi)*x2+(1+exi)*x3+(1-exi)*x4;
    
    j12=-(1-eta)*y1+(1-eta)*y2+(1+eta)*y3-(1+eta)*y4;
    j22=-(1-exi)*y1-(1+exi)*y2+(1+exi)*y3+(1-exi)*y4;
    
    J=[j11 j12;
        j21 j22;];
    
    
end