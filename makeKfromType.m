function K = makeKfromType(sizeOfK,elemTypeNo,elements,nodes,mats,elmats)

K=zeros(sizeOfK,sizeOfK);
if(elemTypeNo ==1)      %for a linear bar element
    c=elemA*globalE/elemL;
    k=c*[1 -1;
        -1 1;];
    noOfElem = size(elements,1);
    for i=1:noOfElem
        K(i:i+1,i:i+1)=K(i:i+1,i:i+1)+k;    %assembly
    end
elseif(elemTypeNo==4)% CPS3
    noOfElem = size(elements,1);
    for i=1:noOfElem
        %get material for this element
        E=mats(elmats(i),1);
        v=mats(elmats(i),2);
        %make property matriex
        C=E/(1-v^2)*[   1 v 0;
                        v 1 0;
                        0 0 (1-v)/2;];
        %get element node locations
        index = binarySearch(nodes,1,sizeOfK,elements(i,2));
        elements(i,2)=index;    %relabel this node number
        x1=nodes(index,2);
        y1=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,3));
        elements(i,3)=index;
        x2=nodes(index,2);
        y2=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,4));
        elements(i,4)=index;
        x3=nodes(index,2);
        y3=nodes(index,3);
        %make J and B matricies
        J= [x1-x3 y1-y3;
            x2-x3 y2-y3];
        B=1/det(J)*[ y2-y3 0 y3-y1 0 y1-y2 0;
            0 x3-x2 0 x1-x3 0 x2-x1;
            x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
        A=0.5*abs(det(J));
        %make element stiffness matrix
        kelem = A*transpose(B)*C*B;
        %assembly into global k
        for j=1:3
            %nodes relabled to their indicies
            x=elements(i,j+1)*2-1;            %x location on global K
            a=j*2-1;                %x location on local k
            for k=1:3  
                y=elements(i,k+1)*2-1;            %y location on global K
                b=k*2-1;                %y location on local k
                K(x,y)=kelem(a,b)+K(x,y);
                K(x+1,y)=kelem(a+1,b)+K(x+1,y);
                K(x,y+1)=kelem(a,b+1)+K(x,y+1);
                K(x+1,y+1)=kelem(a+1,b+1)+K(x+1,y+1);
            end
        end
    end
elseif(elemTypeNo==5)% CPE3
    noOfElem = size(elements,1);
    for i=1:noOfElem
        E=mats(elmats(i),1);
        v=mats(elmats(i),2);
        C=E/(1+v)/(1-2*v)*[ 1-v v 0;
                            v 1-v 0;
                            0 0 .5-v;];
        index = binarySearch(nodes,1,sizeOfK,elements(i,2));
        elements(i,2)=index;        %relable the node number
        x1=nodes(index,2);
        y1=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,3));
        elements(i,3)=index;        %relable the node number
        x2=nodes(index,2);
        y2=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,4));
        elements(i,4)=index;        %relable the node number
        x3=nodes(index,2);
        y3=nodes(index,3);
        J= [x1-x3 y1-y3;
            x2-x3 y2-y3];
        B=1/det(J)*[ y2-y3 0 y3-y1 0 y1-y2 0;
            0 x3-x2 0 x1-x3 0 x2-x1;
            x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
        A=0.5*abs(det(J));
        kelem = A*transpose(B)*C*B;
        %assembly into global k
        for j=1:3
            x=elements(i,j+1)*2-1;
            a=j*2-1;
            for k=1:3
                y=elements(i,k+1)*2-1;
                b=k*2-1;
                
                K(x,y)=kelem(a,b)+K(x,y);
                K(x+1,y)=kelem(a+1,b)+K(x+1,y);
                K(x,y+1)=kelem(a,b+1)+K(x,y+1);
                K(x+1,y+1)=kelem(a+1,b+1)+K(x+1,y+1);
            end
        end
    end
elseif(elemTypeNo==6)% CPS4
    noOfElem = size(elements,1);
    for i=1:noOfElem
        %this is inside loop since each element can have diff mat
        E=mats(elmats(i),1);
        v=mats(elmats(i),2);
        C=E/(1-v^2)*[   1 v 0;
                        v 1 0;
                        0 0 (1-v)/2;];
        %find element coordinates
        index = binarySearch(nodes,1,sizeOfK,elements(i,2));
        elements(i,2)=index;        %relable the node number
        x1=nodes(index,2);
        y1=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,3));
        elements(i,3)=index;        %relable the node number
        x2=nodes(index,2);
        y2=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,4));
        elements(i,4)=index;        %relable the node number
        x3=nodes(index,2);
        y3=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,5));
        elements(i,5)=index;        %relable the node number
        x4=nodes(index,2);
        y4=nodes(index,3);
        %gauss points
        gp1=-.57735;
        gp2=-gp1;
        %gauss point 22
        J=JforQ4(gp2, gp2,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp2,gp2);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B;
        %gauss point 21
        J=JforQ4(gp2, gp1,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp2,gp1);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %gauss point 12
        J=JforQ4(gp1, gp2,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp1,gp2);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %gauss point 11
        J=JforQ4(gp1, gp1,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp1,gp1);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %assembly into global k
        for j=1:4
            %remember elements where renumbered to their indexies
            x=elements(i,j+1)*2-1;
            a=j*2-1;
            for k=1:4
                y=elements(i,k+1)*2-1;
                b=k*2-1;
                K(x,y)=kelem(a,b)+K(x,y);
                K(x+1,y)=kelem(a+1,b)+K(x+1,y);
                K(x,y+1)=kelem(a,b+1)+K(x,y+1);
                K(x+1,y+1)=kelem(a+1,b+1)+K(x+1,y+1);
            end
        end
    end
elseif(elemTypeNo==7)% CPE4
    noOfElem = size(elements,1);
    for i=1:noOfElem
        E=mats(elmats(i),1);
        v=mats(elmats(i),2);
        C=E/(1+v)/(1-2*v)*[ 1-v v 0;
                            v 1-v 0;
                            0 0 .5-v;];
        %find element coordinates
        index = binarySearch(nodes,1,sizeOfK,elements(i,2));
        elements(i,2)=index;        %relable the node number
        x1=nodes(index,2);
        y1=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,3));
        elements(i,3)=index;        %relable the node number
        x2=nodes(index,2);
        y2=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,4));
        elements(i,4)=index;        %relable the node number
        x3=nodes(index,2);
        y3=nodes(index,3);
        index = binarySearch(nodes,1,sizeOfK,elements(i,5));
        elements(i,5)=index;        %relable the node number
        x4=nodes(index,2);
        y4=nodes(index,3);
        %gauss points
        gp1=-.57735;
        gp2=-gp1;
        %gauss point 22
        J=JforQ4(gp2, gp2,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp2,gp2);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B;
        %gauss point 21
        J=JforQ4(gp2, gp1,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp2,gp1);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %gauss point 12
        J=JforQ4(gp1, gp2,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp1,gp2);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %gauss point 11
        J=JforQ4(gp1, gp1,x1,y1,x2,y2,x3,y3,x4,y4);
        G=GforQ4(gp1,gp1);
        A=AforQ4(J);
        B=A*G;
        kelem=det(J)*transpose(B)*C*B+kelem;
        %assembly into global k
        for j=1:4
            x=elements(i,j+1)*2-1;
            a=j*2-1;
            for k=1:4
                y=elements(i,k+1)*2-1;
                b=k*2-1;
                K(x,y)=kelem(a,b)+K(x,y);
                K(x+1,y)=kelem(a+1,b)+K(x+1,y);
                K(x,y+1)=kelem(a,b+1)+K(x,y+1);
                K(x+1,y+1)=kelem(a+1,b+1)+K(x+1,y+1);
            end
        end
    end
end
end


function J=JforQ4(eta, exi,x1,y1,x2,y2,x3,y3,x4,y4)
    j11=-(1-eta)*x1+(1-eta)*x2+(1+eta)*x3-(1+eta)*x4;
    j21=-(1-exi)*x1-(1+exi)*x2+(1+exi)*x3+(1-exi)*x4;
    j12=-(1-eta)*y1+(1-eta)*y2+(1+eta)*y3-(1+eta)*y4;
    j22=-(1-exi)*y1-(1+exi)*y2+(1+exi)*y3+(1-exi)*y4;
    J=1/4*[j11 j12;
        j21 j22;];
end

function G=GforQ4(eta, exi)
    G=1/4*[ eta-1   0       1-eta   0     1+eta 0    -1-eta     0 ;
            exi-1   0       -exi-1  0     1+exi 0    1-exi      0 ;
            0       eta-1   0     1-eta   0     1+eta 0     -eta-1;
            0       exi-1   0     -1-exi  0     1+exi 0     1-exi;];
end

function A=AforQ4(J)
   j11=J(1,1);
   j12=J(1,2);
   j21=J(2,1);
   j22=J(2,2);
   A=1/det(J)*[j22 -j12 0 0;
                0 0 -j21 j11;
                -j21 j11 j22 -j12;];
end