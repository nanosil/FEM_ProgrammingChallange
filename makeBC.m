function [K,F] = makeBC(K,F,BCs,noDOFperNode,nodes)

c = max(max(K))*10^4;

sizeOfBC = size(BCs,1);
sizeOfN = size(nodes,1);
for i=1:sizeOfBC
    index = binarySearch(nodes,1,sizeOfN,BCs(i,1)); %index of node
    x=(index)*noDOFperNode+BCs(i,2)-2;              %index of node in K
    K(x,x) =  K(x,x)+c;                             %add stiffness
    F(x,1)= F(x,1)+BCs(i,3)*c;                      %add to force
end

end