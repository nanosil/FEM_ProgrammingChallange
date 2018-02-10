function graphMesh(nodes, elements,eltype)
%set up and initialize
sizeOfNodes=size(nodes,1);
sizeOfElements=size(elements,1);
AdjencyMatrix=zeros(sizeOfNodes,sizeOfNodes);

x=eltype(size(eltype,2));
x = str2double(x);
if(x==3)
    for i=1:sizeOfElements
        index1 = binarySearch(nodes,1,sizeOfNodes,elements(i,2));
        index2 = binarySearch(nodes,1,sizeOfNodes,elements(i,3));
        index3 = binarySearch(nodes,1,sizeOfNodes,elements(i,4));
        AdjencyMatrix(index1,index2)=1;
        AdjencyMatrix(index2,index3)=1;
        AdjencyMatrix(index3,index1)=1;
    end
elseif(x==4)
    for i=1:sizeOfElements
        index1 = binarySearch(nodes,1,sizeOfNodes,elements(i,2));
        index2 = binarySearch(nodes,1,sizeOfNodes,elements(i,3));
        index3 = binarySearch(nodes,1,sizeOfNodes,elements(i,4));
        index4 = binarySearch(nodes,1,sizeOfNodes,elements(i,5));
        AdjencyMatrix(index1,index2)=1;
        AdjencyMatrix(index2,index3)=1;
        AdjencyMatrix(index3,index4)=1;
        AdjencyMatrix(index4,index1)=1;
    end
end

figure
gplot(AdjencyMatrix,nodes(:,2:3)) %only plot in 2D
axis square
end