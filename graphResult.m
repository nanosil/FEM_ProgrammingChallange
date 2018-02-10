function graphResult(q,nodes,elements,eltype,exaggeration)
sizeofNodes=size(nodes,1);
for i=1:sizeofNodes
    nodes(i,2)=q(i*2-1)*exaggeration+nodes(i,2);
    nodes(i,3)=q(i*2)*exaggeration+ nodes(i,3);
end
graphMesh(nodes, elements,eltype)
end