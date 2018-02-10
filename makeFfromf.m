function F = makeFfromf(noOfNodes,forces,DOFperNode)

F= zeros(noOfNodes,1);  %preAllocate

sizeOfF= size(forces,1);

for i=1:sizeOfF
   F(forces(i,1)*DOFperNode+forces(i,2)-2)= forces(i,3); 
end


end