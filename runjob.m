function q=runjob(nodes,elements,BCs,eltype,elmats,mats,forces)

%initial calculations----------------------------------------------------
[elemTypeNo,~,noDOFperNode] = ElemProp(eltype);%find elem props
noOfNodes = size(nodes,1);          %calc number of global nodes
sizeOfK= noOfNodes*noDOFperNode;    %find size of global K


%make global K and global F matrixies------------------------------------
K = makeKfromType(sizeOfK,elemTypeNo,elements,nodes,mats,elmats);
F = makeFfromf(sizeOfK,forces,noDOFperNode);

%penalty method  for boundry conditions---------------------------------
[K, F] = makeBC(K,F,BCs,noDOFperNode,nodes);

%solving----------------------------------------------------------------
q = mldivide(K,F);
%q = gmres(K,F,[],.0001);
%x0=zeros(sizeOfK,1);
%q = jacobi2(K,F,x0,.001,2);

%[L,U] = lu(K);
%Y= L \ F;
%q = U \ Y;
end

