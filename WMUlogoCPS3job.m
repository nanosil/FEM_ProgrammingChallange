%clear variables
clear
clear,clc

%load mesh
WMUlogoCPS3mesh

%combine BCs and Elements
fixed = [fixedx ; fixedtop  ; fixedbot];
elmats = [ones(size(elements1,1),1); 2*ones(size(elements2,1),1)];
elements = [elements1; elements2];
forces = [];

%define element type
eltype='CPS3';

%define elemment materials
mat1 = [200e3 .25 7000e-12];
mat2 = [70e3 .33 2700e-12];
mats = [ mat1; mat2];

%free space
clear fixedbot fixedtop fixedx mat1 mat2
%uncomment for testing geometry!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%testMeshing
%nodes=sortrows(nodes);              %nodes must be sorted!!!!!!!!!!
graphMesh(nodes, elements,eltype);
tic
q=runjob(nodes,elements,fixed,eltype,elmats,mats,forces);
toc
graphResult(q,nodes,elements,eltype,100);
