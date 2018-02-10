function [elemTypeNo, noOfNodes, noDOFperNode] = ElemProp(elemType)

noOfNodes=-1;
elemTypeNo =-1;
noDOFperNode=-1;

if(strcmp(elemType,'Linear Bar')==1)
    elemTypeNo=1;
    noOfNodes = 2;
    noDOFperNode=1;
elseif(strcmp(elemType,'CST')==1)
    elemTypeNo=2;
    noOfNodes = 3;
    noDOFperNode=2;
elseif(strcmp(elemType,'Q4')==1)
    elemTypeNo=3;
    noOfNodes = 4;
    noDOFperNode=2;
elseif(strcmp(elemType,'CPS3')==1)
    elemTypeNo=4;
    noOfNodes = 3;
    noDOFperNode=2;
elseif(strcmp(elemType,'CPE3')==1)
    elemTypeNo=5;
    noOfNodes = 3;
    noDOFperNode=2;
elseif(strcmp(elemType,'CPS4')==1)
    elemTypeNo=6;
    noOfNodes = 4;
    noDOFperNode=2;
elseif(strcmp(elemType,'CPE4')==1)
    elemTypeNo=7;
    noOfNodes = 4;
    noDOFperNode=2;
end
end
