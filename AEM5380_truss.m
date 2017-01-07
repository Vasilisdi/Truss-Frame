%%
% HA assignment in FEM - part A
% vasileios dimitriou 5380
% Finite element MATLAB program for truss analysis
% real construction of truss with 8 constr. 
%%
clear all;clc;close all;
NodesandElements = load('AEM5380_Geom_truss.txt');  %x, y coordinates & first node, second node, E, A
numnode = NodesandElements(1,1);             %Holds number of nodes
numelem = NodesandElements(1,2);             %Holds number of elements
Nodes=NodesandElements(2:numnode+1,:); Elements=NodesandElements(numnode+3:numnode+3+numelem-1,:);  Constraints=NodesandElements(length(NodesandElements)-(NodesandElements(1,4)-1):length(NodesandElements),:); Forces=NodesandElements(length(NodesandElements)-(NodesandElements(1,4)+NodesandElements(1,3)):length(NodesandElements)-(NodesandElements(1,4)+1),:);  
K = zeros(3*numnode,3*numnode);   strain=zeros(numelem,1);   stress=zeros(numelem,1);   axialforce=zeros(numelem,1);  constr=zeros(1,NodesandElements(length(NodesandElements)-NodesandElements(1,4)));        
F = zeros(3*numnode, 1);                     %gen. forces
U = zeros(3*numnode, 1);                     %F:OXYZ
for ie=1:length(Constraints(:,[1 2]))
 constr(ie)=3*(Constraints(ie,1)-1)+Constraints(ie,2);
end
act = 1:3*numnode;                           %Holds active DOFs
act(constr) = []; ns=zeros(numelem,1); ms=zeros(numelem,1); ls=zeros(numelem,1); le=zeros(numelem,1);
for ie = 1:numelem
DOFs = [3*Elements(ie, 2)-2, 3*Elements(ie, 2)-1, 3*Elements(ie, 2), 3*Elements(ie, 3)-2, 3*Elements(ie, 3)-1, 3*Elements(ie, 3) ];   %Holds element’s DOFs for an elem
X1 = Nodes(Elements(ie,2), 2);  Y1 = Nodes(Elements(ie,2), 3); Z1 = Nodes(Elements(ie,2), 4); X2 = Nodes(Elements(ie,3), 2); Y2 = Nodes(Elements(ie,3), 3); Z2 = Nodes(Elements(ie,3), 4);
le(ie,1) = sqrt((X2-X1)^2+(Y2-Y1)^2+(Z2-Z1)^2);    %Holds length of element
ns(ie,1)=(Z2-Z1)/le(ie,1); ms(ie,1)=(Y2-Y1)/le(ie,1);  ls(ie,1)=(X2-X1)/le(ie,1);
E = Elements(ie,4);                          %Holds modolus of elasticiy of element [N/m2]
A = Elements(ie,5);                          %Holds cross sectional area of element
Te=[ls(ie,1)*ls(ie,1) ms(ie,1)*ls(ie,1) ns(ie,1)*ls(ie,1);...
    ls(ie,1)*ms(ie,1) ms(ie,1)*ms(ie,1) ns(ie,1)*ms(ie,1);...
    ls(ie,1)*ns(ie,1) ns(ie,1)*ms(ie,1) ns(ie,1)*ns(ie,1)];
K(DOFs,DOFs) = K(DOFs,DOFs) + (E*A/le(ie,1))*[Te -Te;-Te Te]; % Calculates the element stiffness matrix and assembles it to the global stiffness matrix
end
Ki=K;
figure(1)
for ie = 1:numelem                           %without stress con.
DOFs = [3*Elements(ie, 2)-2, 3*Elements(ie, 2)-1, 3*Elements(ie, 2), 3*Elements(ie, 3)-2, 3*Elements(ie, 3)-1, 3*Elements(ie, 3) ]; %Holds element’s DOFs for an elem
x=[Nodes(Elements(ie,2), 2) Nodes(Elements(ie,3), 2)]; y=[Nodes(Elements(ie,2), 3) Nodes(Elements(ie,3), 3)]; z=[Nodes(Elements(ie,2), 4) Nodes(Elements(ie,3), 4)];
plot3(x,y,z,'k') 
hold on
grid on
end
coef=100;                                                                  %coefficient of force
[numforce,nvm]=size(Forces);
for ie=1:numforce
DOFs = [3*Forces(ie, 1)-2 ; 3*Forces(ie, 1)-1 ; 3*Forces(ie, 1)];          %Holds element’s DOFs for an elem
F(DOFs)= Forces(ie,2:4)';
end
maxd=(max(diag(K)))*10^4;                                                  %diag enhancement - penalty method
for ie=1:length(constr)
   K(constr(1,ie),constr(1,ie))=K(constr(1,ie),constr(1,ie))+maxd; 
end
%U(act,1) = Calc_inv(Ki(act,act))*F(act,1);                                %elimination method
%U = Calc_inv(K)*F;                                                        %penalty method
F=coef*F;
%if Calc_det(K)~=0 don't utilize this condi. in order 2 save time
if det(K)~=0                                                               %under this cond. => solution  
U = Calc_inv(K)*F;                                                         %[m] - penalty method
Us=reshape(U,3,length(Nodes));  Nodes(:,2)=Nodes(:,2)+Us(1,:)'; Nodes(:,3)=Nodes(:,3)+Us(2,:)'; Nodes(:,4)=Nodes(:,4)+Us(3,:)'; 
f=Ki*U;                                                                    %force spplied on ndoes
for ie = 1:numelem                                                         %stress con.
DOFs = [3*Elements(ie, 2)-2, 3*Elements(ie, 2)-1, 3*Elements(ie, 2), 3*Elements(ie, 3)-2, 3*Elements(ie, 3)-1, 3*Elements(ie, 3) ]; %Holds element’s DOFs for an elem
x=[Nodes(Elements(ie,2), 2) Nodes(Elements(ie,3), 2)]; y=[Nodes(Elements(ie,2), 3) Nodes(Elements(ie,3), 3)]; z=[Nodes(Elements(ie,2), 4) Nodes(Elements(ie,3), 4)];
plot3(x,y,z,'b') 
hold on
if Elements(ie,3)==21
 text(x(1,1),y(1,1),z(1,1),['\bullet \leftarrow ',num2str(Elements(ie,2))]); text(x(1,2),y(1,2),z(1,2),['\bullet \leftarrow B',num2str(Elements(ie,3))]);
elseif  Elements(ie,3)==22
 text(x(1,1),y(1,1),z(1,1),['\bullet \leftarrow ',num2str(Elements(ie,2))]); text(x(1,2),y(1,2),z(1,2),['\bullet \leftarrow A',num2str(Elements(ie,3))]);
else Elements(ie,3)~=22 && Elements(ie,3)~=21;
 text(x(1,1),y(1,1),z(1,1),['\bullet \leftarrow ',num2str(Elements(ie,2))]); text(x(1,2),y(1,2),z(1,2),['\bullet \leftarrow ',num2str(Elements(ie,3))]);
end
grid on
end
title(['coef-of-forces:',num2str(coef)]); xlabel('x[m]');ylabel('y[m]');zlabel('z[m]'); forcendoes=zeros(numelem,6); Fornodes=zeros(numnode,3);
for ie = 1:numelem
    % calc. / element => U [1,6] vector
DOFs = [3*Elements(ie, 2)-2, 3*Elements(ie, 2)-1, 3*Elements(ie, 2), 3*Elements(ie, 3)-2, 3*Elements(ie, 3)-1, 3*Elements(ie, 3) ]; %Holds element’s DOFs for an elem
X1=Nodes(Elements(ie,2),2); Y1=Nodes(Elements(ie,2), 3); Z1=Nodes(Elements(ie,2), 4); X2=Nodes(Elements(ie,3),2); Y2=Nodes(Elements(ie,3),3); Z2=Nodes(Elements(ie,3),4);
T.U = [ls(ie,1) ms(ie,1) ns(ie,1) 0 0 0; 0 0 0 ls(ie,1) ms(ie,1) ns(ie,1)]*U(DOFs);    %preallocated - determined from non-loaded condition - initial condition
B=[-1/le(ie,1) 1/le(ie,1)];
strain(ie,1) = B*T.U;                                                        %calc. / element
stress(ie,1)= Elements(ie, 4)*strain(ie,1);
axialforce(ie,1) = stress(ie,1)*Elements(ie,5);
forcendoes(ie,:) = [-(([ls(ie,1) ms(ie,1) ns(ie,1)'])*axialforce(ie,1)),(([ls(ie,1) ms(ie,1) ns(ie,1)'])*axialforce(ie,1))]; %forces in every node
end
for ie=1:numelem
    Fornodes(Elements(ie,2),[1 2 3])=forcendoes(ie,[1 2 3]);
    Fornodes(Elements(ie,3),[1 2 3])=forcendoes(ie,[4 5 6]);
end
fr=reshape(f,3,numnode);
Fornodes=Fornodes';    %Failure in this method 
SumF=zeros(6,1);
for i=1:numnode
    SumF=SumF+[Fornodes(:,i) ; fr(:,i)];
end
stress=stress*(10^-6);                                                     %conv. from [N/m^2] to [N/mm^2]*(10^-6)
maxstress=max(stress); fconstrnodes=f(constr);
nod=fopen('AEM5380_Results_truss.txt','w+t');
fprintf(nod,  '#NODE    ,Xf-Xi=UX[m] ,Yf-Yi=UY[m] ,Zf-Zi=UZ[m] ,      Fx[N]    ,    Fy[N]    ,    Fz[N]  \n');
fprintf('\n');
for i=1:numnode
fprintf(nod,'%6.0f   %f    %f    %f       %6.2f      %6.2f         %6.2f \n',[i, Us(:,i)'  ,  fr(1,i)  ,  fr(2,i) ,   fr(3,i)] );
fprintf('\n');
end
fprintf('\n');
fprintf(nod,'     opou \n');
fprintf('\n');
fprintf(nod,'  B=21    %f     %f     %f     %6.2f      %6.2f         %6.2f  \n',[ Us(:,21)'  ,  fr(1,21)  ,  fr(2,21) ,   fr(3,21)] );
fprintf('\n');
fprintf(nod,'  A=22    %f     %f     %f     %6.2f      %6.2f         %6.2f  \n',[ Us(:,22)'  ,  fr(1,22)  ,  fr(2,22) ,   fr(3,22)] );
fprintf('\n');
fprintf(nod,  '#Elements    strains             stress [N/mm^2]       axialforce[N]/elem   \n');
for i=1:length(strain)
fprintf(nod, '%6.0f          %f        %6.2f                %6.2f   \n', [ i ,  strain(i,1)  ,  stress(i,1),   axialforce(i,1) ]);
fprintf('\n');
end
fprintf(nod,  '#bottom(reactions)        Fx[N]               Fy[N]           Fz[N]  \n');
fprintf('\n');
for i=1:4
fprintf(nod, '%f                 %6.2f          %6.2f        %6.2f \n', [i  , -fr(1,i)   ,   -fr(2,i)  ,  -fr(3,i) ]);
fprintf('\n');
end
fclose(nod); clc;
else
    fprintf('there is no solution , due to singularity \n')
end