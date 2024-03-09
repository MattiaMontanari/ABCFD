% EXEMPLE FROM ADRIEN
% XY=Pts;
% T20=FACES;
% load('C:\Users\XMAMON\Downloads\mesh.mat')
function [K , M ,B ] = exemple_2D_TRI( XY , T20 )

% X = XY(:,1);
% Y = XY(:,2);
% 
% segments = [1 2; 2 3 ; 3 1];
% 
% [T21,T10,T12] = FEM_Process_T20(T20,segments);

%TXY = table defining entities of dimension X in terms of entities of
%dimension Y
%example: T20 :  2D entities (triangles)
% 
% 
% figure;
% triplot(T20,X,Y);
Nelem = size(T20,1);
Nnodes = size(XY,1);

[Xi,Wi] = hammer_points(2);


% example compute a stiffness and  a mass matrix + a right hand side
f = 1; %(source term)

%initialization
K = sparse( Nnodes ,Nnodes);
M = sparse( Nnodes ,Nnodes);
B = zeros(Nnodes,1);
%loop over the elements
for el=1:Nelem
    %nodes numbers of the element
    map = T20(el,:);
    %nodes coordinates
    XY_el = XY(map,:);
    
    %function to evaluate the shape functions and their gradient at the
    %integration points
    [phi, phi_x,phi_y] = eval_N_TRI_2D(XY_el,Xi);
    %phi: value of the shape fcts size: 3xNintegrationPoints
    %phi_x: dphi/dx, constant over the element. size: 3x1
    %phi_y: dphi/dy, constant over the element. size: 3x1
    
    %function to compute the element jacobian (constant over the element)
    J = TRI_Jac(XY_el);
    
    %initialization of the local matrices
    K_loc = zeros(numel(map),numel(map));
    M_loc = zeros(numel(map),numel(map));
    B_loc = zeros(numel(map),1);

    %loop over the integration points
    for i=1:numel(Wi)      
        K_loc = K_loc + J*Wi(i)*  (phi_x*phi_x' + phi_y*phi_y');
        M_loc = M_loc + J*Wi(i)*  (phi(:,i)*phi(:,i)');
        
        B_loc = B_loc + J*Wi(i)*phi(:,i)*f;
        
    end
    if el == 1
         K(map,map) =  K_loc;
    M(map,map) = M_loc;
    
    B(map) = B_loc;
    else
    %assembly
    K(map,map) = K(map,map)+K_loc;
    M(map,map) = M(map,map)+M_loc;
    
    B(map) = B(map)+B_loc;
    end
end