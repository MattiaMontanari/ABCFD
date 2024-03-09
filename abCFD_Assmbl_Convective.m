% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	ASSEMBLE CONVECTIVE MATRIX FOR A VECTOR AND    %
%  /----\ |  \|    |--  |   |   A SCALAR FIELD GIVEN A SOLUTION ARRAY.         %
% /      \|__/ \__ |    |__/    RUNS FOR LINEAR TRIANGULAR ELEMENTS ONLY, THE  %
%                               CONNECTIVITY MATRICES ARE INCLUDED INTO THE    %
% INPUT STRCUTRES. THIS IS AN OPTIMIZED SOLUTION WHITOUT FOR LOOPS.            %                              
%  * * * calls * * *                                                           %
%               i.  abCFD_GaussQuad2D                                          %
%              ii.  abCFD_round                                                %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ C , K, M] = abCFD_Assmbl_Convective( SOL, mesh, var )
% INPUTS:
%   - SOL   : Solution vector to be included into the matrix
%   - mesh  : Usual structure
%   - var   : Usual structure
% OUTPUT:
%   - C     : Convective matrix for Momentum Equations
%   - V     : Convective matrix for Energy Equation

function [ C , V ] =  abCFD_Assmbl_Convective( SOL, mesh, var, check)

% TEST
% [K_1del1,K_del12,K_2del1,K_22del , Mdel ,Bdel,CTEST, VTEST] = exemple_2D_TIA(var.dofCoord' ,SOL, mesh );
%% INITIALIZE 

% Energy equation variable's tag
V_VarTag = 1;
% Energy equation mapping tags
V_MapTag = [ 3, 4 ];
% Momentum equation variable's tag
C_VarTag = [ 3, 4 ];
% Momentum equation mapping tags
C_MapTag = [ 3, 4 ];
% Integration order
IntPrecision = 2;

%% EXTRACT DATA FROM INPUT STRUCTURES
% In the following the indeces _0,_1,_2 stand for: scalar field, first and
% second variables of vector field, respectively.

% Tot. number of dofs i.e. size of whole output matrices
TotDofs = size( var.dofCoord , 2 );
% Connectivity matrix for Energy Equation leading variable (scalar field)
    V_Conn_T0 = mesh{2}.ele{1}.dof{ V_VarTag }.TAG;
% Connectivity matrix for column mappig for Energy equations
    V_Conn_R1 = mesh{2}.ele{1}.dof{ V_MapTag( 1 ) }.TAG;
    V_Conn_R2 = mesh{2}.ele{1}.dof{ V_MapTag( 2 ) }.TAG;
% Solution scalar field to be included into Energy equation
V_Sol = SOL( V_Conn_T0  ); 
% Connectivity matrices for Momentum Equation leading variables (vector field)
C_Conn_1 = mesh{ 2 }.ele{ 1 }.dof{ C_VarTag( 1 ) }.TAG;
C_Conn_2 = mesh{ 2 }.ele{ 1 }.dof{ C_VarTag( 2 ) }.TAG;
% Solution scalar field to be included into Momentum equation
C_Sol_1 = SOL( C_Conn_1  ); 
C_Sol_2 = SOL( C_Conn_2  ); 

   
%% FINITE ELEMENTS DATA ON LOCAL COORDINATES

% Number of FEM elements and number of Shape functions per element
[NumEle , NumShp] = size( C_Conn_1  );
% Gauss Points coords and weights
[ sGP , tGP, WGP ] = abCFD_GaussQuad2D( mesh{ 2 }.Type , IntPrecision );

% Evaluates the shape functions and their derivatives on an element for some
%   points specified in the parent element
XI = [ sGP ; tGP ]; 
PHI = [1-XI(1,:)-XI(2,:)    ; ...
         XI(1,:)            ; ...
                 XI(2,:)    ]  ;
DPHIDXI = [ -1 -1 ; 1 0 ; 0 1]; % [ NumShp x 2 ]
% Numerical Integral of shape functions' product
IntGtp = PHI * diag(WGP) * PHI';

%% COMMON DATA FOR ENERGY AND MOMENTUM EQUATIONS
% These data should be equal for both Energy and Momentum question since
% the are discretized on the same grid

% Energy leading variable mapping [ 1 x NumShp x NumEle ]. This is used in
%   the next line to obtain the nodes coordinates. So, would make no
%   difference to use the energy leading variable of one of the momentum equation
C_MapVar_1 = reshape( C_Conn_1', 1, 3, []);
% Elements global coordinates [ NumShp x 2 x NumEle ]
COORD = permute(reshape( var.dofCoord( : , permute( C_MapVar_1, [2 1 3] )) , 2, 3,[] ) , [ 2 1 3]);
% Global f.o.r derivative over Local f.o.r. derivatives
DXDXI = reshape( shiftdim( COORD, 1) , [] , 3 ) * DPHIDXI;
% Compute the jacobians
rr( : , 1 ) = + DXDXI(1:2:end-1, 1);
rr( : , 2 ) = - DXDXI(1:2:end-1, 2);
cc = DXDXI(2:2:end  , :);
cc = cc( : , [2 1] );
jacobians = sum(rr.*cc,2); 
% Shape functions' paritial derivatives w.r.t. global f.o.r.
Oodd =   DXDXI(    1      : 2 :  end/2  ) ./ jacobians' ;
Opar =   DXDXI( end/2 + 2 : 2 :   end   ) ./ jacobians' ;
Xodd = - DXDXI( end/2 + 1 : 2 :   end   ) ./ jacobians' ;
Xpar = - DXDXI(    2      : 2 :  end/2  ) ./ jacobians' ;
OX = zeros( size( DXDXI' ) );
OX( 2 , 2 : 2 : end ) = Oodd;
OX( 1 , 1 : 2 : end ) = Opar;
OX( 1 , 2 : 2 : end ) = Xodd;
OX( 2 , 1 : 2 : end ) = Xpar;
DPHIDX = DPHIDXI*OX ;
DPHI_X = DPHIDX( : , 1 : 2 : end );
DPHI_Y = DPHIDX( : , 2 : 2 : end );
% 
shpdopsWgt = PHI*diag( sqrt(WGP) );
order = [ shpdopsWgt ; repmat( zeros(NumShp,1) , NumShp,NumShp )];
% Shape functions product contributions separated by Gauss weithgs
shpdotshp = reshape(order, NumShp*NumShp,[])*reshape(order, NumShp*NumShp,[])';
% Sum of three contributes due to gauss points
sumintopar = sum( reshape( sum( reshape( shpdotshp, [],NumShp, NumShp ),3)' , NumShp,NumShp,NumShp),3);

%% ENERGY EQUATION CONTRIBUTE

% Energy mapping on columns
V_MapRow_T0 = reshape( V_Conn_T0', 1, 3, []);
V_MapCol_R1 = reshape( V_Conn_R1', 1, 3, []);
V_MapCol_R2 = reshape( V_Conn_R2', 1, 3, []);

% Product of shape funtions derivatives time solution  field (constant for each element)
V_N_xSol = sum( DPHI_X' .* V_Sol , 2 ) .* ( jacobians ) ;
V_N_ySol = sum( DPHI_Y' .* V_Sol , 2 ) .* ( jacobians ) ;
% Convective matrices contributes per each element
V_Loc_x = IntGtp(:)*( V_N_xSol' );
V_Loc_y = IntGtp(:)*( V_N_ySol' );
% Assemble the Energy convective matrix
V_iI_x = reshape( repmat( V_MapRow_T0 , [1 ,3 ,1 ]) , [] , 1 );
V_jJ_x = reshape( repmat( V_MapCol_R1( : )' , NumShp, 1 ) , [] , 1 );

V_iI_y = reshape( repmat( V_MapRow_T0 , [1 ,3 ,1 ]) , [] , 1 );
V_jJ_y = reshape( repmat( V_MapCol_R2( : )' , NumShp, 1 ) , [] , 1 );
 
Vx = sparse( V_iI_x(:) , V_jJ_x(:) , V_Loc_x( : ) , TotDofs, TotDofs);

V = sparse( [ V_iI_x(:) ; V_iI_y(:) ], [V_jJ_x(:) ; V_jJ_y(:)] , ...
                        [ V_Loc_x( : ) ; V_Loc_y( : ) ] , TotDofs, TotDofs);
 
%% MOMENTUM EQUATION CONTRIBUTE

% Momentum leading variables mapping [ 1 x NumShp x NumEle ].
C_MapVar_1 = reshape( C_Conn_1', 1, 3, []);
C_MapVar_2 = reshape( C_Conn_2', 1, 3, []);

% shpdopsWgt = PHI*diag( sqrt(WGP) );
% 
% order = [ shpdopsWgt ; repmat( zeros(NumShp,1) , NumShp,NumShp )];
% 
% shpdotshp = reshape(order, NumShp*NumShp,[])*reshape(order, NumShp*NumShp,[])';
% 
% sumintopar = sum( reshape( sum( reshape( shpdotshp, [],NumShp, NumShp ),3)' , NumShp,NumShp,NumShp),3);

% Solution X-component velocity
SolX = reshape( C_Sol_1', NumShp,1, NumEle );
% Global derivative X-component
dpX = reshape( DPHI_X , 1 , NumShp, NumEle);
% Multiply solution X-component velocity with the shape functs' derivative in X 
ok_X = bsxfun( @times, SolX,dpX);

% Solution Y-component velocity
SolY = reshape( C_Sol_2', NumShp,1, NumEle);
% Global derivative Y-component
dpY = reshape( DPHI_Y , 1 , NumShp, NumEle);
% Multiply solution Y-component velocity with the shape functs' derivative in Y
ok_y = bsxfun( @times, SolY,dpY);

% sum contributes X and Y
ok_XY =  ok_X + ok_y;
% Multiply by shape functions
IntSum = sumintopar*reshape( ok_XY, NumShp, [] );
% Local values of the convective matrix per element
% OLD % C_locals =  reshape( IntSum , NumShp*NumShp, [] ) *diag( jacobians);
C_locals= bsxfun( @times, reshape( IntSum , NumShp*NumShp, [] ), jacobians');

C_iI_1 = reshape( repmat( C_MapVar_1 , [1 ,3 ,1 ]) , [] , 1 );
C_jJ_1 = reshape( repmat( C_MapVar_1( : )' , NumShp, 1 ) , [] , 1 );

C_iI_2 = reshape( repmat( C_MapVar_2 , [1 ,3 ,1 ]) , [] , 1 );
C_jJ_2 = reshape( repmat( C_MapVar_2( : )' , NumShp, 1 ) , [] , 1 );
 
C = sparse( [ C_iI_1(:) ; C_iI_2(:) ], [C_jJ_1(:) ; C_jJ_2(:)] , ...
                        [ C_locals( : ) ; C_locals( : ) ] , TotDofs, TotDofs);

%% CHECK SOLUTION INCOMPRESSIBILITY

if strcmp( check,'CHECK')
  
  Incompr = mean( sum( DPHI_X'.*SOL( mesh{2}.ele{1}.dof{3}.TAG ) , 2) + ...
                        sum( DPHI_Y'.*SOL( mesh{2}.ele{1}.dof{4}.TAG ) , 2) );
divergence = sum( DPHI_X'.*SOL( mesh{2}.ele{1}.dof{3}.TAG ) , 2) + ...
                        sum( DPHI_Y'.*SOL( mesh{2}.ele{1}.dof{4}.TAG ) , 2) ;
%   Incompr2 = max(abs(( divergence )));
  
	[value, location] = max(divergence(:));

  fprintf('Worst case of incompressibility %d at node %d \n',value, location)
                    
    if abs( Incompr ) > 10^-5
    
        warning('Solution used is not incompressible for a tollerance of 10^-5')
    else
        disp('Solution Incompressibility verified')
    end
    
end


return

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:   May  2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  09/05/2013 %
% ---------------------------------------------------------------------------- %