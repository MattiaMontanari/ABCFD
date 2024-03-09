% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SOLVES LINEAR POD SYSTEM                       %
%  /----\ |  \|    |--  |   |                                                  %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %


function [ POD ]=abCFD_POD_solver( U,H,R,B,D_T,D_uv,K_T,K_uv,K_Tuv,ph,tSteps,ee,var,mesh)

%% I N I T I A L I Z E

% VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ta1 = mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ; % temperature
ta2 = mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ; % pressure (not used)
ta3 = mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG ; % x-velocity
ta4 = mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG ; % y-velocity
    % Tetha scheme
te = ph.tetha;
    % Number of time steps in current phase
nt = numel( tSteps ); 
    % Numer of ALL dofs
NumDof = sum( var.fieldNDofs );
NumT = double( var.fieldNDofs( 3 ) );
NumU = double(  var.fieldNDofs( 1)/2 );
    % Unknown array
Q = zeros( NumDof, nt );

% MODES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
BT = B( ta1 , :);
BU = B( ta3 , :);
BV = B( ta4 , :);

% MATRICES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    %Constant
D_T     = D_T( ta1, ta1);
D_uv    = D_uv( [ta3, ta4] , [ta3, ta4] );
K_T     = K_T( ta1, ta1);
% K_uv    = K_uv( [ta3, ta4] , [ta3, ta4] ); 
K_uv = [ K_uv(ta3,ta3)+K_uv(ta4,ta4) , zeros(NumU ,NumU );
            zeros(NumU , NumU)       ,  K_uv(ta3,ta3)+K_uv(ta4,ta4) ];
V1      = K_Tuv( ta1 , ta3 );
V2      = K_Tuv( ta1 , ta4 );
% Cx      = Cx(ta3,ta3);
% Cy      = Cy(ta3,ta3);
%% ASSEMBLE REDUCED OREDER MATRICES

% Damping ( constant in time )
D(   1    :   ee   ,    1   :   ee   ) =       BT'    *  D_T *     BT;
D( ee + 1 : 2 * ee , ee + 1 : 2 * ee ) = [ BU ; BV ]' * D_uv * [ BU ; BV];

% Stiffness 
K = zeros( 2*ee , 2*ee);
K(   1    :   ee   ,    1   :   ee   ) =       BT'    *  K_T *     BT;
K( ee + 1 : 2 * ee , ee + 1 : 2 * ee ) = [ BU ; BV ]' * K_uv * [ BU ; BV];


%% INITIALIZE EMPIRACAL COEFFICIENTS
omg = zeros( ee , nt);
alp = zeros( ee , nt);
% a(:,1) = BDB \ (B' * D * H.d1(:,2) ) ;
omg(:,1) = D(1:ee,1:ee) \ (BT' * D_T * H( ta1 ,2) ) ;
alp(:,1) = D( ee + 1 : 2 * ee , ee + 1 : 2 * ee ) \  ([BU ; BV]' * D_uv * H( [ta3;ta4] ,2) ) ;
% alp(:,1) = BDB \ (B' * D * H.d1(:,2) ) ;

%% COMPUTE SOLUTION IN TIME

for id_t = 2 : 1 : nt + 1 - 1
    % Current time step
    dt = tSteps( id_t , 1 ) - tSteps( id_t-1 , 1 );
    
    fprintf('Coeff. %d, %d, %d, %d\n', alp(1,id_t-1),alp(2,id_t-1),...
                    omg(1,id_t-1) ,omg(2,id_t-1) )
    
    
    % ENERGY CONVECTION  
    
    V_new = BT' *( [spdiags(   R( ta1,  id_t  )   , 0 ,NumT,NumT ) * V1 , ...
                    spdiags(   R( ta1,  id_t  )   , 0 ,NumT,NumT ) * V2 ] + ...
                   [spdiags( BT * omg(:,id_t - 1) , 0 ,NumT,NumT ) * V1 , ...
                    spdiags( BT * omg(:,id_t - 1) , 0 ,NumT,NumT ) * V2 ] ) ; 
                
    V_old = BT' *( [spdiags(   R( ta1,  id_t - 1 ), 0 ,NumT,NumT ) * V1 , ...
                    spdiags(   R( ta1,  id_t - 1 ), 0 ,NumT,NumT ) * V2 ] + ...
                   [spdiags( BT * omg(:,id_t - 1) , 0 ,NumT,NumT ) * V1 , ...
                    spdiags( BT * omg(:,id_t - 1) , 0 ,NumT,NumT ) * V2 ] ) ; 
    % MOMENTUM CONVECTION - unkwoen field
    Crom_old = abCFD_AssembleNavier( var, [ mesh{2}.ele{1}.dof{3}.TAG;mesh{2}.ele{1}.dof{3}.TAG],...
                Q( :, id_t - 1 ) , mesh);
    Crom_old = Crom_old( ta3, ta3 ) ;
    % MOMENTUM CONVECTION - reference field
    Cref_new = abCFD_AssembleNavier( var, [ mesh{2}.ele{1}.dof{3}.TAG;mesh{2}.ele{1}.dof{3}.TAG],...
                R( :, id_t  ) , mesh);
    Cref_new = Cref_new( ta3, ta3 );
    
    Cref_old = abCFD_AssembleNavier( var, [ mesh{2}.ele{1}.dof{3}.TAG;mesh{2}.ele{1}.dof{3}.TAG],...
                R( :, id_t - 1  ) , mesh);
	Cref_old = Cref_old( ta3, ta3 ) ;
 
    % ASSEMBLE LINEARIZED SYSTEM
    % Momentum - reference field contribute to R.H.S.
    Mrhs =     te    * ( BU' * (Cref_new+Crom_old) * R( ta3,  id_t  ) + ...
                  BV' * ( Cref_new + Crom_old ) * R( ta4,  id_t  ))  + ...
            ( 1- te )* ( BU' * (Cref_old+Crom_old) * R( ta3,  id_t - 1   ) + ...
                  BV' * ( Cref_old + Crom_old ) * R( ta4,  id_t - 1  ))  ;
    % Energy - reference field contribute to R.H.S.
    Erhs =   te    * V_new * [ R( ta3,    id_t   ) ; R( ta4,    id_t   ) ]  + ...
         ( 1 - te) * V_old * [ R( ta3,  id_t - 1 ) ; R( ta4,  id_t - 1 ) ] ;
    
    LHS = D./dt +   te .*  K;
    LHS( 1 : ee , ee + 1 : 2 * ee) = te * V_new * [ BU ; BV];
    LHS( ee + 1 : 2 * ee , ee + 1 : 2 * ee) = LHS( ee + 1 : 2 * ee , ee + 1 : 2 * ee) + ...
        te * ( BU' * (Cref_new+Crom_old) * BU  + ...
               BV' * (Cref_new+Crom_old) * BV);
           
    RHS = D./dt - (1 - te).* K;
    RHS( 1 : ee , ee + 1 : 2 * ee) = - (1-te) * V_old * [ BU ; BV];
    RHS( ee + 1 : 2 * ee , ee + 1 : 2 * ee) = RHS( ee + 1 : 2 * ee , ee + 1 : 2 * ee) + ...
        (-1)*(1-te) * ( BU' * (Cref_old+Crom_old) * BU  + ...
               BV' * (Cref_old+Crom_old) * BV);
	RHS = RHS * [ alp( : , id_t-1) ; omg( : , id_t-1)] + [ Erhs ; Mrhs];
    
    tmp = LHS \ RHS;
    alp( : , id_t ) = tmp(    1   :   ee   , 1);
    omg( : , id_t ) = tmp( 1 + ee : 2 * ee , 1);
    % Recontrstuct physical solution
    Q( ta1, id_t) = BT * alp( : , id_t );
    Q( [ ta3 ; ta4], id_t ) = [BU ; BV] * omg( : , id_t );
end

POD( [ ta1;ta3;ta4],:)= R( [ ta1;ta3;ta4],:) + Q( [ ta1;ta3;ta4],:);
% POD(  [ ta1;ta3;ta4],:) , 1 ) = U( ( [ ta1;ta3;ta4],:) , 1 );

end

