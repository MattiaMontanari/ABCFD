function test_ReducedOrderModel( var, mesh, REF, SOL, Ee, model, ph,varargin)
 
TRI =  mesh{2}.ele{1}.dof{1}.TAG;

%%   =======================================================================  %%
%% - - - - - - -   C O M P U T E   C O E F F.   M A T R I C E S  - - - - - - - %
%%   =======================================================================  %%

% In this part compute space-discretization matrices of large size. 
% These matrices are used to approximate the modes with linear triangles 
%   elements. 

% Total number of degrees of freedom
TotDofs = sum(var.fieldNDofs);
% Nodes' tags per each dependent variable
ta1 = mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ; % temperature
ta2 = mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ; % pressure (not used)
ta3 = mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG ; % x-velocity
ta4 = mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG ; % y-velocity

% V E R I F I C A T I O N    R O U T I N E S
% [K_11,K_12,K_21,K_22 , Ddel ,Ptmp ] = exemple_2D_TRI( var.dofCoord' , mesh{2}.ele{1}.dof{1}.TAG );
% [ CrefRuv , VrefRT ] =  abCFD_Assmbl_Convective( REF.d1, mesh, var, 'CHECK' );


% Assemble matrices Dpp = Dnn = Dnp == D
[ K11, K12, K21, K22, D ] = abCFD_AssembleNavier( var, TRI );
D = abCFD_round( D ) ;
Dpp = sparse( TotDofs , TotDofs);
Dpp( ta1 , ta1 ) = D( ta1, ta1 );
Dnn = sparse( TotDofs , TotDofs);
Dnn( ta3 , ta3 ) = D( ta1, ta1 );
Dnn( ta4 , ta4 ) = D( ta1, ta1 );
Dnp = sparse( TotDofs , TotDofs);
Dnp( ta3 , ta1 ) = D( ta1, ta1 );
Dnp( ta4 , ta1 ) = D( ta1, ta1 );
% Assembl matrices Kii and Kij
Kpp_jj = sparse( TotDofs , TotDofs);
Kpp_jj( ta1 , ta1 ) = K11( ta1, ta1 ) + K22( ta1, ta1 );
Knn_jj = sparse( TotDofs , TotDofs);
Knn_jj( ta3 , ta3 ) = K11( ta1, ta1 ) + K22( ta1, ta1 );
Knn_jj( ta4 , ta4 ) = K11( ta1, ta1 ) + K22( ta1, ta1 );
Knn_ij = sparse( TotDofs , TotDofs);
Knn_ij( ta3 , ta3 ) = K11( ta1, ta1 ) ;
Knn_ij( ta3 , ta4 ) = K12( ta1, ta1 ) ;
Knn_ij( ta4 , ta3 ) = K21( ta1, ta1 ) ;
Knn_ij( ta4 , ta4 ) = K22( ta1, ta1 ) ;

% Compute convective matrices
[ Crefnnn , VrefTppn ] =  abCFD_Assmbl_Convective( REF.d1 , mesh, var, 'CHECK' );

% Body force terms
% [Kdel , Ddel ,Ptmp ] = exemple_2D_TRI( var.dofCoord' , mesh{2}.ele{1}.dof{1}.TAG );clear Kdel; clear Ddel
Pn = sparse( TotDofs , 1 );
% Pn( ta3 , 1 ) = Ptmp( ta1 , 1 );
% Pn( ta4 , 1 ) = Ptmp( ta1 , 1 );

%% Spy coefficient matrices with elements sorted as [ T,p,U,V]
figure;
subplot( 1 , 3 , 1)
spy( Dpp([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'y')
hold on
spy( Dnn([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'r')
spy( Dnp([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'g')
legend('D_p_p','D_n_n','D_n_p','Location','Best')
title('Damping matrix')

subplot( 1 , 3 , 2)
spy( Kpp_jj([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'m')
hold on
spy( Knn_ij([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'co')
spy( Knn_jj([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'b')
legend('Kpp_j_j','Knn_i_j','Knn_j_j','Location','Best')
title('Stiffness matrix')

subplot( 1 , 3 , 3)
spy( Crefnnn([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'m')
hold on
spy( VrefTppn([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]) ,'c')
legend('Cref_n_n_n','VrefT_p_p_n','Location','Best')
title('Convection matrix')

%%   =======================================================================  %%
%% - - - - - - - - -   S N A P S H O T S    C A P T U R I N G    - - - - - - - %
%%   =======================================================================  %%

% initialize from input
% userLimits = varargin{ 1 };
% userSteps  = varargin{ 2 };

%% EXTRACT DATA FROM FEM COMPUTATION
% All steps taken by the solver
FemTimeSteps = SOL.solinfo.solvals;
% Nodes' tags at outlet
TaOutUx = mesh{1}.ele{3}.dof{ 3 }.uniqTAG ;
TaOutUy = mesh{1}.ele{3}.dof{ 4 }.uniqTAG ;
TaOutT  = mesh{1}.ele{3}.dof{ 1 }.uniqTAG ;
% Define VELOCITY Probes points
samplesUx= TaOutUy( [ round( numel(TaOutUx) /6) ,  round( numel(TaOutUx) /2) ,  round( numel(TaOutUx) *3/4) ]);
samplesUy= TaOutUy( [ round( numel(TaOutUy) /6) ,  round( numel(TaOutUy) /2) ,  round( numel(TaOutUy) *3/4) ]);
samplesT = TaOutT(  [ round( numel(TaOutT) /6) ,   round( numel(TaOutT) /2) ,   round( numel(TaOutT) *3/4) ]);

% Plot the Velocity profile @Outlet in time for three points
figure();
% Plot sampling points
subplot( 2 , 3, [3 6] )
title('Model mesh with edges'' tags and sampling points')
mphmesh(model,'meshFEM','Facealpha',.1,'edgecolor', 'w','Edgelabels','on', 'Edgelabelscolor','r')
hold all
plot( var.dofCoord( 1 , samplesUy ) , var.dofCoord( 2 , samplesUy  ) ,'x','MarkerSize',15)
axis tight
legend( 'Domain', 'Edges', 'Probes' ,'Location','Best')

% Nodes' tag for three probes elements
subplot( 2 , 3, [ 1 2 4 5 ] )   
hold all
plot( FemTimeSteps , max(max(abs(sqrt( SOL.d1( samplesUx , : ).^2 + SOL.d1( samplesUy , : ).^2) )))./sqrt( SOL.d1( samplesUx , : ).^2 + SOL.d1( samplesUy , : ).^2) );
plot( FemTimeSteps , max(max(abs(SOL.d1( samplesT , : ) )))./SOL.d1( samplesT , : ) , ':');
axis tight
title('Normalized Velocity and Temperature at outlet')
xlabel('Time [s]')
grid on
legend( 'U magnitude - Left Probe', 'U magnitude - Center Probe','U magnitude - Right Probe',...
        'Temperature - Left Probe', 'Temperature - Center Probe','Temperature - Right Probe','Location','Best')
% Insert data from plot
disp('Click on snapshots starting capturing point')
% Initial and final time for capturing snapshotss
if numel(varargin) == 0
    [InitialTime, y] = ginput(1);
    pause(.5)
    disp('Click on snapshots ending capturing point')
    % Final time for capturing snapthos
    [FinalTime, y] = ginput(1);
else
    InitialTime = varargin{1}(1);
    FinalTime =   varargin{1}(2);
end


IN_QuasiSteady = find( abCFD_round( FemTimeSteps ) > abCFD_round( InitialTime ) );
EN_SnapCaputur = find( abCFD_round( FemTimeSteps ) > abCFD_round( FinalTime ) );
 % keep important data only
IN_QuasiSteady = IN_QuasiSteady(1);
EN_SnapCaputur = EN_SnapCaputur(1);
% Print to screen
fprintf('Quasi-steady regime detected at t( %d ) = %d \n',IN_QuasiSteady(1),FemTimeSteps(IN_QuasiSteady(1)) )

%% USER DEFINES HOW TO TAKE SNAPSHOTS

SnapMethod = 'Double';

fprintf('Snapshots taken with method %s from t=%d[s] to %d[s] \n',SnapMethod,FemTimeSteps(IN_QuasiSteady(1)),FemTimeSteps(EN_SnapCaputur(1)) )
% Compute time step between snapshots
TimeBetween = sum( [ SOL.solinfo.solvals( 1 : end-1) , SOL.solinfo.solvals(2:end)]*[1 ; -1] , 2);
subplot( 2 , 3, [ 1 2 4 5 ] )
plot( SOL.solinfo.solvals( 2 : end), abs(TimeBetween),'xr','MarkerSize',10)
title( 'Snapshots History' )
ylabel( 'Time step size [s]')
xlabel( 'Physical time [s]')
legend( 'FEM Time steps taken by the solver','Location','Best')
grid on
hold on
AX = axis;
SquareSnapCaputure = [ FemTimeSteps(IN_QuasiSteady(1)), AX(3) ; ...
                       FemTimeSteps(EN_SnapCaputur(1)), AX(3) ; ...
                       FemTimeSteps(EN_SnapCaputur(1)), AX(4) ; ...
                       FemTimeSteps(IN_QuasiSteady(1)), AX(4) ; ...
                       FemTimeSteps(IN_QuasiSteady(1)), AX(3) ];
plot( SquareSnapCaputure(:,1),SquareSnapCaputure(:,2), 'c-', 'LineWidth', 4)
legend( 'FEM Time steps taken by the solver', 'Snapshots Capute interval','Location','Best')


switch SnapMethod
        case 'Double' % Define two intervals to take snapshots with different time steps
            
        disp('Click on snapshots capturing MIDDLE point')
        if numel(varargin)==0
            [MidTime, y] = ginput(1);
            SnapStep_1 = input( 'FIRST interval: consider every: 1st, 2nd, ... snapshot:  ');
            SnapStep_2 = input( 'SECOND interval: Consider every: 1st, 2nd, ... snapshot:  ');
        else
            MidTime = varargin{2}(1);
            SnapStep_1 = varargin{2}(2);
            SnapStep_2 = varargin{2}(3);
        end
        MID_QuasiSteady = find( abCFD_round( FemTimeSteps ) > abCFD_round( MidTime ) );
        SnapShots = (IN_QuasiSteady(1) : SnapStep_1 : MID_QuasiSteady(1))';
        SnapShots = [SnapShots ; (MID_QuasiSteady(1+SnapStep_2) : SnapStep_2 : EN_SnapCaputur(1))'];
        
        TimeBetwSnap = sum( [ FemTimeSteps(SnapShots( 1 : end-1)) , FemTimeSteps(SnapShots(2:end))]*[1 ; -1] , 2);
        semilogx( SOL.solinfo.solvals( SnapShots(2 :end) ), abs(TimeBetwSnap),'yo','MarkerSize',10)
        subplot( 2 , 3, [ 1 2 4 5 ] , 'replace')
        semilogy( SOL.solinfo.solvals( 2 : end), abs(TimeBetween),'xr','MarkerSize',10)
        hold on
        semilogy( SOL.solinfo.solvals( SnapShots(2 :end) ), abs(TimeBetwSnap),'yo','MarkerSize',10)
        grid on
        axis tight
        legend( 'Available solutions', 'Taken Snapshots','Location','Best')
        title( 'Snapshots History' )
        ylabel( 'Time step size [s]')
        xlabel( 'Physical time [s]')
        
end

% Plot FEM Kinetic energy during snapshots capturing and rest of simulation
FC_Kinetic = figure();
hold all
KinEn_FEM = 1/2*diag( SOL.d1(: , SnapShots )' * Dnn * SOL.d1(: , SnapShots ) );
% %         KinEn_POD = 1/2*diag( POD'   * Dnn *   POD   );
% % plot( tSteps, KinEn_POD,'-r','LineWidth',2 )
plot( SOL.solinfo.solvals( SnapShots  ), KinEn_FEM, 'o-.y','LineWidth', .5)
title('Kinetic energy captured by snapshots ','FontSize', 20,'interpreter','latex' );
legend( ['FEM Kinetic. Integral over time = ',num2str(sum(100*round(KinEn_FEM/100)))] ,'Location','Best')
xlabel('time [s]')
grid on

% Plot FEM Thermal energy during snapshots capturing and rest of simulation
FC_thermal = figure();
hold all
TherEn_FEM = diag( SOL.d1(: , SnapShots )' * Dpp * SOL.d1(: , SnapShots ) );

plot( SOL.solinfo.solvals( SnapShots  ), TherEn_FEM, 'o-.b','LineWidth', 2)
title('Thermal energy captured from first snapshots till final time','FontSize', 20,'interpreter','latex' );
legend( ['FEM Therm. Integral over time = ',num2str(sum(100*round(KinEn_FEM/100)))] ,'Location','Best')
xlabel('time [s]')
grid on


%%   =======================================================================  %%
%% - - - -   C O M P U T E     H O M O G E N E O U S     M O D E S   - - - - - %
%%   =======================================================================  %%

% The homogeneous solution isnt' stored but is comptued with:
%   bsxfun( @minus, SOL.d1(:,PodInvervals(1) : 2 : PodInvervals(2) )
% The modes are obtained using snapshots 

% Number of modes
M = numel( SnapShots );
tiin = toc;
ModesFun = 'eigUTU';
% Initialize modes
eeMin = Ee(1);
eeMax = Ee(2);
%% COMPUTE VELOCITY MODES 
ta34 = [ ta3 ; ta4 ];

Usnap = bsxfun( @minus,   SOL.d1([ta3;ta4],SnapShots ) , REF.d1([ta3;ta4],1));

fprintf('Start eigenvalue problem at %s with ''%s'' \n', datestr(now,15), ModesFun)
switch ModesFun
    case 'eigORT'
              
    H = 1/M .* ( (Usnap' * Dnn(ta34,ta34) * Usnap )+(Usnap' * Dnn(ta34,ta34) * Usnap )' );

    [ Vu , du ]= eig( ( H )) ;
    
    if isreal( Vu )
            du = (du( end:-1 :end-eeMax+1 , end:-1 :end-eeMax+1 ));
            Vu =  Vu(:,end:-1:end-eeMax+1);
    end

    case 'eigUTU'
    [ Vu , du ]= eig( 1/M .* ( Usnap' * Usnap ) );
    % For the function 'eig' is necessary to invert the outputs' order if
    % not immaginary
    if isreal( Vu )
            du = (du( end:-1 :end-eeMax+1 , end:-1 :end-eeMax+1 ));
            Vu =  Vu(:,end:-1:end-eeMax+1);
    end

    case 'eigsORT'

    [ Vu , du ]= eigs( (  1/M .* ( ( Usnap' * Dnn(ta34,ta34) * Usnap ) ) ), eeMax ) ;
    if ~isreal( Vu )
        warning( 'Not all velocity modes were real, Apply rounding to U'' * Dnn * U' )
        [ Vu , du ]= eigs( abCFD_round(  1/M .* ( ( Usnap' * Dnn(ta34,ta34) * Usnap ) ) ), eeMax ) ;
    end

    case 'eigsUTU'
    [ Vu , du ]= eigs( 1/M .* ( ( Usnap' * Usnap ) ) , eeMax );
    
    case 'InspectAbility'
        % Ability of the modes to represent the FOM solution
       
        H = 1/M .* ( (Usnap' * Dnn(ta34,ta34) * Usnap )+(Usnap' * Dnn(ta34,ta34) * Usnap )' );

        [ Vu , du ]= eig( ( H )) ;

        if isreal( Vu )
                du = (du( end:-1 :end-M+1 , end:-1 :end-M+1 ));
                Vu =  Vu(:,end:-1:end-M+1);
        end
        
        B = Usnap * Vu ;
        for in = 1 : M
            B( : ,in) =  B(: ,in)/norm(B(:,in));
        end
        
        aplRed = ( B'  * SOL.d1([ta3;ta4], IN_QuasiSteady:EN_SnapCaputur ))';
        
        URed = B *  aplRed' ;
%         URed = bsxfun(@plus, URed, REF.d1( [ta3;ta4] , : ) );
      
        varepsilon = 0;
        for in = 1 : M
            varepsilon = varepsilon + norm( Usnap(:,in) - URed(:,in) ) ;     
        end
        
        varepsilon = varepsilon/M;
        fprintf('The average least-square truncation error of %d modes is %d \n',M,varepsilon)
        
        error('not an error')
        
    %   end     da_PROJCT = [ omgTil , aplTil ];
    
    % For the function 'eig' is necessary to invert the outputs' order if
    % not immaginary

end
% Procect empirical eigenvalues onto solutions
B = zeros( TotDofs , eeMax );
B( ta34 , : ) = Usnap * ( Vu / sqrt( ( 1 ) ) );
% B( ta34 , : ) = Usnap * ( Vu / sqrt( ( du ) ) );
for in = 1 : eeMax
    B(ta34,in) =  B(ta34,in)/norm(B(ta34,in));
end

% Check if modes are real
if ~isreal( Vu )
    warning( 'Not all velocity modes are real! ' )
end
% Check if modes are orthogonal

if abCFD_round( sum(diag( B'* Dnn * B ) )/M, 1e-1 ) ~= eeMax
    warning('Velocity modes are not orthogonal')
end

%% COMPUTE TEMPERATURE MODES 
Tsnap = bsxfun( @minus,   SOL.d1( ta1 ,SnapShots ) , REF.d1( ta1 ,1));

fprintf('Start eigenvalue problem at %s with ''%s'' \n', datestr(now,15), ModesFun)
switch ModesFun
    
    case 'eigORT'
        
    Ht =  1/M .* ( (Tsnap' * Dpp( ta1 , ta1 ) * Tsnap)+(Tsnap' * Dpp( ta1 , ta1 ) * Tsnap)' );
%     Ht = 1/2 * 1/M .* ( (Tsnap' * Dpp( ta1 , ta1 ) * Tsnap)+(Tsnap' * Dpp( ta1 , ta1 ) * Tsnap)' );
%     [ Vu , du ]= eig( ( H )) ;
    [ VT , dT ]= eig( Ht );
    % For the function 'eig' is necessary to invert the outputs' order if
    % not immaginary
    if isreal( VT )
            dT = (dT( end:-1 :end-eeMax+1 , end:-1 :end-eeMax+1 ));
            VT =  VT(:,end:-1:end-eeMax+1);
    end

    case 'eigUTU'
    [ VT , dT ]= eig( 1/M .* ( Tsnap'* Tsnap ) );
    % For the function 'eig' is necessary to invert the outputs' order if
    % not immaginary
    if isreal( VT )
            dT = (dT( end:-1 :end-eeMax+1 , end:-1 :end-eeMax+1 ));
            VT =  VT(:,end:-1:end-eeMax+1);
    end

    case 'eigsORT'

    [ VT , dT ]= eigs( 1/M .* ( Tsnap' * Dpp( ta1 , ta1 ) * Tsnap ) , eeMax );

    case 'eigsUTU'
    [ VT , dT ]= eigs( 1/M .* ( Tsnap' * Tsnap ) , eeMax );

end

% Procect empirical eigenvalues onto solutions
BTx = Tsnap * ( VT / sqrt( ( dT ) )  );

% Check if modes are real
if ~isreal( VT )
    warning( 'Not all velocity modes are real! ' )
end
% Check if modes are orthogonal
if abCFD_round( sum(diag( BTx'* Dpp( ta1 , ta1 ) *BTx ) / M ) , 1e-2) ~= eeMax
    warning('Temperature modes are not orthogonal')
end
%% ASSEMBLE MODES

% B = zeros( TotDofs , eeMax );
B( ta1 , : ) = BTx;
%TEST NORMALIZATION
for in = 1 : eeMax
    B(ta1,in) =  B(ta1,in)/norm(B(ta1,in));
end
 
%% PLOT ENERGY PER MODE FOR VELOCITY AND TEMPERATURE FIELDS
figure();
NrgXmode_u = 100 ./ sum(diag( du )) .* diag(du)';
semilogy( NrgXmode_u  ,'x:r','LineWidth',4 );
hold on
NrgXmode_T = 100 / sum( diag(dT) ) * diag(dT)';
semilogy( NrgXmode_T ,'x:b','LineWidth',4 );
title('Energy Percentage per Mode','FontSize', 20,'interpreter','latex' );
plot( [eeMin+.5 eeMin+.5], [ 1e-8 1e2] ,'--g*')
legend( char('Velocity Modes','Temperature Modes',['Limit for ',num2str(sum( NrgXmode_u(1:eeMin)))...
    ,'% of kinetic nrg and ',num2str(sum( NrgXmode_T(1:eeMin))),'% for thermal nrg' ]),'Location','Best')
pause(.005)
axis tight
xlabel('Index of POD modes','FontSize', 10); grid on

%% PLOT ALL MODES FOR VELOCITY TEMPERATURE AND INCOMPRESSIBILITY STATISTICS

% Initialize coefficient matrices to solve the continuity equation locally
% on each element of the mesh
DPHIDXI = [ -1 -1 ; 1 0 ; 0 1];
% Connectivity matrix for Energy Equation leading variable (scalar field)
V_Conn_0 = mesh{2}.ele{1}.dof{ 2  }.TAG;
V_MapVar_0 = reshape( V_Conn_0', 1, 3, []);
% Elements global coordinates [ NumShp x 2 x NumEle ]
COORD = permute(reshape( var.dofCoord( : , permute( V_MapVar_0, [2 1 3] )) , 2, 3,[] ) , [ 2 1 3]);
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

% SOLVE CONTIUITY EQUATION FOR THE FIRST SNAPSHOT

SolTest = SOL.d1( : , IN_QuasiSteady(1) );
Udivergence =  sum( DPHI_X'.*SolTest( mesh{2}.ele{1}.dof{3}.TAG ) , 2) + ...
                    sum( DPHI_Y'.*SolTest( mesh{2}.ele{1}.dof{4}.TAG ) , 2) ;
% Element-wise divergence field
DivField = full( sparse( mesh{2}.ele{1}.dof{3}.TAG( :), linspace(1,1,numel(mesh{2}.ele{1}.dof{3}.TAG( :))), reshape(repmat( Udivergence,1,3)',[],1) ,TotDofs,1) );
%plot divergence field
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 3, cat(2,abs(DivField)) ,1)
hold on
mphgeom( model,'geomFEM', 'Facemode','off','Edgecolor','r')
colorbar
view(-30,28)
title([ 'Diverge of velocity at time t = ',num2str(FemTimeSteps(IN_QuasiSteady(1))),' [s]'],'FontSize',20,'interpreter','latex' );

% PLOT ONLY THE FIRST 'eeMin' MODES
ModesToPlot = eeMax;
WorstDiv = zeros( ModesToPlot , 1 );
MeanDiv  = zeros( ModesToPlot , 1 );
% Incompressibility = zeros( ModesToPlot , 1 );
% InlTag = mesh{1}.ele{ 2 }.dof{ 4 }.TAG;
% OutTag = mesh{1}.ele{ 3 }.dof{ 4 }.TAG;

for id_m = 1 : eeMax
    % check magnitude for velocity modes
    SolTest = B( : , id_m );
    Udivergence =  sum( DPHI_X'.*SolTest( mesh{2}.ele{1}.dof{3}.TAG ) , 2) + ...
                    sum( DPHI_Y'.*SolTest( mesh{2}.ele{1}.dof{4}.TAG ) , 2) ;
    WorstDiv(id_m) = max( abs(Udivergence ));
    MeanDiv(id_m) = mean( abs(Udivergence ));

    if id_m >= eeMin
        figure()
        % Plot magnitude for temperature modes    
        ah1=subplot( 2, 3 , [ 3 6 ]);
        abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, cat(2,B) ,id_m, 'no' )
        title([ 'Temperature magnitude for mode ',num2str(id_m)],'FontSize',15,'interpreter','latex' );   
        view(2)
%         colorbar
        colormap(hot)
        hold on
        mphgeom( model,'geomFEM', 'Facemode','off','Edgecolor','black')
        axis tight
        caxis([ -.03 , 0.03 ])
        freezeColors
        
        ah2=subplot( 2, 3 , [ 1 4 ]);
        abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,B,B) ,id_m,'no' )
        title([ 'Velocity magnitude for mode ',num2str(id_m)],'FontSize',15,'interpreter','latex' );
        view(2)
%         colorbar
        caxis([ 0 , 0.025 ])
        hold on
        colormap(jet)
        mphgeom( model,'geomFEM', 'Facemode','off','Edgecolor','w')
        axis tight
       
        % Plot divergence of velocity modes
        ah3=subplot( 2, 3 , [ 2 5 ]);
%         SolTest = B( : , id_m );
%         Udivergence =  sum( DPHI_X'.*SolTest( mesh{2}.ele{1}.dof{3}.TAG ) , 2) + ...
%                             sum( DPHI_Y'.*SolTest( mesh{2}.ele{1}.dof{4}.TAG ) , 2) ;
%         WorstDiv(id_m) = max( abs(Udivergence ));
%         MeanDiv(id_m) = mean( abs(Udivergence ));
        % Element-wise divergence field
        DivField = full( sparse( mesh{2}.ele{1}.dof{3}.TAG( :), linspace(1,1,numel(mesh{2}.ele{1}.dof{3}.TAG( :))), reshape(repmat( Udivergence,1,3)',[],1) ,TotDofs,1) );
        %plot divergence field
        abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 3, cat(2,abs(DivField)) ,1,'no')
        hold on
        mphgeom( model,'geomFEM', 'Facemode','off','Edgecolor','w')
        freezeColors
        colorbar;
        colormap(jet)
        view(-30,28)
        title( char(['Worst value = ',num2str( WorstDiv(id_m)  ),'           ']),'FontSize',15,'interpreter','latex' );
        view(2)
        caxis([ .0045 , 0.025 ]/1000)
%         caxis( [1e-14 1e-4] )
        freezeColors
        
        posimody=get(ah3,'Position');
        positarget=get(ah1,'Position');
        set(ah3,'Position',[posimody(1),positarget(2:4)])
        
    end
    
end

figure();
semilogy( 1:ModesToPlot, WorstDiv ,'>--m', 'Linewidth', 2,'MarkerEdgeColor','b'...
                        ,'MarkerFaceColor','m','MarkerSize',9)
hold on
semilogy( 1:ModesToPlot, MeanDiv ,'--^c', 'Linewidth', 2,'MarkerEdgeColor','b'...
                        ,'MarkerFaceColor','c','MarkerSize',9)
% semilogy( 1:ModesToPlot, Incompressibility , 'xw')
title( 'Residual on Continuity equation for each mode','FontSize',15,'interpreter','latex' );   
xlabel('Mode number [-]','FontSize',12)
ylabel('Residual','FontSize',12)
axis tight
grid on
legend( char('Worst value','Mean value'),'Location','Best')
pause(.5)
%% CLOSING 

fprintf('Modes computed in %f minutes \n', (toc-tiin)/60 )

%%   =======================================================================  %%
%% - -  A S S E M B L E    R E D U C E D    O R D E R    S Y S T E M   - - - - %
%%   =======================================================================  %%


% Simulation time steps from Comsol
tSteps = SOL.solinfo.solvals( IN_QuasiSteady : end );

% Function handles to evaluate mass flow rate and heat flow rate in time
MInT = @(t) t./t;     % Constant and equal to one
%                       SINUIDAL
% MInT = @(t) sin(  (t - tSteps(1)) * (2*pi/(tSteps(end)-tSteps(1))  )); 

HInT = @(t) t./t;     % Constant and equal to one

% Material properties
kkon = ph.mat.k; % = kon
kref = ph.mat.k; % = kon.reference
expa = ph.mat.thex; % thermal expansion coeff
rhoo = ph.mat.rho; % density
Tblk = 0; % Bulk temperature
mudd = ph.mat.mud; % Dynamic viscosity
lamb = kkon / ( rhoo * ph.mat.Cp );

% Initialize reference solution
RT = zeros( TotDofs , 1 );
RT( ta1, : ) = REF.d1( ta1, : );
Ru = zeros( TotDofs , 1 );
Ru( [ ta3; ta4 ] , : ) = REF.d1( [ ta3; ta4 ] , : );

% Loop using different number of modes
for ee = eeMin : 1 : eeMax

    % Assemble mode
    BT = zeros( TotDofs, ee );
    BT(ta1, : ) = B(ta1, 1:ee );
    Bu = zeros( TotDofs, ee );
    Bu( [ ta3; ta4 ], : ) = B( [ ta3; ta4 ], 1:ee );
 
    %% C O M P U T E   R E D U C E D   O R D E R    O P E R A T O R S
    % Do the projection of the modes to obtain the reduced operators. Each
    %   governing equation brings reduced operators, and these will be computed
    %   in the following.


    %% COMPUTE FIRST ORDER REDUCED OPERATORS
    % initialize
    tiin = toc;
    % On energy equation
    BTDTRT      = BT' *    Dpp      * RT    ;
    BTVrefRu    = BT' * VrefTppn    * Ru    ;

    % On momentum equations
    BuDuRu      = Bu' *     Dnn     *   Ru  ;
    BuCrefRu    = Bu' *  Crefnnn    *   Ru  ;
    BuGTRT      = Bu' *     Dnp     *   RT  ;
    BuPu        = Bu' *     Pn              ;

    fprintf(' FIRST ORDER OPERATORS computed in %f minutes \n', (toc-tiin)/60 )

    %% COMPUTE SECOND ORDER COEFFICIENT MATRICES
    % Initialize
    tiin = toc;
    if ee == eeMin
        start = 1;
        
    CpodBu = cell( 1 , ee ); 
    CpodBu= cellfun(@(x) sparse(TotDofs,TotDofs),CpodBu,'UniformOutput',false);
    VpodBT = cell( 1 , ee );
    VpodBT = cellfun(@(x) sparse(TotDofs,TotDofs),VpodBT ,'UniformOutput',false);
 
    else
        start = ee;
    end
    
    
    for l = start : ee

    [ CpodBu{ l } , VpodBT{ l } ] =  abCFD_Assmbl_Convective( B( : , l ), mesh, var, 'DontCHECK' );
     
    end
  
        % On energy equation
    BTDTBT      = BT' *     Dpp     * BT    ;
    BTVTpodRu   = cell2mat( cellfun(@(x) BT' * x * Ru, VpodBT ,'UniformOutput',0) );%CORRECT! TESTD with 'target'
    BTKTBT      = BT' *   Kpp_jj    * BT    ;
    BTVrefBu    = BT' *  VrefTppn   * Bu    ;

        % On momentum equations
    BuDuBu      = Bu' *     Dnn     * Bu    ;
    BuCrefBu    = Bu' *  Crefnnn    * Bu    ;
    BuCpodRu    = cell2mat( cellfun(@(x) Bu' * x * Ru, CpodBu ,'UniformOutput',0) );%CORRECT! TESTD with 'target'
    BuKuBu      = Bu' *   Knn_jj    * Bu    ;
    BuKvBv      = Bu' *   Knn_ij    * Bu    ;
    BuGTBT      = Bu' *   Dnp       * BT    ;
    
    fprintf(' SECOND ORDER OPERATORS computed in %f minutes \n', (toc-tiin)/60 )

    %% COMPUTE THIRD ORDER COEFFICIENT MATRICES 
    % Initialize
    tiin = toc;

        % On energy equation
    BTVTpodBu    = reshape( full( cell2mat( cellfun(@(x) BT' * x * Bu, VpodBT ,'UniformOutput',0) ) ) , ee,ee,ee); % TESTED!
        % On momentum equations
    BuCpodBu    =  reshape( full( cell2mat( cellfun(@(x) Bu' * x * Bu, CpodBu ,'UniformOutput',0) ) ) , ee,ee,ee) ;

    fprintf(' THIRD ORDER OPERATORS computed in %f minutes \n', (toc-tiin)/60 )
    
%%   =======================================================================  %%
%% - - - - - -    S O L V E   R E D U C E D    O R D E R    M O D E L    - - - %
%%   =======================================================================  %%
    % Solve the Reduced order system with different solvers
    
    %% DATA COMMON TO ALL SOLVERS
    % Initialize
    tiin = toc;

    % Initial conditions for Reduced Model unkonwns
    aInit = [ BTDTBT \ ( BT' * Dpp * bsxfun( @minus, SOL.d1( :, IN_QuasiSteady ), RT( : , : ) ) ) ;  ...
              BuDuBu \ ( Bu' * Dnn * bsxfun( @minus, SOL.d1( :, IN_QuasiSteady ), Ru( : , : ) ) ) ] ;
    
%     aInit = [ ( BT' * Dpp * bsxfun( @minus, SOL.d1(: , IN_QuasiSteady ), REF.d1 )  ) ;
%               ( Bu' * Dnn * bsxfun( @minus, SOL.d1(: , IN_QuasiSteady ), REF.d1 )  ) ];
         
    
    % Idendity tags of unknowns for different solvers
    omg = 1 : ee;       % Energy equation coeff
    alp = ee+1 : 2*ee;  % Momentum equation coeff.

    
%% USE BUILT-IN SOLVER   ODE45
    % Define product function for third order operators

Prod3D = @(Q_klm, a_l, a_m) reshape( reshape( Q_klm , ee*ee, ee) * a_l , ee ,ee) * a_m;

da45 = @(t,a) [-(  0                                * BTDTRT    *   1     + ...
                  HInT(t)*( MInT(t) - kkon/kref )   * BTVrefRu  *   1     + ...
                  MInT(t)                           * BTVTpodRu * a(omg)  + ...
                  HInT(t)                           * BTVrefBu  * a(alp)  + ...
    reshape( reshape( BTVTpodBu , ee*ee, ee) * a(omg) , ee ,ee) * a(alp)  + ...
HInT(t)                           * BTKTBT    * a(omg) )...
                   ;
                 -( 0                                * BuDuRu    *   1     + ...
                   MInT(t).^2                       * BuCrefRu  *   1     + ...
                   MInT(t)                          * BuCrefBu  *  a(alp) + ...
                   MInT(t)                          * BuCpodRu  *  a(alp) + ...
                  Prod3D( BuCpodBu , a(alp), a(alp) )           *   1     + ...
                  mudd/rhoo                         * BuKuBu    *  a(alp) + ...
                  mudd/rhoo                         * BuKvBv    *  a(alp) + ...
                 - 0                                * BuGTRT    *   1     + ...
                 - 0                                * BuGTBT    *  a(omg) + ...
                 + 0                                * BuPu      *   1       ) ];
% DAMPING = [ BTDTBT , zeros(ee,ee); zeros(ee,ee), BuDuBu ];
DAMPING = diag( [diag(BTDTBT); diag(BuDuBu )]);

tollerance = 1e-7;
opts = odeset('RelTol',tollerance,'AbsTol',tollerance*10,'Stats','off','Mass',DAMPING);
 
[tSteps_ode45 , da_ODE45 ] = ode45( da45 , tSteps , aInit ,opts);
  
    fprintf('S O L V E R:   ode45   for %d modes finished in %f sec. \n',ee,(toc-tiin))

%% USE BUILT-IN SOLVER   ODE15i
 
Prod3D = @(Q_klm, a_l, a_m) reshape( reshape( Q_klm , ee*ee, ee) * a_l , ee ,ee) * a_m;

res15i = @(t,a,da) [ (  0                           * BTDTRT    *   1     + ...
                   1                                * BTDTBT    * da(omg) +...
                  HInT(t)*( MInT(t) - kkon/kref )   * BTVrefRu  *   1     + ...
                  MInT(t)                           * BTVTpodRu * a(omg)  + ...
                  HInT(t)                           * BTVrefBu  * a(alp)  + ...
    reshape( reshape( BTVTpodBu , ee*ee, ee) * a(omg) , ee ,ee) * a(alp)  + ...
                  lamb                              * BTKTBT    * a(omg) )...
                   ;
                 ( 0                                * BuDuRu    *   1     + ...
                   1                                * BuDuBu    * da(alp) + ...
                   MInT(t).^2                       * BuCrefRu  *   1     + ...
                   MInT(t)                          * BuCrefBu  *  a(alp) + ...
                   MInT(t)                          * BuCpodRu  *  a(alp) + ...
                  Prod3D( BuCpodBu , a(alp), a(alp) )           *   1     + ...
                  mudd/rhoo                         * BuKuBu    *  a(alp) + ...
                  mudd/rhoo                         * BuKvBv    *  a(alp) + ...
                 - 0                                * BuGTRT    *   1     + ...
                 - 0                                * BuGTBT    *  a(omg) + ...
                 + 0                                * BuPu      *   1       ) ];
             
opts = odeset('RelTol',tollerance,'AbsTol',tollerance*10,'Stats','off');

% Compute consistent initial conditions.
da = zeros( size( aInit ) );
[y0,yp0] = decic( res15i ,tSteps(1),aInit,[],da,[],opts);

% Solve the problem.
[tSteps_ode15i , da_ODE15i ] = ode15i( res15i ,tSteps', y0, yp0,opts);

fprintf('S O L V E R:   ode15i   for %d modes finished in %f sec. \n',ee,(toc-tiin))
 
%% USE THETA-FAMILY SCHEME FOR ALPHA AND OMEGA MONOLITIC APPROACH
    tiin = toc;
    % Select scheme
    tht = 1;  % For tht = 1 we have fully implicit method
    % Initialize LHS and RHS
    LHS = zeros( 2*ee , 2*ee );
    RHS = zeros( 2*ee ,   1  );
    % Initialize unknowns
    da_MONOi = zeros( 2*ee , numel( tSteps ) );
    da_MONOi( : , 1 ) = aInit;
    
    for id_t = 2 : 1 : numel( tSteps )
        
        dt = tSteps( id_t , 1 ) - tSteps( id_t-1 , 1 );
      % ASSEMBLT LEFT-HAND SIDE
        %  f i r s t    q u a d r a n t s 
    LHS( 1:ee, 1 : ee ) =           1/dt .* ( BTDTBT  )                     + ...
                          ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * tht * BTVTpodRu + ...
                          lamb * BTKTBT * tht ;
        %  S e c o n d    q u a d r a n t s
    LHS(    1:ee     , 1+ee : 2*ee ) =                         ...
                         (tht * HInT(id_t) + ( 1-tht) * HInT(id_t-1) ) * tht * BTVrefBu  + ...
                         reshape( reshape( BTVTpodBu , ee*ee, ee ) * da_MONOi(  1  :  ee  , id_t-1 ) , ee ,ee);
        %  t h i r d    q u a d r a n t s 
    LHS( 1+ee : 2*ee ,     1:ee    ) =  0 * tht * BuGTBT                     ;
        %  F o u r t h    q u a d r a n t s 
    LHS( 1+ee : 2*ee , 1+ee : 2*ee ) =           1/dt .* ( BuDuBu  )        + ...
                         + ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * tht * BuCrefBu  + ...
                         + ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * tht * BuCpodRu  + ...
                         + mudd / rhoo * tht * BuKuBu                        + ...
                         + mudd / rhoo * tht * BuKvBv                        + ...
 + reshape( reshape( BuCpodBu , ee*ee, ee) * da_MONOi( 1+ee : 2*ee , id_t-1 ) , ee ,ee) ;

    % ASSEMBLE RIGHT-HAND SIDE
        % E n e r g y   e q u a t i o n
    RHS(    1:ee     , 1 ) = - ( HInT(id_t) - HInT(id_t-1) ) / dt * BTDTRT                + ...
  - ( tht * HInT(id_t) +(1-tht) * HInT(id_t-1) ) * ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) -kkon/kref) * BTVrefRu  + ...
  + 1/ dt * BTDTBT * da_MONOi(  1  :  ee  , id_t-1 )                            + ...
  - ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * BTVTpodRu * (1-tht) * da_MONOi(   1  :  ee  , id_t-1 ) + ...
  - lamb * BTKTBT * (1-tht) * da_MONOi(  1  :  ee  , id_t-1 )                   + ...             
  - ( tht * HInT(id_t) +(1-tht) * HInT(id_t-1) ) * (1-tht) * BTVrefBu * da_MONOi( 1+ee : 2*ee , id_t-1 ) + ...
 - (1-tht) * reshape( reshape( BTVTpodBu , ee*ee, []) * da_MONOi(  1  :  ee  , id_t-1 ) , ee ,ee) * da_MONOi( 1+ee : 2*ee , id_t-1 ) ;

        % M o m e n t u m    e q u a t i o n
    RHS( 1+ee : 2*ee , 1 ) = - ( MInT(id_t) - MInT(id_t-1) ) / dt * BuDuRu  + ...
  - ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) )^2 * BuCrefRu              + ...
  + 0 * ( tht * HInT(id_t) +(1-tht) * HInT(id_t-1) ) * BuGTRT               + ...
  - 0 * Tblk * BuPu                                                         + ...
  + 1/ dt * BuDuBu * da_MONOi( 1+ee : 2*ee , id_t-1 )                           + ...
  - ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * (1-tht) * BuCrefBu * da_MONOi( 1+ee : 2*ee , id_t-1 ) + ...
  - ( tht * MInT(id_t) + (1-tht) * MInT(id_t-1) ) * (1-tht) * BuCpodRu * da_MONOi( 1+ee : 2*ee , id_t-1 ) + ...
  - mudd/rhoo * (1-tht) * BuKuBu * da_MONOi( 1+ee : 2*ee , id_t-1 )             + ...
  - mudd/rhoo * (1-tht) * BuKvBv * da_MONOi( 1+ee : 2*ee , id_t-1 )             + ...
  + 0 * (1-tht) * BuGTBT * da_MONOi(   1  :  ee  , id_t-1 )                     + ...
  - (1-tht) * reshape( reshape( BuCpodBu , ee*ee, ee) * da_MONOi( 1+ee : 2*ee , id_t-1 ) , ee ,ee) * da_MONOi( 1+ee : 2*ee , id_t-1 ) ;
     
    da_MONOi( 1 : 2*ee , id_t ) = LHS \ RHS;
    
    end % End time loop
 
    fprintf('S O L V E R:   theta   for %d modes finished in %f sec. \n',ee,(toc-tiin))
    
    
    
  
%% RECONSTRUCT PHYSICAL SOLUTION
    if ee > 5 
        
     % FROM SOLVER ODE45    
    Upod_ODE45 = BT *  da_ODE45(:,omg)' ;
    Upod_ODE45 = Upod_ODE45 + Bu *  da_ODE45(:,alp)' ;
    Upod_ODE45 = bsxfun(@plus, Upod_ODE45, REF.d1(:,end) );
 
    % FROM SOLVER ODE    i
    Upod_ODE15 = BT *  da_ODE15i(:,omg)' ;
    Upod_ODE15 = Upod_ODE15 + Bu *  da_ODE15i(:,alp)' ;
    Upod_ODE15 = bsxfun(@plus, Upod_ODE15, REF.d1(:,end) );
    

     FG_geom = figure;
    
    mphgeom(model); hold on

    probes =[ 29452, 32040 , 51632 , 72501, 59741];
    ErrBar = zeros( numel(probes) , 4);
    
        for i_prob = 1 : numel( probes )
            U_prob=mesh{2}.ele{1}.dof{3}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{3}.uniqTAG) ==  ...
                                        var.dofCoord(1,probes(i_prob))) ); %#ok<*FNDSB>
            V_prob=mesh{2}.ele{1}.dof{4}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{4}.uniqTAG) ==  ...
                                        var.dofCoord(1,probes(i_prob))) );
            T_prob=mesh{2}.ele{1}.dof{1}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{1}.uniqTAG) ==  ...
                                        var.dofCoord(1,probes(i_prob))) );
                figure(FG_geom)
            plot( var.dofCoord(1,U_prob), var.dofCoord(2,U_prob) ,'sr',...
                'MarkerSize',10,'LineWidth',2)
            plot( var.dofCoord(1,V_prob), var.dofCoord(2,V_prob) ,'xb',...
                'MarkerSize',10,'LineWidth',2)
            plot( var.dofCoord(1,T_prob), var.dofCoord(2,T_prob) ,'+y',...
                'MarkerSize',10,'LineWidth',2)


            text( var.dofCoord(1,U_prob), var.dofCoord(2,U_prob) , num2str( i_prob )  , ...
                    'VerticalAlignment'     ,   'bottom'   ,...
                    'HorizontalAlignment'   ,   'right' ,...
                    'color'     , 'k'   , 'FontWeight', 'bold',...
                    'Margin',1,'FontSize',20);
            axis tight
            xlabel(['x-coord [m]'],'FontSize',12)
            ylabel(['y-coord [m]'],'FontSize',12)
 

             
            ErrBar( i_prob , ee ) = ...
                norm( ( SOL.d1(V_prob, IN_QuasiSteady:end )-Upod_ODE45(V_prob, : )));

            StoreErrPlot_Ux{i_prob}(ee,:)=Upod_ODE45(U_prob, : );
            StoreErrPlot_Uy{i_prob}(ee,:)=Upod_ODE45(V_prob, : );
            StoreErrPlot_T{i_prob}(ee,:)=Upod_ODE45(T_prob, : );
            
                figure()
                subplot(311)
                title(['Probe # ',num2str(i_prob),' for ',num2str(ee),' modes'],'FontSize',15)

                hold all; grid on
                plot(  SOL.d1(U_prob, IN_QuasiSteady:end ) ,'LineWidth',4)
%                 plot(  Upod_ODE45(U_prob, : ), '--' , 'LineWidth',2)  
                ylabel([' Ux [m/s]'],'FontSize',10)

                subplot(312)
                hold all; grid on
                plot(  SOL.d1(V_prob, IN_QuasiSteady:end ) ,'LineWidth',4)
%                 plot(  Upod_ODE45(V_prob, : ), '--' , 'LineWidth',2) 
                ylabel([' Uy [m/s]'],'FontSize',10)

                subplot(313)
                hold all; grid on
                plot(  SOL.d1(T_prob, IN_QuasiSteady:end ) ,'LineWidth',4)
%                 plot(  Upod_ODE45(T_prob, : ), '--' , 'LineWidth',2) 
                ylabel([' T [?]'],'FontSize',10)
%                  legend(char('Full Model','Garlerkin-POD'))

        end
    elseif ee == 10 || ee == 15 || ee == 25
    
    if ee >=1
    % FROM SOLVER ODE45    
    Upod_ODE45 = BT *  da_ODE45(:,omg)' ;
    Upod_ODE45 = Upod_ODE45 + Bu *  da_ODE45(:,alp)' ;
    Upod_ODE45 = bsxfun(@plus, Upod_ODE45, REF.d1(:,end) );
 
    % FROM SOLVER ODE    i
    Upod_ODE15 = BT *  da_ODE15i(:,omg)' ;
    Upod_ODE15 = Upod_ODE15 + Bu *  da_ODE15i(:,alp)' ;
    Upod_ODE15 = bsxfun(@plus, Upod_ODE15, REF.d1(:,end) );
 
    
%     % MONOLITIC IMPLICIT SOLVER
%     Upod_MONOi = BT *  da_MONOi(omg,:) ;
%     Upod_MONOi = Upod_MONOi + Bu *  da_MONOi(alp,:) ;
%     Upod_MONOi = bsxfun(@plus, Upod_MONOi, REF.d1(:,end) );
%     
    
    
    for i_prob = 1 : numel( probes )
    U_prob=mesh{2}.ele{1}.dof{3}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{3}.uniqTAG) ==  ...
                                var.dofCoord(1,probes(i_prob))) ); %#ok<*FNDSB>
    V_prob=mesh{2}.ele{1}.dof{4}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{4}.uniqTAG) ==  ...
                                var.dofCoord(1,probes(i_prob))) );
    T_prob=mesh{2}.ele{1}.dof{1}.uniqTAG( find( var.dofCoord( 1 , mesh{2}.ele{1}.dof{1}.uniqTAG) ==  ...
                                var.dofCoord(1,probes(i_prob))) );
    
    ErrBar( i_prob , ee ) = ...
        norm( ( SOL.d1(V_prob, IN_QuasiSteady:end )-Upod_ODE45(V_prob, : )));

%         
%     figure( FG_errorV )
%     markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
%     subplot(5,1,i_prob)
%     semilogy( tSteps , ...
%   (100*abs( SOL.d1(V_prob, IN_QuasiSteady:end )-Upod_ODE45(V_prob, : )))./ (( SOL.d1(V_prob, IN_QuasiSteady:end ))  ), ...
%         [':',markers{ee}],'LineWidth',2, 'MarkerSize',8)
%     hold all
%     axis([ min(tSteps), max(tSteps) ,1e-10, 1 ])
%     grid on
    
    
            StoreErrPlot_Ux{i_prob}(ee,:)=Upod_ODE45(U_prob, : );
            StoreErrPlot_Uy{i_prob}(ee,:)=Upod_ODE45(V_prob, : );
            StoreErrPlot_T{i_prob}(ee,:)=Upod_ODE45(T_prob, : );

    
   end
    legend(char('lowR','lowC','medC','topC','topL'))
 

    end
    
    elseif ee==1
%%   =======================================================================  %%
%% - - - - - - - - - -    A N A L Y I Z E     R E S U L T S    - - - - - - - - %
%%   =======================================================================  %%
%% INITIALIZE

        %  PROJECT HOMOGENEOUS SOLUTION ONTO REDUCED BASIS
%        omgTil = ( BT' * Dpp *  bsxfun( @minus, SOL.d1(: , IN_QuasiSteady:end ), REF.d1 )  )';
%        aplTil = ( Bu' * Dnn *  bsxfun( @minus, SOL.d1(: , IN_QuasiSteady:end ), REF.d1 )  )';
       omgTil = ( BT'  *  bsxfun( @minus, SOL.d1(: , IN_QuasiSteady:end ), REF.d1 )  )';
       aplTil = ( Bu'  *  bsxfun( @minus, SOL.d1(: , IN_QuasiSteady:end ), REF.d1 )  )';
       da_PROJCT = [ omgTil , aplTil ];

        %% COMPARE MONO vel ODE45
     figure()
     subplot(3,1,[1 2])
     hold on;plot(linspace(0,27.6,numel(tSteps)), aplTil(:,1),'-go','LineWidth',2)
     plot(linspace(0,27.5,numel(tSteps)), da_ODE45(:,ee+1),'-kx','LineWidth',.02)
     grid on
     ylabel('\alpha_{ 1}','FontSize',14)
     legend( char('Projection','Galerkin-POD'),'Location','Best')
     axis( [0 25 -70 70])
     
     subplot(3,1,3)
     plot(linspace(0,26,numel(tSteps)), aplTil(:,1)-  da_ODE45(:,ee+1),'-b','LineWidth',.02)
     ylabel('\alpha_{ 1}^{PRJ} - \alpha_{ 1}^{POD}','FontSize',14)
     grid on
     axis( [0 25 -70 70])
     
         errMONO= sqrt(mean( abs( aplTil(:,1))-abs(  da_MONOi(ee+1,:)' )).^2 );
  fprintf('\nAverage time integrated relative errMONO of %d modes is %d \n',ee,errMONO)
     
     errODE= sqrt(mean( abs( aplTil(:,1))-abs(  da_ODE45(:,ee+1) )).^2 );
  fprintf('Average time integrated relative errODE of %d modes is %d  \n',ee,errODE)

       errODEi= sqrt(mean( abs( aplTil(:,1))-abs(  da_ODE15i(:,ee+1) )).^2 );
  fprintf('Average time integrated relative errODE of %d modes is %d \n\n',ee,errODEi)

%% ERROR ON VELOCITY

%         varepsilon = 0;
%         for in = 1 : numel( tSteps )
%             varepsilon = varepsilon + norm( SOL.d1(: , IN_QuasiSteady + in-1 )...
%                                                 - Upod_ODE15(:,in) ) ;     
%         end        
%         varepsilon = varepsilon/ee;
%         fprintf('\n Average square truncation error of %d modes is %d \n\n',ee,varepsilon)
        
        
%% FEW EMPIRICAL COEFFICIENTS EVOLUTION IN TIME
    if ee == eeMin % run for first ee-cycle
        
       FG_TIME_Alp = figure(); 

       hold on
       tildLines = plot( tSteps',   aplTil( : ,  1:4   )'   , 'r', 'LineWidth', 2)     ;
       MonoLines = plot( tSteps', da_MONOi( ee+1:ee+1+3 , :)', 'm','LineWidth', 1.1 );
       ode15Line = plot( tSteps, da_ODE15i( : , ee+1:ee+1+3) , 'b','LineWidth', 1.1 );
       ode45Line = plot( tSteps, da_ODE45( : , ee+1:ee+1+3),'-.g','LineWidth', 1.1 );
       axis tight       
       title('First 4 Velocity Empirical Coefficients Comparison','FontSize', 20,'interpreter','latex' );
       xlabel('Time [s]')
       ylabel('Empirical coeff. values')
       grid on
       TildGroup = hggroup;
       MonoGroup = hggroup;
       ode15Group = hggroup;
       ode45Group = hggroup;
       set(tildLines,'Parent',TildGroup)
       set(MonoLines,'Parent',MonoGroup)
       set(ode15Line,'Parent',ode15Group)
       set(ode45Line,'Parent',ode45Group)
        % Include these hggroups in the legend:
        set(get(get(TildGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on'); 
        set(get(get(MonoGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on'); 
        set(get(get(ode15Group,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on'); 
        set(get(get(ode45Group,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on'); 
        lgn_fewcomp_u = char('ModesProjection','MONOi Solver','ODE15i Solver','ODE45 Solver');
        legend(lgn_fewcomp_u,'Location','Best')
       
       

       FG_TIME_Omg = figure(); 
       hold all
       plot( tSteps,   omgTil( : ,  1:4   )   , 'r', 'LineWidth', 2)  
       plot( tSteps, da_MONOi( 1:4 , :), 'm','LineWidth', 1.1 )
       plot( tSteps, da_ODE15i(  :,1:4 ), 'b','LineWidth', 1.1 )
       plot( tSteps, da_ODE45( :,1:4), '-.g','LineWidth', 1.1 )
       axis tight       
%        legend('ModesProjection','MONOi Solver','SEGRi solver','ODE23s Solver','Location','Best')
       title('First 4 Temperature Empirical Coefficients Comparison','FontSize', 20,'interpreter','latex' );
       xlabel('Time [s]')
       ylabel('Empirical coeff. values')
       grid on

    else
        
       figure( FG_TIME_Alp )
       pause(0.5)
       plot( tSteps, da_MONOi( ee+1:ee+1+3 , :), 'm','LineWidth', 1.1 )
       plot( tSteps, da_ODE15i( : , ee+1:ee+1+3) , 'b','LineWidth', 1.1 )
       plot( tSteps, da_ODE45( : , ee+1:ee+1+3),'-.g','LineWidth', 1.1 )
       pause(.5)

       figure( FG_TIME_Omg )
       pause(0.5)
       plot( tSteps, da_MONOi( 1:4 , :), 'm','LineWidth', 1.1 )
       plot( tSteps, da_ODE15i(  :,1:4 ), 'b','LineWidth', 1.1 )
       plot( tSteps, da_ODE45( :,1:4), '-.g','LineWidth', 1.1 )
       pause(0.5)
       
    end

    pause( 6 )
%% EMPIRICAL COEFFICIENTS TIME INTEGRAL
% Project solution onto modes

    if ee == eeMin % run for first ee-cycle
       FG_INT_Alp = figure(); 
       hold all
       plot( tSteps, sum(   aplTil( : ,  1:end   ) ,2) , 'xr' )
       plot( tSteps, sum( da_MONOi( 1+ee:end , :) ,1) , ':', 'LineWidth', 3 )
%        plot( tSteps, sum( da_SEGRi( : , :)        ,1) , '.-', 'LineWidth', 1.5 )

       axis tight
%        plot( tSteps( 1 : size(da_ODE23,1) ), sum( da_ODE23( : , 1+ee:end ) ,2) , '-', 'LineWidth', 1.5 )
%        legend('ModesProjection','Implicit Solver','ODE23s Solver','Location','Best')
       title('Velocity Empirical Coefficients - Time Integral','FontSize', 20,'interpreter','latex' );
       xlabel('Time [s]')
       ylabel('Empirical coeff. values')
       grid on
    else
        figure( FG_INT_Alp )
        pause(.5)
       plot( tSteps, sum(   aplTil( : ,  1:end   ) ,2) , 'xr' )
       plot( tSteps, sum( da_MONOi( 1+ee:end , :) ,1) , ':', 'LineWidth', 3 )
%        plot( tSteps, sum( da_SEGRi( : , :)        ,1) , '.-', 'LineWidth', 1.5 )

%        plot( tSteps( 1 : size(da_ODE23,1) ), sum( da_ODE23( : , 1+ee:end ) ,2) , '-', 'LineWidth', 1.5 )

    end

    pause( 5 )
    
%% RELATIVE ERROR OF EMPIRICAL COEFFICENT IN TIME

if ee == eeMin
    
    FG_RelErrTime = figure();
    
    CoefErr_T = sum( abs( da_PROJCT(:,omg) - da_ODE45(:,omg)  ), 2 ) ./ (100 - (NrgXmode_T(ee)));
    CoefErr_u = sum( abs( da_PROJCT(:,alp) - da_ODE45(:,alp)  ), 2 ) ./ (100 - (NrgXmode_u(ee)));

    semilogy( tSteps, CoefErr_T , 'LineWidth', 3 )
    hold all
    semilogy( tSteps, CoefErr_u , 'LineWidth', 3 )

    grid on
    axis tight
    title('RELATIVE ERROR OF EMPIRICAL COEFFICENT IN TIME','FontSize', 20,'interpreter','latex' );
    xlabel('time [s]')
    ylabel('sum( abs( da_PROJCT - da_ODE45) ) ./ (100 - (NrgXmode_T)')
    LG_thermal_1 = char(['Thermal rel.err with ',num2str(ee),' modes'],['Momentum rel.err with ',num2str(ee),' modes']);
    legend( LG_thermal_1 )
    grid on
    pause( 2 )

else
    figure( FG_RelErrTime );
    
    CoefErr_T = sum( abs( da_PROJCT(:,omg) - da_ODE45(:,omg)  ), 2 ) ./ (100 - (NrgXmode_T(ee)));
    CoefErr_u = sum( abs( da_PROJCT(:,alp) - da_ODE45(:,alp)  ), 2 ) ./ (100 - (NrgXmode_u(ee)));

    semilogy( tSteps, CoefErr_T , 'LineWidth', 3 )
    hold all
    semilogy( tSteps, CoefErr_u , 'LineWidth', 3 )
    
    LG_thermal_1 = char( LG_thermal_1, ['Thermal rel.err with ',num2str(ee),' modes'],['Momentum rel.err with ',num2str(ee),' modes']);
    legend( LG_thermal_1 )
    pause( 2 )
    
end


%% POD-COEFFICIENT MUTUAL REPRESENTATION
  
        
        % THREE
        Plotinterval = 1:100;
%         for id_mode = 1 : ee : ee+2 
                id_mode =  1;
%             id_mode =  4;
            figure()
            subplot(121)
            hold on
            % Projected values
            plot3( da_PROJCT(Plotinterval,id_mode+1) , da_PROJCT(Plotinterval,id_mode+1+1) , ...
                da_PROJCT(Plotinterval,id_mode+1+2) , 'og','MarkerSize',12,'LineWidth', 3 )
            % Rom values

            plot3( da_ODE45(:,id_mode+1)  , da_ODE45(:,id_mode+1+1)  , da_ODE45(:,id_mode+1+2)  ,'-', 'LineWidth', .5 )
            plot3( da_ODE15i(:,id_mode+1)  , da_ODE15i(:,id_mode+1+1)  , da_ODE15i(:,id_mode+1+2)  ,'--r', 'LineWidth', .5 )
            grid on
            
            xlabel(['\omega _{',num2str(id_mode+1),'}'],'FontSize',15)
            ylabel(['\omega _{',num2str(id_mode+1+1),'}'],'FontSize',15)
            zlabel(['\omega _{',num2str(id_mode+1+2),'}'],'FontSize',15)
%             title('\alpha. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','ODE45e','ODE15i'),'Location','Best')
            view( 3 )
            pause( .5 )
            axis equal
            
            
            id_alp = id_mode+1+1;
            subplot(122)
            hold on
            % Projected values
            plot3( da_PROJCT(Plotinterval,id_alp) , da_PROJCT(Plotinterval,id_alp+1) , ...
                da_PROJCT(Plotinterval,id_alp+2) , 'og','MarkerSize',12,'LineWidth', 3 )
            % Rom values

            plot3( da_ODE45(:,id_alp)  , da_ODE45(:,id_alp+1)  , da_ODE45(:,id_alp+2)  ,'-', 'LineWidth', .5 )
            plot3( da_ODE15i(:,id_alp)  , da_ODE15i(:,id_alp+1)  , da_ODE15i(:,id_alp+2)  ,'--r', 'LineWidth', .5 )
            grid on
            
            xlabel(['\alpha _{',num2str(id_mode),'}'],'FontSize',15)
            ylabel(['\alpha _{',num2str(id_mode+1),'}'],'FontSize',15)
            zlabel(['\alpha _{',num2str(id_mode+2),'}'],'FontSize',15)
%             title('\alpha. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','ODE45e','ODE15i'),'Location','Best')
            view( 3 )
            pause( .5 )
            axis equal
%         end
        
        da_MONOi = da_MONOi';
           id_mode =  1;
%             id_mode =  4;
            figure()
            subplot(121)
            hold on
            % Projected values
            plot3( da_PROJCT(Plotinterval,id_mode+1) , da_PROJCT(Plotinterval,id_mode+1+1) ,...
                da_PROJCT(Plotinterval,id_mode+1+2) , 'og','MarkerSize',12,'LineWidth', 3 )
            % Rom values

            plot3( da_MONOi(:,id_mode+1)  , da_MONOi(:,id_mode+1+1)  , da_MONOi(:,id_mode+1+2)  ,'r-', 'LineWidth', 2)            
            grid on
            
            xlabel(['\omega _{',num2str(id_mode+1),'}'],'FontSize',15)
            ylabel(['\omega _{',num2str(id_mode+1+1),'}'],'FontSize',15)
            zlabel(['\omega _{',num2str(id_mode+1+2),'}'],'FontSize',15)
%             title('\alpha. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Euler Implicit'),'Location','Best')
            view( 3 )
            pause( .5 )
            axis equal
            
            
            id_alp = id_mode+1+1;
            subplot(122)
            hold on
            % Projected values
            plot3( da_PROJCT(Plotinterval,id_alp) , da_PROJCT(Plotinterval,id_alp+1) ,...
                da_PROJCT(Plotinterval,id_alp+2) , 'og','MarkerSize',12,'LineWidth', 3 )
            % Rom values

            plot3( da_MONOi(:,id_alp)  , da_MONOi(:,id_alp+1)  , da_MONOi(:,id_alp+2)  ,'r-', 'LineWidth', 2 )
            grid on
            
            xlabel(['\alpha _{',num2str(id_mode),'}'],'FontSize',15)
            ylabel(['\alpha _{',num2str(id_mode+1),'}'],'FontSize',15)
            zlabel(['\alpha _{',num2str(id_mode+2),'}'],'FontSize',15)
%             title('\alpha. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Euler Implicit'),'Location','Best')
            view( 3 )
            pause( .5 )
            axis equal
  
        da_MONOi = da_MONOi';
        
        
        % TWO
    if ee == 1
        for id_mode = 1 : 2 : ee - 2

            figure()
            subplot(121)
            hold on
            % Projected values
            plot( da_PROJCT(:,id_mode) , da_PROJCT(:,id_mode+1) , 'og' )
            % Rom values
            plot( da_ODE45(:,id_mode)  , da_ODE45(:,id_mode+1)   ,'-', 'LineWidth', 1 )
            plot( da_ODE15i(:,id_mode)  , da_ODE15i(:,id_mode+1) ,'--r', 'LineWidth', 1 )
            grid on
            xlabel(['POD coeff. ',num2str(id_mode)])
            ylabel(['POD coeff. ',num2str(id_mode+1)])
            title('POD coeff. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )
            
            id_alp = id_mode+1;
            subplot(122)
            hold on
            % Projected values
            plot( da_PROJCT(:,id_alp) , da_PROJCT(:,id_alp+1) , 'og' )
            % Rom values
            plot( da_ODE45(:,id_alp) , da_ODE45(:,id_alp+1)   ,'-', 'LineWidth', 1 )
            plot( da_ODE15i(:,id_alp) , da_ODE15i(:,id_alp+1)  ,'--r', 'LineWidth', 1 )
            grid on
            xlabel(['POD coeff. ',num2str(id_alp)])
            ylabel(['POD coeff. ',num2str(id_alp+1)])
            title('POD coeff. Time Variation','FontSize', 20,'interpreter','latex' );
            legend( char('Project','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )
            
        end

        % ONE NORMALIZED

%         for id_mode = 1 : ee 
% 
%             figure()
%             subplot(211)
%             hold on
%  
%             plot( da_PROJCT(:,2)./max(da_PROJCT(:,2)) , da_PROJCT(:,id_mode)./max(da_PROJCT(:,id_mode)) , 'or', 'LineWidth', 2)  
%             axis fill
%             plot( da_MONOi(2,end/2:end)./max(da_MONOi(2,end/2:end))  , da_MONOi(id_mode,end/2:end)./max(da_MONOi(id_mode,end/2:end)), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,2)./max(da_ODE15i(end/2:end,2)) , da_ODE15i(end/2:end,id_mode)./max( da_ODE15i(end/2:end,id_mode)) ,'b','LineWidth', 1.1 )
%             plot( da_ODE45(end/2:end,2)./max(da_ODE45(end/2:end,2)) , da_ODE45(end/2:end,id_mode)./max(da_ODE45(end/2:end,id_mode)) , '-.g','LineWidth', 1.1 )
% 
%             grid on
%             xlabel( 'POD coeff. 2' )
%             ylabel(['POD coeff. ',num2str(id_mode)])
%             title('Temp POD coeff. Time Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
%             legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
%             pause(.5)
%             
%             id_alp = id_mode+ee;
%             subplot(212)
%             hold on
% 
%             plot( da_PROJCT(:,2+ee)./max(da_PROJCT(:,2+ee)) , da_PROJCT(:,id_alp)./max( da_PROJCT(:,id_alp)) , 'or', 'LineWidth', 2)  
%             plot( da_MONOi(2+ee,end/2:end)./max(da_MONOi(2+ee,end/2:end))  , da_MONOi(id_alp,end/2:end)./max(da_MONOi(id_alp,end/2:end)), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,2+ee)./max(da_ODE15i(end/2:end,2+ee)) , da_ODE15i(end/2:end,id_alp)./max(da_ODE15i(end/2:end,id_alp)) ,'b','LineWidth', 1.1 )
%             plot( da_ODE45(end/2:end,2+ee)./max(da_ODE45(end/2:end,2+ee)) , da_ODE45(end/2:end,id_alp)./max(da_ODE45(end/2:end,id_alp)) , '-.g','LineWidth', 1.1 )
%            
%             
%             grid on
%             xlabel( 'POD coeff. 2' )
%             ylabel(['POD coeff. ',num2str(id_mode)])
%             title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
%             legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
%             pause( .5 )
%            
%         end
%         
%         
        % ONE 

        for id_mode = 1 : ee 

            figure()
            subplot(211)
            hold on
            
            plot( da_PROJCT(:,2) , da_PROJCT(:,id_mode) , 'or', 'LineWidth', 2)  
            plot( da_MONOi(2,end/2:end)  , da_MONOi(id_mode,end/2:end), 'm','LineWidth', 1.1 )
            plot( da_ODE15i(end/2:end,2) , da_ODE15i(end/2:end,id_mode) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(end/2:end,2) , da_ODE45(end/2:end,id_mode) , '-.g','LineWidth', 1.1 )

            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Temp POD coeff. Time Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause(.5)
            
            id_alp = id_mode+ee;
            subplot(212)
            hold on
            plot( da_PROJCT(:,2+ee) , da_PROJCT(:,id_alp) , 'or', 'LineWidth', 2)  
            plot( da_MONOi(2+ee,end/2:end)  , da_MONOi(id_alp,end/2:end), 'm','LineWidth', 1.1 )
            plot( da_ODE15i(end/2:end,2+ee) , da_ODE15i(end/2:end,id_alp) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(end/2:end,2+ee) , da_ODE45(end/2:end,id_alp) , '-.g','LineWidth', 1.1 )
              
            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )
           
        end
        
        % ONE FROM #Bergmann et all #2009
        figure
                    subplot(221)
        hold on
            plot( da_PROJCT(:,2+ee) , da_PROJCT(:,3+ee) , 'or', 'LineWidth', 2)  
%             plot( da_MONOi(2+ee,end/2:end)  , da_MONOi(3+ee,end/2:end), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,2+ee) , da_ODE15i(end/2:end,3+ee) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(end/2:end,2+ee) , da_ODE45(end/2:end,3+ee) , '-.g','LineWidth', 1.1 )
              
            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )

                    subplot(222)
        hold on
            plot( da_PROJCT(:,2+ee) , da_PROJCT(:,4+ee) , 'or', 'LineWidth', 2)  
%             plot( da_MONOi(2+ee,end/2:end)  , da_MONOi(4+ee,end/2:end), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,2+ee) , da_ODE15i(end/2:end,4+ee) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(end/2:end,2+ee) , da_ODE45(end/2:end,4+ee) , '-.g','LineWidth', 1.1 )
              
            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )
                    subplot(223)
        hold on
            plot( da_PROJCT(:,2+ee) , da_PROJCT(:,5+ee) , 'or', 'LineWidth', 2)  
%             plot( da_MONOi(2+ee,end/2:end)  , da_MONOi(5+ee,end/2:end), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,2+ee) , da_ODE15i(end/2:end,5+ee) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(round(end/2:end),2+ee) , da_ODE45(round(end/2:end),5+ee) , '-.g','LineWidth', 1.1 )
              
            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )

                    subplot(224)
        hold on
            plot( da_PROJCT(:,4+ee) , da_PROJCT(:,5+ee) , 'or', 'LineWidth', 2)  
%             plot( da_MONOi(4+ee,end/2:end)  , da_MONOi(5+ee,end/2:end), 'm','LineWidth', 1.1 )
%             plot( da_ODE15i(end/2:end,4+ee) , da_ODE15i(end/2:end,5+ee) ,'b','LineWidth', 1.1 )
            plot( da_ODE45(round(end/2:end),4+ee) , da_ODE45(round(end/2:end),5+ee) , '-.g','LineWidth', 1.1 )
              
            grid on
            xlabel( 'POD coeff. 2' )
            ylabel(['POD coeff. ',num2str(id_mode)])
            title('Velo POD coeff. Variation - HALF SIMULATION CONSIDERED','FontSize', 20,'interpreter','latex' );
            legend( char('Project','Monoi','ODE45e','ODE15i'),'Location','Best')
            pause( .5 )
            
            
      
     
%% COMPARE KINETIC ENERGYte
        if ee == eeMin

            figure( FC_Kinetic );
            KinEn_FEM = 1/2*diag( SOL.d1(: , IN_QuasiSteady: end  )' * Dnn * SOL.d1(: , IN_QuasiSteady: end ) );
            KinEn_POD = 1/2*diag( Upod_ODE15'   * Dnn *   Upod_ODE15   );
            % % plot( tSteps, KinEn_POD,'-r','LineWidth',2 )
            plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), KinEn_FEM, 'oy','LineWidth',1.2)
            plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), KinEn_POD, ':ob','LineWidth',.5)
            axis tight
            title('Kinetic energy comparison ','FontSize', 20,'interpreter','latex' );
            legend( ['FEM Kinetic. captured by modes';
                     'FEM Kinetic.                  ';
                     'ROM Kinetic.                  ' ] )
            xlabel('time [s]')
            ylabel('Kinetic energy')

            pause( 5 )

        else
            figure( FC_Kinetic );
            KinEn_POD = 1/2*diag( Upod_ODE15'   * Dnn *   Upod_ODE15   );
            plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), KinEn_POD', ':','LineWidth',.5)

            pause( 5 )

        end
    end
    
%% COMPARE THERMAL ENERGY

if ee == eeMin
    
    figure( FC_thermal ); clf ; hold all
    TherEn_FEM = diag( SOL.d1(ta1 , IN_QuasiSteady: end  )' * Dpp(ta1,ta1) * SOL.d1(ta1 , IN_QuasiSteady: end ) );
    TherEn_POD = diag( Upod_ODE15(ta1,:)'   * Dpp(ta1,ta1) *   Upod_ODE15(ta1,:)   );

    plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), TherEn_FEM, ':ob','LineWidth',  2)
    plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), TherEn_POD, ':or','LineWidth', 2)
    axis tight
    title('Thermal energy comparison ','FontSize', 20,'interpreter','latex' );
    xlabel('time [s]')
    ylabel('Kinetic energy')
    LG_thermal_2 = char('FEM Thermal energy',['POD Thermal energy ',num2str(ee),' modes']);
    legend( LG_thermal_2 )
    grid on
    pause( 5 )

else
    figure( FC_thermal ); hold all
    TherEn_POD = diag( Upod_ODE15(ta1,:)'   * Dpp(ta1,ta1) *   Upod_ODE15(ta1,:)   );
    plot( SOL.solinfo.solvals( IN_QuasiSteady: end ), TherEn_POD', '-','LineWidth', 2)
    LG_thermal_2 = char( LG_thermal_2, ['POD Thermal energy ',num2str(ee),' modes']);
    legend( LG_thermal_2 )
    pause( 5 )
    
end



 
    %% RELATIVE ERROR ON THERMAL AND KINETIC ENERGY


    %% PLOT SOLUTION
        if ee == eeMax
    %         close all
    figure()
    abCFD_plot_CompareUinTime(cat(3,SOL.d1(:,IN_QuasiSteady:end),SOL.d1(:,IN_QuasiSteady:end)),cat(3,Upod_ODE15,Upod_ODE15), da_ODE15i(:,alp), ...
                                        model, mesh, var, 2, 1, [3,4] , tSteps )
    % abCFD_plot_CompareUinTime(cat(2,SOL.d1),cat(2,POD), da(:,omg), ...
    %                                     model, mesh, var, 2, 1,2 , tSteps)
        end
    end   

end
        figure(  )
        clf
        bar( ErrBar(:,[5,10,15,25]) )
         grid on
        xlabel(['Probe number'],'FontSize',12)
        ylabel(['|| Uy^{ PRJ}-Uy^{ POD} || _2'],'FontSize',12)
        legend(char('5 modes','10 modes','15 modes','25 modes'))
       
        
         for i_prob = 1 : numel( probes )
       
        figure()
%             StoreErrPlot_Ux{i_prob}(ee,:)=Upod_ODE45(U_prob, : );
%             StoreErrPlot_Uy{i_prob}(ee,:)=Upod_ODE45(V_prob, : );
%             StoreErrPlot_T{i_prob}(ee,:)=Upod_ODE45(T_prob, : );
            tested = [5,15];
            markers = {'.','--' ,'-' };
            LineWidth=[3,2];
            for i_ee = 1 : length(tested)

                subplot(311)
                title(['Probe # ',num2str(i_prob),' for ',num2str(ee),' modes'],'FontSize',15)

                hold all; grid on
                plot(  StoreErrPlot_Ux{i_prob}(tested(i_ee),:),markers{i_ee} , 'LineWidth',LineWidth(i_ee)) 
    %                 plot(  Upod_ODE45(U_prob, : ), '--' , 'LineWidth',2)  
                    ylabel([' Ux [m/s]'],'FontSize',10)
                axis tight
                    subplot(312)
                    hold all; grid on
                plot(  StoreErrPlot_Uy{i_prob}(tested(i_ee),:),markers{i_ee} , 'LineWidth',LineWidth(i_ee)) 
    %                 plot(  Upod_ODE45(V_prob, : ), '--' , 'LineWidth',2) 
                    ylabel([' Uy [m/s]'],'FontSize',10)
                axis tight
                    subplot(313)
                    hold all; grid on
                plot(  StoreErrPlot_T{i_prob}(tested(i_ee),:),markers{i_ee} , 'LineWidth',LineWidth(i_ee))
    %                 plot(  Upod_ODE45(T_prob, : ), '--' , 'LineWidth',2) 
                    ylabel([' T [?]'],'FontSize',10)
                    axis tight
                    xlabel(['Time step'],'FontSize',12)
    %                  legend(char('Full Model','Garlerkin-POD'))
            end
            legend(char('Full Model','ROM 5 modes','ROM 15 modes'))
        end
return % Close function


