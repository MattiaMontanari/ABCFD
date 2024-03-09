% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO COMPUTE, ORDER AND PLOT THE     %
%  /----\ |  \|    |--  |   |   EMPIRICAL MODES OF THE REDUCED BASIS.          %
% /      \|__/ \__ |    |__/    PROPER ORGHOGONAL DECOMPOSITION METHOD         %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[eigval,B] = abCFD_ComputeModes( SNAP,model,ee,dof, ele,elements,fun,A)
%   The array SNAP contains solutions fields for time dependent problem, in
%   the array B a number of ee modes is stored. Each one corresponds to a
%   certain rate energy into eigval.

function [ eigval , B ] = abCFD_ComputeModes(SNAP,model,ee,dof, ele,elements,fun,A)

%% Solve eigenvalue problem
switch fun
    case 'eig'
        tic; [ B , eigval ] = eig( SNAP*SNAP' ); toc
    case 'eigs'
        tic; [ B , eigval ] = eigs( SNAP*SNAP' ); toc
end

% Rearrange eigenvalues
eigval = diag(eigval);
EnergyTot = sum( eigval );
eigval = sort( eigval( end - ee + 1 : end ,1) ,'descend') ;
% Compute rate of energy per eigenmode
EnergyEig = 100 / EnergyTot * eigval;
% Reassange eigenmodes
order = ee:-1:1;
B = B(: , end - ee + 1 : end );
B = B( : , order);

%% Plots

% Plot Energy per modes rate
if strcmp( A , 'Y') == 1
    figure()
    bar( EnergyEig )
    title('Energy per modes rate')
end

% Plot each single mode
for id_e = 0 : ee-1
    abCFD_PlotExtMesh(model, A ,dof.dofCoord', ele.DOFToNode', ...
                                            B(:,end-id_e), elements)
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   1.0 - kick-off                                                  14/02/2013 %     
%  
% ---------------------------------------------------------------------------- %