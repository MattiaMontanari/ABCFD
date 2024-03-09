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

function [ eigval , B, EnergyEig ] = abCFD_EigProblem( sol , ee, fun )

%% INITIALIZE
if isstruct( sol )
    SNAP = sol.d1';
else
    SNAP = sol';
end

%% Solve eigenvalue problem
switch fun
    case 'eig'
        [ B , eigval ] = eig( SNAP'*SNAP ); 
    case 'eigs'
        [ B , eigval ] = eigs( SNAP'*SNAP );
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

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  April 2013             % 
% ---                                                                      --- %
%   1.0 - Modify to accept as input both: structures or arrays      06/04/2013 %     
%   1.0 - kick-off                                                  28/02/2013 %
% ---------------------------------------------------------------------------- %