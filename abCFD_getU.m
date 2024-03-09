% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    EXTRACTS THE UNKNOWN VECTOR COMPUTED BY COMSOL %
%  /----\ |  \|    |--  |   |   WHILE SOLVING THE LINEAR SYSTEM.               %
% /      \|__/ \__ |    |__/    NO BUILT-IN FUNCTIONS ARE USED, IT IS SLOW     % 
%                               SINCE EXTRACTS THE SOLUTION ARRAIS             % 
% ONE-BY-ONE, BUT IT KEEPS THE CORRECT DOF' NUMBERING AND THUS USEFUL FOR POD  % 
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ U ] = abCFD_ 

function [ U ] = abCFD_getU( model, Solver, mesh )

tiin = toc;
%% Extract general information
U.solinfo = mphsolinfo(model,'nu','on');
last = U.solinfo.sizesolvals;
 
%% EXTRACT UNKNOWN ARRAY

for id_t = 1 : last
   U.d1( : , id_t ) = model.sol( Solver ).getU( id_t ,'Sol',1);
end

%% ASSEMBLE SEPARATELY DEPENDENT VARIABLE FIELDS

% Assemble temparture
U.T.d1 =  U.d1( mesh{2}.ele{1}.dof{1}.uniqTAG , : );
U.T.expr = 'T';
% Assemble pressure
U.p.d1 =  U.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : );
U.p.expr = 'p';
% Assemble velocity
U.u.d1 =  U.d1( mesh{2}.ele{1}.dof{3}.uniqTAG , : );
U.u.expr = 'T';
U.v.d1 =  U.d1( mesh{2}.ele{1}.dof{4}.uniqTAG , : );
U.v.expr = 'T';

fprintf('Solution extraction completed in %f sec \n', (toc-tiin) )



end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  April 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off.                                                 06/04/2013 %
% ---------------------------------------------------------------------------- %