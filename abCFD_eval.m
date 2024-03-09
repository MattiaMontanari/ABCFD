% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    EVALUATES SEPARATELY EACH UNKNOWN FIELD USING  %
%  /----\ |  \|    |--  |   |   'mpheval' WHICH IS FAST BUT CREATES DISORDER   %
% /      \|__/ \__ |    |__/    IN NODE NUMBERING, THUS NOT USEFUL FOR POD.    %
%                               ON THE OTHER HAND, THIS HELPS PLOTTING.        %  
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ U ] = abCFD_eval( model, T_order, V_order, P_order ) evaluates the
% solution fields: u,v,T,p with thieris relative orders and stores all the
% solutions in time in the structure U

function [ U ] = abCFD_eval( model, T_order, V_order, P_order, mesh )

tiin = toc;
%% Extract general information
U.solinfo = mphsolinfo(model,'nu','on');
last = U.solinfo.sizesolvals;
 
%% Extract separately solution fields
% Both velocities
uv = mpheval( model, { 'u','v' } ,'dataset','dset1','refine',V_order,...
                            'smooth','none','solnum', 1:last );
% x-velocity component
U.u = rmfield(uv, 'd2');
U.u.expr = U.u.expr{1};
U.u.unit = U.u.unit{1};
% y-velocity component
U.v = rmfield(uv, 'd1');
oldField = 'd2';
newField = 'd1';
[U.v.(newField)] = U.v.(oldField);
U.v = rmfield(U.v,oldField);
U.v.expr = U.v.expr{2};
U.v.unit = U.v.unit{1};

% Pressure
U.p = mpheval( model, { 'p' } ,'dataset','dset1','refine',P_order,...
                            'smooth','none','solnum', 1:last );
% Temperature
U.T = mpheval( model, { 'T' } ,'dataset','dset1','refine',T_order,...
                            'smooth','none','solnum', 1:last );

fprintf('Solution extraction completed in %f sec \n', (toc-tiin) )

%% ASSEMBLE SOLUTION VECTOR IN TIME

% % Assemble temperature contribue
% % U.d1( mesh{2}.ele{1}.dof{1}.uniqTAG , : ) = U.T.d1';
% % Assemble pressure contribute
% % U.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = U.p.d1';
% % Assemble velocity contribute
% % U.d1( mesh{2}.ele{1}.dof{3}.uniqTAG , : ) = U.u.d1';
% % U.d1( mesh{2}.ele{1}.dof{4}.uniqTAG , : ) = U.v.d1';

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  April 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Assembly of whole solution array in time                  06/04/2013 %
%   1.0 - kick off.                                                 28/02/2013 %
% ---------------------------------------------------------------------------- %