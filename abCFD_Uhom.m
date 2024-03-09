% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	COMPUTE THE HOMOGENEOUS SOLUTION FOR THE       %
%  /----\ |  \|    |--  |   |   FIELDS 'T','u','v' OF A STRCUTRE 'U'. IF ONE   %
% /      \|__/ \__ |    |__/    OF THESE FIELDS DO NOT EXIST, THIS WILL BE NOT %
%                               CONSIDERED.                                    %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%H = abCFD_Uhom( U, R, mesh )

function [ H ] = abCFD_Uhom( U, R, mesh )
 
%% COMPUTE WHOLE SOLUTION ARRAY
H.d1 = U.d1 - R.d1;

%% ASSEMBLE SEPARATELY DEPENDENT VARIABLE FIELDS

% Assemble temparture
H.T.d1 =  H.d1( mesh{2}.ele{1}.dof{1}.uniqTAG , : );
H.T.expr = 'T';
% Assemble pressure
H.p.d1 =  H.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : );
H.p.expr = 'p';
% Assemble velocity
H.u.d1 =  H.d1( mesh{2}.ele{1}.dof{3}.uniqTAG , : );
H.u.expr = 'T';
H.v.d1 =  H.d1( mesh{2}.ele{1}.dof{4}.uniqTAG , : );
H.v.expr = 'T';


%% OBSOLETE
% if isfield( U , 'T' )
% % Temperature
% H.T.d1 = U.T.d1 - R.T.d1 ;
% H.T.p = U.T.p;
% H.T.t = U.T.t;
% H.T.ve = U.T.ve;
% H.T.unit = U.T.unit;
% H.T.expr = U.T.expr;
% end
% 
% if isfield( U , 'u' )
% % velocity u-component
% H.u.d1 = U.u.d1 - R.u.d1 ;
% H.u.p = U.u.p;
% H.u.t = U.u.t;
% H.u.ve = U.u.ve;
% H.u.unit = U.u.unit;
% H.u.expr = U.u.expr;
% end
% 
% % velocity v-component
% if isfield( U , 'v' )
% H.v.d1 = U.v.d1 - R.v.d1 ;
% H.v.p = U.v.p;
% H.v.t = U.v.t;
% H.v.ve = U.v.ve;
% H.v.unit = U.v.unit;
% H.v.expr = U.v.expr;
% end


fprintf('Homogeneous solution computed at %s  \n', datestr(now,15) )



end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   0.2 - Due to numbering issue the whole routine has been changed 06/04/2013 %
%   1.0 - kick-off                                                  28/02/2013 %
% ---------------------------------------------------------------------------- %