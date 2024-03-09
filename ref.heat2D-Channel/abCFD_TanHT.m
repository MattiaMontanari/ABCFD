% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO CALCULATE THE TANGENTIAL HEAT   %
%  /----\ |  \|    |--  |   |   FLUX TO A PARTICULAR EDGE. IT ALSO PLOTS THE   %
% /      \|__/ \__ |    |__/    FLUX ALONG THE REAL COORDINATES OF THE EDGE    %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[HT_tan] = abCFD_TanHT( U, ph, dof, edg, gamma,A ) Computes and plots the 
% tangential heat flux 'hT_tan' to the edg \gamma of the solution field 
% solution. It works only for LINEAR ELEMENTS!

% EXAMPLE : abCFD_TanHT( U.fem, ph1, dof, edg, edg.gamma4 )
function [] = abCFD_TanHT( solution , ph, dof, edg, gamma , A )

if 2 == numel(edg.localCoor) && 1 == strcmp( A, 'Y')

    % Compute the flux at time:
    timing = ph.nt;
    % Extarct correct time-solution
    T = solution(:,timing);
    % num ele on \gamma
    NumEle = size( gamma.EdgToDof,1);

    Beb = [ 0 1 ];
    gradT = 0;
    figure(); hold all
    for id_el = 1 : NumEle
    %   DOFs'ä tags of the considered element id_el
        DofTag = gamma.EdgToDof( id_el, :);
    %     Temperture at DofTag
        Tel = T( DofTag , 1);
    %     Y-Cordinates of DofTag
        C = [ ones( size(edg.localCoor,2) , 1) , dof.dofCoord(2,DofTag)' ];
        gradT( id_el , : ) = Beb * inv( C ) * Tel ;

        plot( dof.dofCoord(2,DofTag), gradT( id_el , : ) * [1 ,1],'*')
    end
    
else
    return

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kickoff. Implementation for linear elements only.         25/02/2013 %
% ---------------------------------------------------------------------------- %