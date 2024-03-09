% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	DEFINES GAUSS POINTS AND WEIGHTS FOR 2D GAUSS  %
%  /----\ |  \|    |--  |   |   QUADRATURE RULE                                %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ xGP , yGP, WGP ] = abCFD_GaussQuad2D( ElType , IntPrecision )
% INPUTS:
%   - ElType: Geometry of element: 'tri' or 'quad'
%   - IntPrecision: integration precision
%           IF == 2 uses 3 points on element's 2D domain
%           IF == 5 uses 7 points on element's 2D domain
% OUTPUT:
%   - xGP:  Gauss points' coord in x-local coord
%   - yGP:  Gauss points' coord in y-local coord
%   - WGP:  Weights associated to Gauss points
%    % ref. BATHE p.467

function [ xGP , yGP, WGP ] = abCFD_GaussQuad2D( ElType , IntPrecision )


switch ElType 
    case 'tri'
        switch IntPrecision
            case 2
        % Natuar coordinates 
        xGP = [ 2/3, 1/6, 1/6 ];
        yGP = [ 1/6, 2/3, 1/6 ];
        % Weights
        WGP = [1 1 1]/6;
            
            case 5
        % Natural coordinates
        xGP = [ 0.1012865073235,...
                0.7974269853531,...
                0.1012865073235,...
                0.4701420641051,...
                0.4701420641051,...
                0.0597158717898,...
                0.3333333333333];
        yGP = [     xGP( 1 )  ,...
                    xGP( 1 )  ,...
                    xGP( 2 )  ,...
                    xGP( 6 )  ,...
                    xGP( 4 )  ,...
                    xGP( 4 )  ,...
                    xGP( 7 )  ];
                
         % Weights
         WGP = [ .1259391505448,...
                 .1259391505448,...
                 .1259391505448,...
                 .1323941527885,...
                 .1323941527885,...
                 .1323941527885,...
                 .225           ]/2;
        end
    case 'quad'
         
	  error('Check Gauss point integration specification')
        % Refer to BATHE p.466
end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  22/03/2013 %
% ---------------------------------------------------------------------------- %