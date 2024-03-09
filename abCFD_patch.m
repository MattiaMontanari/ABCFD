% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO PLOT BY MEANS OF THE MATLAB     %
%  /----\ |  \|    |--  |   |   BUILT-IN FUNCTION PATCH. HERE THIS PLOT        %
% /      \|__/ \__ |    |__/    IS SET UP                                      %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%pat =  abCFD_patch( sol, NumSol)
%pat = abCFD_patch( sol, NumSol, d1 )
%pat = abCFD_patch( sol, NumSol, d1, p ,t )

function pat = abCFD_patch( data, NumSol , varargin)

if isfield( data, 't' )
    
    switch nargin 
        case 2

            pat = patch('faces',data.t'+1,'Vertices',[data.p',data.d1( NumSol ,:)'],...
                        'CData',data.d1( NumSol ,:) ,'LineWidth',1,...
                        'EdgeColor','interp','FaceColor','interp');

        case 3
            data.d1 = varargin{1}';


            pat = patch('faces',data.t'+1,'Vertices',[data.p',data.d1( NumSol ,:)'],...
                        'CData',data.d1( NumSol ,:) ,'LineWidth',1,...
                        'EdgeColor','interp','FaceColor','interp');
    end
    
else
    %% Check input
    XY = varargin{1};  % Coordinates
    d1 = data;          % Data
    % reshape coordinate arrays
    [ row, col ] = size( XY );
    if row < col
        XY = XY';
        [ row, col ] = size( XY );
    end
    % reshape data arrays
    if size( d1, 2) ~= row
        d1 = d1';
    end
    
    X = XY( : , 1 );
    Y = XY( : , 2 );
    t = delaunay( X , Y);
        
    pat = patch('faces', t ,'Vertices',[ XY, d1( NumSol ,:)'],...
                        'CData', d1( NumSol ,:) ,'LineWidth',1,...
                        'EdgeColor','interp','FaceColor','interp');
                    
end
axis equal

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  April 2013             % 
% ---                                                                      --- %
%   0.3 - Case 'data' isn't a stucture and input check added        08/04/2013 % 
%   0.2 - This sub-ruotine is called also for mesh plotting.        09/03/2013 % 
%         Added a new output different from the classic FC.                    %
%   1.0 - kick-off                                                  28/02/2013 %     
%  
% ---------------------------------------------------------------------------- %