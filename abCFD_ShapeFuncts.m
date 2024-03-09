% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	COMPUTES ISOPARAMETRIC ELEMENTS' SHAPE         %
%  /----\ |  \|    |--  |   |   FUNCTIONS AND THERIS DERIVATIVES WITH RESPECT  %
% /      \|__/ \__ |    |__/    TO THE NATURAL(LOCAL) COORDINATES. THE OUTPUT  %
%                               IS A STRUCTURE WHICH DEFINES A POLYNOMIAL      %
% FOR CERTAIN SHAPE FUNCTIONS A POLYNOMIAL INTERPOLATION IS USED               %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ ShpFun , ShpDer] = abCFD_ShapeFuncts( LocalXY, MeshType )
% INPUT:
%   - LocalXY local coordinates of N points stored in [2,N] array.
%
% LocalXY = mesh{2}.localDofCoords( : , 1:6);

function [ ShpFun , ShpDer] = abCFD_ShapeFuncts( LocalXY, MeshType)

% INITIALIZE

% Number shape functions ( nodes ) & Number independent variables(x,y,z)
[ VarNum , ShpNum ] = size( LocalXY );
% 

%% SHAPE FUNCTIONS 
if strcmp( MeshType , 'tri' )  && ShpNum == 3
    switch ShpNum
        case 3
            ShpFun{ 1 }.ModelTerms = [  1   0 ; ...
                                        0   1 ; ...
                                        0   0];
            ShpFun{ 1 }.Coefficients = [    -1  -1  1   ];
            ShpFun{ 1 }.VarNames = {'X1'  'X2'};
            
            ShpFun{ 2 }.ModelTerms = [  1   0 ; ...
                                        0   1 ; ...
                                        0   0];
            ShpFun{ 2 }.Coefficients = [    1   0   0   ];
            ShpFun{ 2 }.VarNames = {'X1'  'X2'};

            ShpFun{ 3 }.ModelTerms = [  1   0 ; ...
                                        0   1 ; ...
                                        0   0];
            ShpFun{ 3 }.Coefficients = [    0   1   0   ];
            ShpFun{ 3 }.VarNames = {'X1'  'X2'};

        case 6
    end
    
elseif strcmp( MeshType , 'quad' ) 
    switch ShpNum
        case 4
            
     	case 9
    end
else
    % CASES NOT INCLUDED BEFORE
    % Initialize Store shape functions is vector cell
    ShpFun = cell( ShpNum  , 1 );
    % Initiaze unitare shape values
     ZallNs = eye( ShpNum );

    % Loop over shape functions
    for i_shp = 1 : ShpNum
        % Interpolate polynomials
        ShpFun{ i_shp }  = polyfitn( LocalXY', ZallNs( : , i_shp ), 1 );
    end
end
%% SHAPE FUNCTIONS DERIVATIVES

% if strcmp( MeshType , 'tri' ) 
%     switch ShpNum
%         case 3
%             
%      	case 6
%     end
%     
% elseif strcmp( MeshType , 'quad' ) 
%     
% else
            
    ShpDer = cell( ShpNum , VarNum );
    % Loop over shape functions
    for i_shp = 1 : ShpNum
        ShpDer( i_shp , : )  =  abCFD_PolynDer( ShpFun{ i_shp } ) ;
    end

% end
end


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  22/03/2013 %
% ---------------------------------------------------------------------------- %