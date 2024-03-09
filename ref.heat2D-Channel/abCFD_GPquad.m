% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE PERFORMING NUMERICAL INTEGRATION,  %
%  /----\ |  \|    |--  |   |   WITH GAUSS QUADRATURE RULE, OF SHAPE FUNCTIONS %
% /      \|__/ \__ |    |__/    (AND THEIRS DERIVATIVES)OVER A SINGLE ELEMENT  %                     
%                               ELEMENT AND GIVES THE ASSEMBLED LOCAL MATRIX   % 
% WITH LOCAL CONNECTIVITY MATRIX READY TO BE MAPPED INTO THE GLOBAL  MATRICES  %
% THE INTEGRATION MAY INVOLVE A SINGLE SHAPE FUNCTION OF THE PRODUCT OF TWO.   %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ Q ] = abCFD_GPquad( ShpNum, a, b,ele, ShapInt, A )
%where: ShpNum  - is the number of the sphape functions to be integrated
%       a,b     - integral's interval [a,b] on a straight line
%       ele     - is the structure variable concering all the elements
%       ShapInt - defines if evaluating the integral of a single test
%               function(ShapInt = 1), or the product of two(ShapInt = 2)
%       A       - Answer for plotting which is = 'n' by default not not
%               user-controlled!

function [ Q ] = abCFD_GPquad( ShpNum, a, b,ele, ShapInt, A )
%% INITIALIZE VALUES
% Number shape functions per 1D element
ord = ShpNum - 1;
% Initialize X-coordinates of nodes
Xcoord = linspace( a, b, ShpNum);
% Value of the each shape function on each node
ShpVal = eye( ShpNum );
%% COMPUTE GAUSS POINTS AND WEIGHTS
[ xi , wi ] = lgwt( ShpNum   , a , b );
%% LOOP OVER THE SHAPE FUNCTIONS
% Initialize loop
Ycoord = zeros( ShpNum) ;     % Function evaluated at points xi
p = zeros( ShpNum );          % Polynomial coeff. vector

for id_shp = 1 : ShpNum
    % Get coefficients of the polynomial in the form:
 % Ycoord = p(1)*Xcoord^ord + p(2)*Xcoord^(ord-1) +...+ p(N)*Xcoord + p(N+1)
    p( id_shp , : ) = polyfit( Xcoord, ShpVal( id_shp , : ) , ord);  
    Ycoord( id_shp , : )  = polyval( p( id_shp , : ) , xi);
    
end

switch ShapInt
    case 1
        Q = Ycoord * diag(wi) * ones(ShpNum,1);
    case 2
        Q = Ycoord * diag(wi) * Ycoord';
end


% plot( Xcoord , idShpeFunction  , 'bo')
    
if strcmp( A, 'Y' ) == 1
     abCFD_PlotElement( 'Y' , ele.localCoor ); hold all
     view( 49 , 29 )
     title('Legendre-Gauss Quadrature rule applied to')

    for id_shp = 1 : ord + 1
        hold all
        plot3( linspace(0,1) , zeros(100,1),polyval( p( id_shp , : ) , linspace(a,b)) ...
            , 'LineWidth',2 )

    end
    grid on
    axis tight
end


%%  - - - - - - - - -  NESTED FUNCTION   - - - - - - - - - - - - - - - - - - - %

function [x,w]=lgwt(N,a,b)

%     SOURCE:
% http://www.mathworks.com/matlabcentral/fileexchange/
%                           4540-legendre-gauss-quadrature-weights-and-nodes
% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
% 
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice,
%  this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
%  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
%  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
    N=N-1;
    N1=N+1; N2=N+2;

    xu=linspace(-1,1,N1)';

    % Initial guess
    y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

    % Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);

    % Derivative of LGVM
    Lp=zeros(N1,N2);

    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method

    y0=2;

    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0))>eps


        L(:,1)=1;
        Lp(:,1)=0;

        L(:,2)=y;
        Lp(:,2)=1;

        for k=2:N1
            L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end

        Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   

        y0=y;
        y=y0-L(:,N2)./Lp;

    end

    % Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;      

    % Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

    end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Support for integral of both a SINGLE shape function and  20/02/2013 %
%         the product of TWO shap functions. Implemented also the nested       %
%         function by Greg von Winckel.                                        %
%   0.1 - kick-off                                                  19/02/2013 %     
% ---------------------------------------------------------------------------- %