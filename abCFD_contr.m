% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	MASS AND STIFFNESS CONTRIBUTES FROM A SINGLE   %
%  /----\ |  \|    |--  |   |   ELEMENT DISCRETIZING A SCALAR FIELD            %
% /      \|__/ \__ |    |__/    NO PHYSICAL VALUES ARE INCLUDED                %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_contr( mesh ,var, I )
% INPUTS:
%   - I: structure of the integrand. It's a functions of the global
%       coordinates and is the product of two shape functions (or theirs
%       derivatives)

function [ Ke, De] =abCFD_contr( I , ElType, IntPrecision)

%% COMPUTE GAUSS POINTS AND WEIGHTS
[ xGP , yGP, WGP ] = abCFD_GaussQuad2D( ElType , 2 );
% [Xi,WGP] = hammer_points( 2 );
% xGP = Xi( 1 , :);
% yGP = Xi( 2 , :);

% %% EVALUATE SHAPE FUNCTIONS AND DERIVATIVES
% % Evaluate Shape functions at i_GP-th Gauss point
% LocalShp = cellfun( @(x) polyvaln( x, [xGP',yGP'] ) , I{ 1 }.ShpFun ,...
%                                                     'UniformOutput', false);
% % Evaluate Shape functions derivatives at i_GP-th Gauss point
% LocalDeriv = cellfun( @(x) polyvaln( x, [xGP',yGP'] ) , I{ 1 }.ShpDer ,...
%                                                     'UniformOutput', false);
% %% EXPAND GLOBAL COORDINATE VECTOR
% ic = ones( I{1}.ShpNum ,1);
% for i_ex = 1 : I{1}.ShpNum - 1
%     
%     ic = [ ic ; ic( (end - I{1}.ShpNum + 1):end ,: ) + ones( I{1}.ShpNum,1)]; %#ok<AGROW>
%    
% end
% % IWAN = I{ 1 }.GlbCor( ic , : );
% 
% 
% 
%% EVALUATE INTEGRALS
% %% COMPUTE JACOBIAN
% detJak = det( cell2mat( LocalDeriv)'* IWAN / I{1}.ShpNum ); % It certainly works..but maybe
% %           only for linear mapping with constan jacobian!
% %% COMPUTE COEFF MATRICES
% tic
% De = cell2mat( cellfun( @(x) (detJak/2)*x'*diag((WGP))*...
%     cell2mat(cellfun( @(k) k, LocalShp,'UniformOutput',0)), ...
%     LocalShp,'UniformOutput',0)' ); 
% toc

Ke = zeros( I{1}.ShpNum ,I{1}.ShpNum);
De = zeros( I{1}.ShpNum ,I{1}.ShpNum);
Ce = zeros( I{1}.ShpNum ,I{1}.ShpNum);

for i_GP = 1 : numel( xGP )
    
%     Ke = Ke + detJak * WGP( i_GP ) *
    % Evaluate i_GP-th point with following coords and weight
    xiGP = xGP ( i_GP );
    yiGP = yGP ( i_GP );
    WiGP = WGP ( i_GP );
    
    % Evaluate Shape functions at i_GP-th Gauss point
    LocalShp = cellfun( @(x) polyvaln( x, [xiGP,yiGP] ) , I{ 1 }.ShpFun ,...
                                                    'UniformOutput', false);
    % Evaluate Shape functions derivatives at i_GP-th Gauss point
    LocalDeriv = cellfun( @(x) polyvaln( x, [xiGP,yiGP] ) , I{ 1 }.ShpDer ,...
                                                    'UniformOutput', false);
    
    % Compute Jacobian and its inverse
    % MIND: A smaterter way would be to compute directly its inverse. this
    %       is possible mapping using linear triangles
    Jak = I{ 1 }.GlbCor'*cell2mat( LocalDeriv );
    invJak = inv( Jak );
    GloblDeriv = cell2mat( LocalDeriv ) * invJak ;
    
    % DO NOT Compute strain-displacement commutation matrix
%     B = zeros( 2 , I{1}.ShpNum  );
%     B( 1 , 1 : end ) = GloblDeriv(: , 1 )';
%     B( 2 , 1 : end ) = GloblDeriv(: , 2 )';
%     B( 3 , 1 : I{1}.ShpNum ) = GloblDeriv(: , 2 )';
%     B( 3 , I{1}.ShpNum + 1 : end ) = GloblDeriv(: , 1 )';
    
    % Element Stifness matrix
    Ke = Ke + GloblDeriv * GloblDeriv' * WiGP *det( Jak );
    
    % Element Damping matrix
%     De = De + cell2mat( LocalShp )' * cell2mat( LocalShp ) * WiGP *det( Jak );
    
    % Convective matrix
    De = De + cell2mat( LocalShp )' * cell2mat( LocalShp ) * [ 1 ; 1 ; 1] *...
                GloblDeriv( : , 1)' * WiGP *det( Jak );
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  25/03/2013 %
% ---------------------------------------------------------------------------- %