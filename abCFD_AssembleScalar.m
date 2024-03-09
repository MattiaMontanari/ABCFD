% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	COMPUTE MASS AND STIFFNESS MATRICES FOR THE    %
%  /----\ |  \|    |--  |   |   TEMPERATURE (SCALAR) FIELD                     %
% /      \|__/ \__ |    |__/    NO PHYSICAL VALUES ARE INCLUDED                %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
% 
% 
function [Kvar , Dvar, Bvar] = abCFD_AssembleScalar( ...
            dofCoord  , GlbConnectivity, LocalXY, IntOrd, MeshType, TotDof )

%% INITIALIZE

% Number of elements in the mesh
NumEle = max( size( GlbConnectivity ) );
% Damping matrix
Dvar = sparse( TotDof , TotDof);
% Mass matrix
Kvar = sparse( TotDof , TotDof);
% Forcing vector
Bvar = sparse( TotDof ,    1  );

% Shape functions and theirs derivatives 
[ ShpFun , ShpDer] = abCFD_ShapeFuncts( LocalXY, MeshType );



for i_ele = 1 : NumEle
    % Global coordinates
    I{ 1 }.GlbCor = dofCoord( GlbConnectivity( i_ele ,:) , :);
    I{ 1 }.GlbCnct = GlbConnectivity( i_ele ,:);
    I{ 1 }.ShpFun = ShpFun; % Polynomial of shape functions
    I{ 1 }.ShpDer = ShpDer; % Polynomial of shape functions
    I{ 1 }.ShpNum = numel(ShpFun);
    [ Ke, De ] = abCFD_contr( I , MeshType, IntOrd);

    % Mapping 
    Kvar( I{ 1 }.GlbCnct , I{ 1 }.GlbCnct) = Kvar( I{ 1 }.GlbCnct , I{ 1 }.GlbCnct) + Ke;
    Dvar( I{ 1 }.GlbCnct , I{ 1 }.GlbCnct) = Dvar( I{ 1 }.GlbCnct , I{ 1 }.GlbCnct) + De;
end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  25/03/2013 %
% ---------------------------------------------------------------------------- %