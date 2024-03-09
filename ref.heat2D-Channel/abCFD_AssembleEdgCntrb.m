% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO ASSEMBLE A SPARCE ARRAY, IN     %
%  /----\ |  \|    |--  |   |   GLOBAL DOF NUMBERING, CONTAINING THE INTEGRALS %
% /      \|__/ \__ |    |__/    OF ONE TEST FUNCTION OR THE PRODUCT OF TWO     %                     
%                               BASIS FUNCTIONS. THE INTEGRAL IS EVALUATED     %
% ALONG A PARTICULAR EDGE OF THE GEOMETRY. PHYSICAL QUANTITIES ARE NOT         %
% INCLUDED IN THIS INEGRALS!                                                   %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ Contribut ] = abCFD_AssembleEdgCntrb( ShapInt,gamma ,dof, ele,A )
%where: ShapInt - defines if evaluating the integral of a single test
%               function(ShapInt = 1), or the product of two(ShapInt = 2)
%       gamma   - is the structure variable concering a particuar edge
%       dof     - is the structure variable concering all the dofs
%       ele     - is the structure variable concering all the elements
%       A       - Answer to some plot request
% TEST: assembWithRobin.K - assembNoRobin.K - Contribut =?= ZERO
%       used to validate the Version 0.2 of this routine

function [ Contribut ] = abCFD_AssembleEdgCntrb(ShapInt,gamma ,dof, ele,A)

%% Initialize
% Obtain number of elements and connectivity for elements on \Gamma
[NumEle , NumLocalDof] = size( gamma{1}) ;

%% Loop over each element
% Initialized output
Contribut = sparse( dof.number_V1 , 1+( (ShapInt-1)*(dof.number_V1-1) ) );
for id_ele = 1 : NumEle
    % Map row, i.e. test functions
    map_i = gamma{1}( id_ele, : );
    % Map columns, i.e. shape functions
    map_j = gamma{1}( id_ele, : );
    % Coordinates of nodes of id_ele-th element on \Gamma
    coord =  dof.dofCoord( :, gamma{1}( id_ele, :)  );
    
    % Verify coordinates
    dmax = pdist( [coord(:,1) , coord(:,NumLocalDof)]' );
    for id_int = 1 : NumLocalDof-2
        d1 = pdist( [coord(:, id_int+1 ) , coord(:, 1 )]' );
        d2 = pdist( [coord(:, id_int+1 ) , coord(:,NumLocalDof)]' );
        if d1 > dmax || d2 > dmax
             error('Distance not computed correctly');
        end
    end
    
    % Define integral boundaries
    a = 0;
    b = dmax;
    % Compute local coeff. matrix
    [ Q ] = abCFD_GPquad( NumLocalDof , a, b, ele, ShapInt,A);
    % Map local coeff. matrix
    switch ShapInt
        case 1
            Contribut( map_i , 1 ) = Contribut( map_i , 1 ) + Q;
        case 2
            Contribut( map_i , map_j ) = Contribut( map_i , map_j ) + Q;
    end
end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Support for integral of both a SINGLE shape function and  20/02/2013 %
%         the product of TWO shap functions. Control made through the
%         variable 'ShapInt'. The variable 'Contribut' is initialized with
%         a formula which gives: 1 if ShapInt=1 or dof.number_V1 if ShapInt=2
%         Few changes have been made in order to deal with cells instead of
%         structures concerning 'gamma'
%   0.2 - Routine tested for elements:'M2','T1','T2','T5','Q5','Q1' 19/02/2013 %
%         TEST: assembWithRobin.K - assembNoRobin.K - Contribut =(check)= 0    %
%         Further versions should include:                                     %
%             1. General formulation to calculate and surface integral of      %
%                both shape functions and theirs derivatives. Consider also    %
%                the case in which a function = 1 to assembly the load         %
%                vector. (see Adrien's code)                                   %
%             2. Round values before exporting                                 %
%   0.1 - kick-off                                                  18/02/2013 %     
%  
% ---------------------------------------------------------------------------- %