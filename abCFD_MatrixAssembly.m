% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	ASSEMBLES THE SPARSE MATRICES AND VECTORS BY   %
%  /----\ |  \|    |--  |   |   COMSOL.                                        %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[assemb,FC_spy]=abCFD_MatrixAssembly(model,assemb, UsefulSolTag)
% Enriches the structure 'assemb' with the assembled matrices and vectors.
% Inputs:
%       - model : Comsol's obj
%       - assemb: assembly structure
%       - SolTag: if == 1 transient solver. (Solver #1)
%                 if == 2 stationary solver. (Solver #2)

function [ assemb ] = abCFD_MatrixAssembly( model, assemb, SolTag  )

tiin = toc;
%% Initialize 
PhyTag = cell( model.physics.tags );
SinPhasFlow = PhyTag{1};
% HeatTransfe = PhyTag{2};
SolvTag = cell( model.sol.tags );
SolvTag = SolvTag{ SolTag };

%% Assembly matrices one by one 
for id_m = 1 : size( assemb.matrices , 2 )
    v = model.sol( SolvTag ).feature('asmbl').getSparseMatrixVal( assemb.matrices{ id_m } );
    c = model.sol( SolvTag ).feature('asmbl').getSparseMatrixCol( assemb.matrices{ id_m } );
    r = model.sol( SolvTag ).feature('asmbl').getSparseMatrixRow( assemb.matrices{ id_m } );
    assembly = sparse( double( r )+1 ,double( c )+1 , v );
    assemb = setfield(assemb, assemb.matrices{ id_m } , assembly );
end
%% ASSEMBLE VECTORS
for id_v = 1 : size( assemb.vector , 2 )
    vect = model.sol( SolvTag ).feature('asmbl').getVector( assemb.vector{ id_v } );
    assemb = setfield(assemb, assemb.vector{ id_v } , vect );
end

 
fprintf('Matrices and vector assembled in %f sec \n', (toc-tiin) )

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.7                                 date:  APRIL 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.7 - Adapted to use solvers with differnt tags                 05/04/2013 %
%   0.6 - Excluded plot and added case for: steady study and UNI    21/03/2013 %
%   0.5 - General review                                            08/03/2013 %
%   0.4 - All the vector and matrices are extracted directly here   18/02/2013 %
%   0.3 - Added third matrix concerning the NODES' numbering        14/02/2013 %
%   0.2 - OPTIMIZED THE OUTPUT IF 'A'~='Y'                          13/02/2013 %
%   0.1 - kick-off. Only the matrices are extracted and assembled   11/02/2013 %     
%         the vectors are not treated but listed in assemb.vector
% ---------------------------------------------------------------------------- %