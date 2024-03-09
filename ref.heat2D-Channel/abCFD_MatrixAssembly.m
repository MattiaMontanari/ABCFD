% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE ASSEMBLE THE SPARSE MATRICES FROM  %
%  /----\ |  \|    |--  |   |   COMSOL.                                        %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ assemb , FC_spy ] = abCFD_MatrixAssembly(assemb, model, dof, A ) enriches the 
% structure variable 'assemb' with the sparse coefficient matrices as new fileds


function [ assemb , FC_spy ] = abCFD_MatrixAssembly(assemb, model, dof, A )
%% Assembly one by one matrix
for id_m = 1 : size( assemb.matrices , 2 )
    v = model.sol('Solver1').feature('asmbl').getSparseMatrixVal( assemb.matrices{ id_m } );
    c = model.sol('Solver1').feature('asmbl').getSparseMatrixCol( assemb.matrices{ id_m } );
    r = model.sol('Solver1').feature('asmbl').getSparseMatrixRow( assemb.matrices{ id_m } );
    assembly = sparse( double( r )+1 ,double( c )+1 , v );
    assemb = setfield(assemb, assemb.matrices{ id_m } , assembly );
end

for id_v = 1 : size( assemb.vector , 2 )
    vect = model.sol('Solver1').feature('asmbl').getVector( assemb.vector{ id_v } );
    assemb = setfield(assemb, assemb.vector{ id_v } , vect );
end

if strcmp( A , 'Y')
    figure()
    FC_spy = gcf;
    spy(assemb.Kc,2,'og');
    hold on
    del = assemb.K;
    del = del( dof.DofToNode , dof.DofToNode);
    spy( del , 1.5 , '*b' )
    clear del
    spy(assemb.K,2,'xr')
    title('Coefficient Matrices Sparsity Pattern')
    legend('Removed','Assembled','Reordered')
else
    FC_spy = [];
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   1.1 - All the vector and matrices are extracted directly here   18/02/2013 %
%   1.0 - Added third matrix concerning the NODES' numbering        14/02/2013 %
%   0.2 - OPTIMIZED THE OUTPUT IF 'A'~='Y'                          13/02/2013 %
%   0.1 - kick-off. Only the matrices are extracted and assembled   11/02/2013 %     
%         the vectors are not treated but listed in assemb.vector
% ---------------------------------------------------------------------------- %