% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	ASSEMBLE SPARCE ARRAYS WITH CONTRIBUTES FROM   %
%  /----\ |  \|    |--  |   |   ANY BOUNDARY CONDITION.                        %
% /      \|__/ \__ |    |__/                                                   %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %

% DofNum = U.solinfo.size

function [ assemb ] = abCFD_AssembleAllBC( mesh , assemb, DofNum, var )

% Find whic cell containts info on 'edg' elements
indexT = 0;
it = 1;
while indexT == 0 && it <= numel( mesh )
    indexT = strcmp(mesh{ it }.Type,'edg');
    indexC = it;    % Cell containing 'vtdx'
    it = it + 1;
end

%% INITIALIZE
EdgNum = numel( mesh{ indexC }.ele );
DofsNum = numel( var.DofsNames );

%% Loop over the variables
unit = cell( 1 , DofsNum);
for i_var = 1 : DofsNum
    
    % Local dofs' coordinates on the edge 
    LocalDofName = char(mesh{ indexC }.localDofNames);
    LocalDofName = strfind( char(LocalDofName(:,6))', var.DofsNames{i_var}( 1,6) );
    LocalCoord = mesh{indexC}.localDofCoords( LocalDofName );
    
    % Initialize
    toto = cell( 1 , EdgNum);
    
    for i_edg = 1 : EdgNum

        % Evaluate Contribute on Domain       
        Double = abCFD_AssembleEdgCntrb( 2 , ...
            mesh{indexC}.ele{ i_edg }.dof{i_var}.TAG, DofNum, var.dofCoord );
    
        % Evaluate Robin's contribute on Stiffness matrix
        Single = abCFD_AssembleEdgCntrb( 1 , ...
            mesh{indexC}.ele{ i_edg }.dof{i_var}.TAG, DofNum, var.dofCoord );
        
        toto{ i_edg } = struct( 'single' , Single, 'double', Double );
        
    end
%      toto = struct( 'single' , Single, 'double', Double);
    % Store the data for i_var into the structure 'assemb'
    unit{i_var}= toto;

    
end

assemb.unit = cell2struct( unit', var.DofsNames,  1);

end



% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off. Included null output in case of non-specified   20/02/2013 %
%         secondary boundary conditions and fields in 'assemb'.
%                                                                              %
% ---------------------------------------------------------------------------- %