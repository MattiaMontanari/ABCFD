% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	MAIN ROUTINE FOR ASSEMBLING COEFFICIENT FEM    %
%  /----\ |  \|    |--  |   |   MATRICES FOR EACH:                             %
% /      \|__/ \__ |    |__/     - MESH TYPE (e.g. 'tri' or 'quad')            %
%                                - VARIABLE  (e.g. 'mod_T','mod_u','mod_v'...) %
%                                - MESH ELEMENT                                %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_AssembleMain( mesh ,var, Smeshingo )
% INPUTS are smiply common strcutures
% 
% NOTES:
%       - running on Adrien's code for scalar field (temperature)


function [ Kvar, Dvar, B, K,D] = abCFD_AssembleMain( mesh ,var, SmeshInfo )

% Initialize loop over the mesh type
NumMesh = numel( SmeshInfo.stat.types);
% Total number of dofs
TotDof = sum( var.fieldNDofs );
% Initialize coefficient matrices
K = sparse( TotDof , TotDof);
D = sparse( TotDof , TotDof);
B = zeros( TotDof , 1);
% Loop over the mesh types
for i_mesh = 1 : NumMesh
    % Switch for different mesh type
    MeshType = char( SmeshInfo.stat.types( i_mesh ) );
    switch MeshType( 1 , 1:3)
        case 'tri'
            
            % Initialize loop over the variables
            NumVars = numel( var.DofsNames );
%             NumVars = 1; %  - - - - - - - WRONG, I JUST WANNA DEAL WITH T - - - -
            
            
            % Loop over the variables
            for i_var = 1 : NumVars
                
                    % Current variable to be treated
                CurrentVar = char( var.DofsNames( i_var ) );
                    % Find local variables for i_var
                LocalVar = strfind( mesh{ i_mesh }.localDofNames,  CurrentVar );
                LocalVar = find( cellfun( @(x) isempty(x) , LocalVar) == 0)';
                    % Local Coordinates
                LocalXY = mesh{ i_mesh }.localDofCoords( : , LocalVar );
                    % dofs' tags
                GlbConnectivity = mesh{ i_mesh }.ele{ 1 }.dof{ i_var }.TAG ;
                
                % Switch for different variables (scalar o vector fields)
                switch CurrentVar( 1 , end )      
                    case 'T' % It's a scalar field
                        
                        %% ADRIEN'S CODE
%     [ KADR , DADR ,Bvar ] = exemple_2D_TRI( var.dofCoord'  , GlbConnectivity );
                        
                        %% MINE
                        % Define integration order
                        IntOrd = 2;
                        % Compute contributes                        
 [ Kvar , Dvar, Bvar] = abCFD_AssembleScalar(var.dofCoord', GlbConnectivity, LocalXY, IntOrd, MeshType,TotDof);
                        
                    case 'u' || 'v'
                        Kvar = sparse( TotDof , TotDof);
                        Dvar = sparse( TotDof , TotDof);
                        Bvar= sparse( TotDof , TotDof);
                        
                    case 'p'
                        Kvar = sparse( TotDof , TotDof);
                        Dvar = sparse( TotDof , TotDof);
                        Bvar= sparse( TotDof , 1);
                end
                
                K = K + Kvar;
                D = D + Dvar;
                B = B + Bvar;
                
            end

            
            
        case 'qua'
            
        case 'edg'
            
        case 'vtx'
            
    end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  26/03/2013 %
% ---------------------------------------------------------------------------- %