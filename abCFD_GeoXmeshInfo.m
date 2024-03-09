% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO EXTRAC DATA OF EXTENDED MESH    %
%  /----\ |  \|    |--  |   |   AND GEOMETRY. TWO COMMANDS ARE USED:           %
% /      \|__/ \__ |    |__/    XmeshInfo and mphmeshstats                     %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ SmeshInfo, var, mesh ] = abCFD_GeoXmeshInfo( model ) 

% Gives three outputs structures:
% DESCRIPTION OF OUTPUT STRUCTURES
% * * SmeshInfo has two sub-strucures * * * *
%   << stat >> among the others fields:
%    [ ] - types: cell array with mesh element names
%    [ ] - numelem: vector with number of element of each type
%    [ ] - qualitydistr: vector with elements' qualiti information
%   << data >> all it contains is:
%     [ ]    - vertex: coordinates of all mesh's vertexes
%     [ ]    - elem: cell array with connectivity matrix for SIMPLE MESH 
%     [ ]    - elementity: correlation between mesh's elements and geometry 
% * * var has 11 sub-strucures/cells * * * * *
%   <<  fieldNames >> {cell}
%   <<  fieldNDofs >> [int32]
%   <<   meshTypes >> {cell}
%   <<   DofsNames >> {cell}
%   <<    nameInds >> [int32]
%   <<    dofCoord >> [double]
%   <<   NodeToDof >> [int32]
%   <<     geomNum >> [int32]
%   <<  NodesNames >> {cell}
%   << NodesToDofS >> [int32]
%   <<   nodeCoord >> [double]
% * * mesh has three sub-strucures * * * * * * *

function [ SmeshInfo, var, mesh ] = abCFD_GeoXmeshInfo( model )

tiin = toc;

%% Initialize 

SolvTag = cell( model.sol.tags );
UsefulSolTag = 1;
SolvTag = SolvTag{ UsefulSolTag };
MeshTag = char( model.mesh.tags );
%(CHECK IF xmiStudy and xmiVarV1 are equal in this case)
grid = model.mesh( MeshTag );
solver = model.sol( SolvTag );
% xmiStudy = solver.feature('st1').xmeshInfo();
xmiVarV1 = solver.feature('v1').xmeshInfo(); 

%% EXTRACT ALL DATA FROM xmiVarV1

% All variables names
var.fieldNames = cell( xmiVarV1.fieldNames);
% All number dofs
var.fieldNDofs = xmiVarV1.fieldNDofs;
% All element types
var.meshTypes = cell( xmiVarV1.meshTypes );

%  EXTRACT DATA ONLY FROM SUBSTRUCTURE xmiVarV1.dofs 
% DOF names
var.DofsNames = cell( xmiVarV1.dofs.dofNames );
% Dof's index with respect to the DOF's name
var.nameInds = xmiVarV1.dofs.nameInds + 1 ;
% All DOF's COORD
var.dofCoord = double( xmiVarV1.dofs().coords() );
% All Nodes to Dof
var.NodeToDof = xmiVarV1.dofs().nodes() + 1;
% Geomety corrispondance
var.geomNum = xmiVarV1.dofs.geomNums;

%  EXTRACT DATA ONLY FROM SUBSTRUCTURE  xmiVarV1.nodes 
% NODE names
var.NodesNames = cell( xmiVarV1.nodes.dofNames );
% node to dof
var.NodesToDofS = xmiVarV1.nodes.dofs + 1;
% All nodes' coord
var.nodeCoord = double( xmiVarV1.nodes().coords() );
  
%% CONSTRUCT VAR OUTPUT STRUCTURE
% Loop over single variable
for id_dof = 1 :numel( var.DofsNames )
    % Make the names DOFS' and NODEs' names usable in matlab
    var.DofsNames{id_dof}(1,5) = '_';
    var.NodesNames{id_dof}(1,5) = '_';
end

%% CONSTRUCT MESH OUTPUT STRUCTURE
% initialize structures
mesh = cell( numel( var.meshTypes ) , 1) ;

% PARTICULAR - MESH TYPE % Loop over mesh types
for t = 1 : numel( var.meshTypes ) % eg: i = 1 : 3 ['edg',tri',vtx']
    mesh{ t }.Type = var.meshTypes{ t };
    % tag's geo-verteces onto ordered mesh verteces
    mesh{ t }.getElemEntity = grid.getElemEntity( mesh{ t }.Type );
    mesh{ t }.getElemEntity = double( mesh{ t }.getElemEntity );
    % ordered mesh elements onto tags of mesh nodes
    mesh{ t }.getElem = grid.getElem( mesh{ t }.Type ) + 1  ;
    mesh{ t }.getElem = double( mesh{ t }.getElem );
    % Xinfo for each DOF on edges (mind internal edges!)
    mesh{ t }.elementsDOFS = xmiVarV1().elements( mesh{ t }.Type ).dofs + 1;
    mesh{ t }.elementsDOFS = double( mesh{ t }.elementsDOFS );
    % Xinfo for each DOF on edges  
    mesh{ t }.elementsNODE = xmiVarV1().elements( mesh{ t }.Type ).nodes + 1;
    mesh{ t }.elementsNODE = double( mesh{ t }.elementsNODE );
    % Local coordinates on edge
    mesh{ t }.localCoor = xmiVarV1().elements( mesh{ t }.Type ).localCoords;
    % Local dof's names
    mesh{ t }.localDofNames = ...
        cell( xmiVarV1().elements( mesh{ t }.Type ).localDofNames );   
    for id_dof = 1 :numel( mesh{ t }.localDofNames )
        mesh{ t }.localDofNames{id_dof}(1,5) = '_';
    end
    % Local coordinates on edge
    mesh{ t }.localDofCoords = ...
        xmiVarV1().elements( mesh{ t }.Type ).localDofCoords;
        % GEOMETRY's PARTICULAR BOUNDARY of PARTICULAR MESH TYPE       
    for e = 1 : max( mesh{ t }.getElemEntity)         
        % Xinfo on edges on \Gamma
        mesh{ t }.ele{e}.elementsDOFS = ...
            mesh{ t }.elementsDOFS( : ,mesh{ t }.getElemEntity == e )';
%         mesh{ t }.ele{e}.dofCoord = ...
%             var.dofCoord( : , mesh{ t }.ele{e}.elementsDOFS)'; 
        mesh{ t }.ele{e}.elementsNODE = ...
            mesh{ t }.elementsNODE( : ,mesh{ t }.getElemEntity == e )';
        
        % All dofs on e-th element of t-type
        AllDof = mesh{ t }.elementsDOFS(:, mesh{ t }.getElemEntity == e);
        
        % PARTICULAR VARIABLE on PARTICULAR BOUNDAR of PARTICULAR MESH TYPE
        for v = 1 : numel( var.DofsNames )
            % Variable name
            mesh{ t }.ele{e}.dof{ v }.name = var.DofsNames(v,:);
            % DOF's tag of v-variable (T,p,u,v) on e-th element of t-type
                % Initialize
                mesh{ t }.ele{e}.dof{ v }.TAG = [];
                for i_vet = 1 : size( AllDof, 2) % Loop over the elements
                    vet_tag = var.nameInds( AllDof(:,i_vet), 1) == v ;
                    PartDof = AllDof( vet_tag,i_vet ) ;
                    mesh{ t }.ele{e}.dof{ v }.TAG = ...
                                [ mesh{ t }.ele{e}.dof{ v }.TAG ;   PartDof'     ];
                end
                mesh{ t }.ele{e}.dof{ v }.uniqTAG = ...
           unique( reshape(mesh{ t }.ele{e}.dof{ v }.TAG', numel( mesh{ t }.ele{e}.dof{ v }.TAG) ,[]) ,'stable');
            % Numver of degrees of freedom
            mesh{ t }.ele{e}.dof{ v }.numb = size( unique(mesh{ t }.ele{e}.dof{ v }.TAG ),1) ;
            % Types of elements into the mesh
            mesh{ t }.ele{e}.dof{ v }.types = var.meshTypes{ t } ;             
        end
    end
end
         
fprintf('Xmesh data extracted in %f minutes \n', (toc-tiin)/60 )


%% CLEAN MEMORY

model.sol( SolvTag ).feature('st1').clearXmesh();
model.sol( SolvTag ).feature('v1').clearXmesh();

%           ======================== 
%% -------  USE MPHMESHSTATS COMMAND  ------------------------- ref. LLMatLab
%           ========================

[ mesh_stat , mesh_data ] = mphmeshstats( model,  MeshTag );


% Find if among the elements you have 'quad'. If yes invert two rows
indexC = strfind( mesh_stat.types , 'quad');
index = find(not(cellfun('isempty', indexC)));

if numel( index ) ~= 0
    mesh_data.elem{ index }( [4,3] , :) = mesh_data.elem{ index }(  [3,4],:);
end
    
SmeshInfo = struct('stat',mesh_stat,'data',mesh_data);

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.7                                 date:  APRIL 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.7 - deleted field .elementsDOFS                               15/04/2013 %
%   0.6 - Properly restablished the field uniqTag                   06/04/2013 %
%   0.5 - Restablished the field uniqTag                            05/04/2013 %
%   0.4 - Overall review                                            27/03/2013 %
%   0.3 - Written for general mesh and solver                       20/03/2013 %
%   0.2 - Added proper row invertion in case of quadrilateral elem. 10/03/2013 %
%         Added also few loops in order to change Comsol's standard
%         notation for the dofs' names.
%   0.1 - kick-off. Partially described the outputs                 28/02/2013 %
% ---------------------------------------------------------------------------- %