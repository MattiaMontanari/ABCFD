% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO EXTRAC DATA OF EXTENDED MESH    %
%  /----\ |  \|    |--  |   |   AND GEOMETRY. TWO COMMANDS ARE USED:           %
% /      \|__/ \__ |    |__/    XmeshInfo and mphmeshstats                     %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[gen,dof,node,elem] = abCFD_GeoXmeshInfo( model ) provides four structure 
%   variables about EDGES, DOF, VERTEX and ELEMENTS

function [ SmeshInfo, dof, ele , edg ] = abCFD_GeoXmeshInfo( model )
%% DESCRIPTION OF OUTPUT STRUCTURES
% edg.
%  - .GeoToMesh
%  - .EdgToNode
%  - .EdgToDofs

%% Initialize (xmiStudy and xmiVarV1 should be equal in this case)
mesh = model.mesh('mesh1');
solver = model.sol('Solver1');
xmiStudy = solver.feature('st1').xmeshInfo();
xmiVarV1 = solver.feature('v1').xmeshInfo(); 

%%   WORK OUT THE STRUCTURE:  dof

% dof.number        -> number of dofs
dof.number_V1 = xmiVarV1.fieldNDofs();     
dof.number_V1 = double( dof.number_V1 );
% dof.name          -> variable name
dof.name_V1 = xmiVarV1.fieldNames; 
dof.name_V1 = char(dof.name_V1);
% dof.types         -> Types of elements into the mesh
dof.types_V1 = xmiVarV1.meshTypes;
dof.types_V1 = cell ( dof.types_V1 );          
% dof.NodeToDof     -> from NODES to DOF's 
dof.NodeToDof_V1 = xmiVarV1.dofs().nodes() + 1;
dof.NodeToDof_V1 = double( dof.NodeToDof_V1 );
% dof.coord         -> Coordinates of DOF's
dof.dofCoord =  xmiVarV1.dofs().coords(); 
% dof.nodeCoord     -> Coordinates of DOF's on NODE's reference
dof.nodeCoord =  xmiVarV1.nodes().coords();
% dof.DofToNode     -> From DOFs to NODEs
dof.DofToNode = xmiVarV1.nodes().dofs()+1; 
dof.DofToNode = double( dof.DofToNode );
% dof.GeoToMesh     -> tag's geo-verteces onto ordered mesh verteces
dof.GeoToNode = mesh.getElemEntity( dof.types_V1{3} );
dof.GeoToNode = double( dof.GeoToNode );
% dof.DofToVtx      -> ordered mesh elements onto tags of mesh nodes
dof.DofToVtx = mesh.getElem( dof.types_V1{3} ) + 1;
dof.DofToVtx = double( dof.DofToVtx );

%%   WORK OUT THE STRUCTURE:  edg

% edg.GeoToMesh     -> tag's geo-elements onto ordered mesh edges
edg.GeoToMesh = mesh.getElemEntity( dof.types_V1{1} ); 
edg.GeoToMesh = double( edg.GeoToMesh );
% edg.EdgToDof      -> Xinfo for each DOF on edges (mind internal edges!)
edg.EdgToDof = xmiStudy().elements( dof.types_V1{1} ).dofs+1;
edg.EdgToDof = double( edg.EdgToDof);
% edg.EdgToNode     -> Xinfo for each NODE on edges (mind internal edges!)
edg.EdgToNode = xmiStudy().elements( dof.types_V1{1} ).nodes+1;
edg.EdgToNode = double( edg.EdgToNode);
% edg.localCoor     -> Local coordinates on edge
edg.localCoor = xmiStudy().elements( dof.types_V1{1} ).localCoords;
% edg.gamma1.EdgToDof   -> Xinfo on edges on \Gamma1
edg.gamma1.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==1)';
edg.gamma1.coordDof = dof.dofCoord(:,edg.gamma1.EdgToDof)';
edg.gamma1.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==1)';
% edg.gamma2.EdgToDof   -> Xinfo on edges on \Gamma2
edg.gamma2.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==2)';
edg.gamma2.coordDof = dof.dofCoord(:,edg.gamma2.EdgToDof)';
edg.gamma2.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==2)';
% edg.gamma3.EdgToDof   -> Xinfo on edges on \Gamma3
edg.gamma3.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==3)';
edg.gamma3.coordDof = dof.dofCoord(:,edg.gamma3.EdgToDof)';
edg.gamma3.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==3)';
% edg.gamma4.EdgToDof   -> Xinfo on edges on \Gamma4
edg.gamma4.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==4)';
edg.gamma4.coordDof = dof.dofCoord(:,edg.gamma4.EdgToDof)';
edg.gamma4.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==4)';
% edg.gamma5.EdgToDof   -> Xinfo on edges on \Gamma5
edg.gamma5.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==5)';
edg.gamma5.coordDof = dof.dofCoord(:,edg.gamma5.EdgToDof)';
edg.gamma5.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==5)';
% edg.gamma6.EdgToDof   -> Xinfo on edges on \Gamma6
edg.gamma6.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==6)';
edg.gamma6.coordDof = dof.dofCoord(:,edg.gamma6.EdgToDof)';
edg.gamma6.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==6)';
% edg.gamma7.EdgToDof   -> Xinfo on edges on \Gamma7
edg.gamma7.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==7)';
edg.gamma7.coordDof = dof.dofCoord(:,edg.gamma7.EdgToDof)';
edg.gamma7.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==7)';
% edg.gamma8.EdgToDof   -> Xinfo on edges on \Gamma8
edg.gamma8.EdgToDof = edg.EdgToDof(:,edg.GeoToMesh==8)';
edg.gamma8.coordDof = dof.dofCoord(:,edg.gamma8.EdgToDof)';
edg.gamma8.EdgToNode = edg.EdgToNode(:,edg.GeoToMesh==8)';

%%   EXTEND THE STRUCTURE:  'dof'  FOR THE EDGES
% dof.gamma1.DofTags    -> Xinfo on DOF laying on \gamma1
dof.gamma1.DofTags = unique(edg.gamma1.EdgToDof');
% dof.gamma2.DofTags    -> Xinfo on DOF laying on \gamma2
dof.gamma2.DofTags = unique(edg.gamma2.EdgToDof');
% dof.gamma3.DofTags    -> Xinfo on DOF laying on \gamma3
dof.gamma3.DofTags = unique(edg.gamma3.EdgToDof');
% dof.gamma4.DofTags    -> Xinfo on DOF laying on \gamma4
dof.gamma4.DofTags = unique(edg.gamma4.EdgToDof');
% dof.gamma5.DofTags    -> Xinfo on DOF laying on \gamma5
dof.gamma5.DofTags = unique(edg.gamma5.EdgToDof');
% dof.gamma6.DofTags    -> Xinfo on DOF laying on \gamma6
dof.gamma6.DofTags = unique(edg.gamma6.EdgToDof');
% dof.gamma7.DofTags    -> Xinfo on DOF laying on \gamma7
dof.gamma7.DofTags = unique(edg.gamma7.EdgToDof');
% dof.gamma8.DofTags    -> Xinfo on DOF laying on \gamma8
dof.gamma8.DofTags = unique(edg.gamma8.EdgToDof');

%%   WORK OUT THE STRUCTURE:  vtx


%%   WORK OUT THE STRUCTURE:  ele
% XINFO for EACH ELEMENT ( type_el_mesh{2} )OF THE EXTENDED
ele.NodeToDOF = xmiStudy().elements( dof.types_V1{2} ).nodes+1;
ele.NodeToDOF = double( ele.NodeToDOF );
ele.localCoor = xmiStudy().elements( dof.types_V1{2} ).localCoords();
ele.localDOF = xmiStudy().elements( dof.types_V1{2} ).localDofNames;
ele.localDOF = cell( ele.localDOF);
ele.DOFToNode = xmiStudy().elements( dof.types_V1{2} ).dofs+1;
ele.DOFToNode = double( ele.DOFToNode);

% ele.GeoToMesh     -> tag's geo-domains onto ordered mesh elments
ele.GeoToNode = mesh.getElemEntity( dof.types_V1{2} );
ele.GeoToNode = double( ele.GeoToNode );

% Add element number
ele.num = size(ele.NodeToDOF ,2 );

% CLEAN MEMORY
model.sol('Solver1').feature('st1').clearXmesh();
model.sol('Solver1').feature('v1').clearXmesh();

%           ======================== 
%% -------  USE MPHMESHSTATS COMMAND  ------------------------- ref. LLMatLab
%           ========================

[ mesh_stat , mesh_data ] = mphmeshstats( model, 'mesh1');
% The two output strucures included:
%  << mesh_stat >> among the others fields:
%    [ ]    - TYPES: cell array with type names
%    [ ]    - NUMELEM: vector with number of elemetns of each type
%  << mesh_data >> all it contains is:
%    [ ]    - VERTEX: coordinates of all the nodes into the mesh
%    [ ]    - ELEM: cell array with connectivity matrix for each element type
%                   gives same results xmiStudy().elements( type_el_mesh{x} ) 
%    [ ]    - ELEMENTITY: 'entity information' of each element type - i.e.
%           correlation between mesh's elements and geometry entities
%                   gives same results as GeoToMesh.EDG
 SmeshInfo = struct('stat',mesh_stat,'data',mesh_data);





%% OBSOLETE COMMANDS OR FOR THE SIMPLE MESH (Smesh) THESE MAY BE INCLUDED INTO ANOTHER .m FILE!

%          ===================== 
% -------  USE XMESHINFO COMMAND  ----------------------------- ref. API p.484
%          =====================
% dof.GeoNum_V1 = xmiVarV1.dofs().geomNums();  % dof belowging to each geo's tag ---% Multiphysics purpuse
% dof.gCoord_V1 = xmiVarV1.dofs().gCoords();      % Similiar but in geom. unit  
% dof.idNames_V1 = xmiStudy.dofs().nameInds + 1 ; % Match DOFs' tags and names ---% Multiphysics purpuse 
% node.NODE_Coord_V1 = xmiVarV1.nodes().coords();      % Coordinates of nodes' -----% DELETE#  

end


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Added file ele.num to count mesh's elements               25/02/2013 %
%   0.2 - Corrected a bug the section # EXTEND THE STRUCTURE DOF    21/02/2013 %
%         for which I had problem while a geo edg had only a single element 
%   1.0 - Big change in outputs structures, the overall output      14/02/2013 %     
%         strategy has changed and reduced to few well organized variables.
%         A large section 'Obsolete command' is added as comment and
%         includes working commands which provided a redoundant results, or
%         simply results not need so far. Details should be added to
%         describe the output structures and the structures 'vtx' and 'ele'
%         are not completed yet.
%   0.1 - kick-off                                                  06/02/2013 %     
%  
% ---------------------------------------------------------------------------- %