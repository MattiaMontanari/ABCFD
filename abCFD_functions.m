% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SET COMSOL's FUNCTIONS FOR A    %
%  /----\ |  \|    |--  |   |   GIVEN PHASE. SETS THREE FUNCTIONS, EACH OF A   %
% /      \|__/ \__ |    |__/    SINGLE BOUNDARY CONDITIONS TYPE, PLUS OTHER    %
%                               FOUR INTERPOLATION FUNCTIONS, TO ASSESS THE    %
% VARIABLES' GRADIENT AT THE BOUNDARIES.                                       %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_functions( model )

function abCFD_functions( model, ModelName  )

%% RAMP FUNCTIONS FOR BOUNDARY CONDITIONS

% Initialize 
% physics(1) = ph.ht;
physicsname(1,:) = 'HT';
% physics(2) = ph.ht;
physicsname(2,:) = 'fd';

for id_p = 1 : 2
    name_Robin = [  physicsname(id_p,:) , '_RobinRamp' ];
    name_Diriclet = [ physicsname(id_p,:) , '_DiricletRamp'] ;
    name_Neumann = [ physicsname(id_p,:) , '_NeumannRamp'];

    DiricletRamp = model.func.create( name_Diriclet , 'Ramp');
    DiricletRamp.name( name_Diriclet );
    DiricletRamp.set('funcname', name_Diriclet );
    DiricletRamp.set('location', [physicsname(id_p,:),'_DIR_','t_initial'] );             
    DiricletRamp.set('cutoff', [physicsname(id_p,:),'_DIR_','cutoff'] );
    DiricletRamp.set('slope', [physicsname(id_p,:),'_DIR_','slope'] );

    RobinRamp = model.func.create( name_Robin , 'Ramp');
    RobinRamp.name( name_Robin );                 % Name of the function shown
    RobinRamp.set('funcname', name_Robin);    
    RobinRamp.set('location', [physicsname(id_p,:),'_ROB_','t_initial'] ) ;  
    RobinRamp.set('cutoff', [physicsname(id_p,:),'_ROB_','cutoff'] ) ;
    RobinRamp.set('slope', [physicsname(id_p,:),'_ROB_','slope'] );

    NeumannRamp = model.func.create( name_Neumann , 'Ramp');
    NeumannRamp.name( name_Neumann );                       
    NeumannRamp.set('funcname', name_Neumann );                  
    NeumannRamp.set('location', [physicsname(id_p,:),'_NEU_','t_initial'] ) ;           
    NeumannRamp.set('cutoff', [physicsname(id_p,:),'_NEU_','cutoff'] ) ;
    NeumannRamp.set('slope',[physicsname(id_p,:),'_NEU_','slope'] );
end

%% INTERPOLATION FUNCTIONS FOR BOUNDARY FLUX EVALUATIONS

% On \Gamma1
IntFun = model.func.create('InterpGamma1', 'Interpolation');
IntFun.model( ModelName );
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma1', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');   
model.func('InterpGamma1').active(false);

% On \Gamma2
IntFun = model.func.create('InterpGamma2', 'Interpolation');
IntFun.model( ModelName );
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma2', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma2').active(false);

% On \Gamma3
IntFun = model.func.create('InterpGamma3', 'Interpolation');
IntFun.model( ModelName );
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma3', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma3').active(false);

% On \Gamma4
IntFun = model.func.create('InterpGamma4', 'Interpolation');
IntFun.model( ModelName );
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma4', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma4').active(false);

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  27/02/2013 %     
% ---------------------------------------------------------------------------- %