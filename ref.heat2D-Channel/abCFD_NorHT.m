% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO CALCULATE THE NORMAL HEAT FLUX  %
%  /----\ |  \|    |--  |   |   TO ALL PARTICULAR EDGES OF 2D RECTANGULAR      %
% /      \|__/ \__ |    |__/    DOMAIN USING COMSOL. INTERPOLATION             %                     
%                               FUNCTIONS ARE USED ONCE THE SOLUTION FIELD AT  %
% A PARTICULAR TIME FRAME IS SENT TO COMSOL BY MODIFYING THE MODEL OBJECT      %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_NorHT( model, ph, dof, edg, ele, solution ) Produced n plots for n
% edges with secondary variable prescribed. 
function [] = abCFD_NorHT( model, ph, dof, edg, ele, solution , A )

if 1 == strcmp( A, 'Y')

    %% INITIALIZE
    % Particular time to compute graT
    timinig = 1;
    Gammas = 1 : 4 ;
    Gammas ( Gammas == ph.BC.Diriclet ) = [];

    %% PLOT HEAT FLUX ON EACH non-DIRICLET BOUNDARY
    % Loop among the edges

    for id_g = 1 : numel(Gammas)
        % Initialize
        gamma = num2str( Gammas( id_g ) );

        switch Gammas( id_g )
            case 1
            % Number of Elements on \gamma    
            [ EleNum , Dofknown ] = size( edg.gamma1.EdgToDof );   
            % Order all the DOF on \ gamma
            Dof = edg.gamma1.EdgToDof ;    
            case 2
            [ EleNum , Dofknown ] = size( edg.gamma2.EdgToDof );   
            Dof = edg.gamma2.EdgToDof ;
            case 3
            % Number of Elements on \gamma    
            [ EleNum , Dofknown ] = size( edg.gamma3.EdgToDof );   
            % Order all the DOF on \ gamma
            Dof = edg.gamma3.EdgToDof ;    
            case 4
            [ EleNum , Dofknown ] = size( edg.gamma4.EdgToDof );   
            Dof = edg.gamma4.EdgToDof ;     
        end

        % Number of DOF per element
        DofNum = size( ele.localCoor , 2);
        % List all elements:
        Ele = ele.DOFToNode' ;
        %% IDENTIFY ELEMENTS' DOFS ON \gamma
        % Loop among the edges of \gamma to identify all the DOF belonging to each
        % element

        % Initilize 
    %     eleTag = zeros( ele.num , Dofknown * DofNum);
        coonnec = zeros( EleNum, DofNum);
        for id_ele = 1 : EleNum 

            eleTag = Ele == Dof( id_ele , 1 );
            for id_dof = 2 : Dofknown
                eleTag = [ eleTag, Ele == Dof( id_ele,id_dof) ];
            end

            so = sum( eleTag, 2);
            [ row, col , val ] = find( so == Dofknown);

            coonnec( id_ele , : ) = ele.DOFToNode( : , row);
        end
        clear eleTag
        %% X,Y,Z COORDINATES PER DOF PER ELEMENT
        XYZ = [];
        for id_ele = 1 : EleNum

            XY = dof.dofCoord( :,  coonnec( id_ele , : ) )';
            tem = [ XY , solution( coonnec( id_ele , : ) , timinig) ];
            XYZ = [ XYZ ; tem];
        end

        %% EXPORT TXT FILE
        % Write in the current fodler a txt file to be sent to Comsol
        dlmwrite(['XYZgamma',gamma,'.txt'] , XYZ, 'delimiter', '\t', 'precision', 6);

        %% MODIFY MODEL OBJECT
        % Deactivate the physics at the first iteration and two plots
        if id_g == 1
            model.physics('ht').active(false);
            model.result('pg1').active(false);
            model.result('pg2').active(false);
        end
        % Activate interpolation functinos
        model.func( [ 'InterpGamma',gamma ] ).active(true);
        model.func( [ 'InterpGamma',gamma ] ).set('filename', [pwd , '\XYZgamma',gamma,'.txt'] );

        % Active interpolation varibles
        model.variable( ['var',gamma] ).active(true);
        % Activat plot and initialize
        PlotGr = model.result('eval_interp_hT');
        PlotGr.name( ['GradGamma_', gamma] );  
        LineUnique = PlotGr.feature( 'lngr1' );
        LineUnique.active(true);
        % Set up for different Gammas( id_g )

        if Gammas( id_g ) == 1 || Gammas( id_g ) == 4
            LineUnique.selection.set( Gammas( id_g ) );
            LineUnique.set( 'expr', ['d( IntGamma',gamma,' ,x)'] );
            LineUnique.set( 'xdataexpr', 'y');

        elseif Gammas( id_g ) == 2 || Gammas( id_g ) == 3
            LineUnique.selection.set( Gammas( id_g ) );
            LineUnique.set( 'expr', ['d( IntGamma',gamma,' ,y)'] );
            LineUnique.set( 'xdataexpr', 'x');
        end 

        % RUN INTERPOLATION
        model.sol('Solver1').runAll;
        % Run plot
        PlotGr.run; 
        % Export plot to Matlab
        figure()
        mphplot(model,'eval_interp_hT')
        % Modify properties
        grid on
        title( ['Heat flux normal to ',gamma] );
    
    end

    % Restore modelobject
    model.physics('ht').active(true);
    model.result('pg1').active(true);
    model.result('pg2').active(true);
    LineUnique.active(false);

end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Rountine tested and completed including plots for heat    26/02/2013 %
%         flux normal to all edges of a particular solution field, eg a        %
%         snapshot or a single mode                                            %
%   0.1 - kickoff. What is missing:                                 25/02/2013 %
%              1. Definition of plots                                          %
%              2. To speed up this routine is better to work with sparse       %
%              matrices in the two loops where working for each element        %
% ---------------------------------------------------------------------------- %