% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	PLOTS MATRICES' PATTERN IN SUBPLOTS            %
%  /----\ |  \|    |--  |   |                                                  %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_plot_MatrixPattern( assemb, var, mesh, A ) 
% where:
%       - 'var','mesh','assemb' are the common structures
%       - A is a string. If A ='KDC' the in the new figure you have:
%          i.   Stiffness matrix ( reordered & deleted)
%          ii.  Damping matrix ( reordered & deleted )
%          iii. Both matrices in consecutive ordering of the unknowns
%               (eg. [ T ; U ; V ; P]

function abCFD_plot_MatrixPattern( assemb, var, mesh, A )
% Initialize subplot counter
count = 0 ;
if 1 == strcmp( A(1) , 'K') 
%% STIFFNESS MATRIX
    count = 1 + count;
    figure();
    subplot( 2,1, count );     hold on
    % Plot deleted stiffness matrix
    spy(assemb.Kc, 5 ,'og');
 
    % Plot reordered stiffness matrix
    spy(assemb.K, 3 ,'xr')
    title(' STIFFNESS  MATRICES SPARSITY PATTERN')
    legend( [ 'Removed   (',abCFD_IsSymm( assemb.Kc),')'] ,...
            [ 'Reordered (',abCFD_IsSymm( assemb.K), ')'] );    
%             [ 'Assembled (',abCFD_IsSymm( del    ),  ')'] ,...
            
end

%% DAMPING MATRIX
if 1 == strcmp( A(2) , 'D') 
    count = 1 + count;
    subplot( 2,1, count );     hold on
    % Plot deleted stiffness matrix
    spy(assemb.Dc,5,'og');
     
    % Plot reordered stiffness matrix
    spy(assemb.D, 2 ,'.r')
    title(' DAMPING MATRICES SPARSITY PATTERN')
    legend( [ 'Removed   (',abCFD_IsSymm( assemb.Dc),')'] ,...
            [ 'Reordered (',abCFD_IsSymm( assemb.D), ')'] );    
end
%% PLOT COEFFCIENT MATRICES WITH CONSECUTIVE REORDERING
%   Consecutive numbering on unknowns [ T , p , u , v ]
if 1 == strcmp( A(3), 'C')
    
    % Extract all fields
    FiledS = char( var.DofsNames);
    % This has been not generalized yet, so suppose a particular variables config
    if strcmp( FiledS, ['mod1_T';'mod1_p';'mod1_u';'mod1_v'] ) == 1
        
        VarDofs =[ unique(mesh{2}.ele{1}.dof{1}.TAG);...
                   unique(mesh{2}.ele{1}.dof{2}.TAG);...
                   unique(mesh{2}.ele{1}.dof{3}.TAG);...
                   unique(mesh{2}.ele{1}.dof{4}.TAG)];
    elseif strcmp( FiledS, ['mod1_p';'mod1_u';'mod1_v'] ) == 1

        VarDofs =[ unique(mesh{2}.ele{1}.dof{1}.TAG);...
                   unique(mesh{2}.ele{1}.dof{2}.TAG);...
                   unique(mesh{2}.ele{1}.dof{3}.TAG)];
    else
        VarDofs = []; % Works for any general variables config but without a
                      %   predefined arrangement of the matrices
        for i_var = 1 : size( FiledS , 1 )

            VarDofs = [VarDofs; unique( mesh{2}.ele{1}.dof{ i_var }.TAG ) ]; %#ok<AGROW>
        end
    end

    % PLOT MATRICES

    % Stiffness Matrix
    figure()
    subplot( 2,1, 1)
    spy( assemb.K(VarDofs,VarDofs) , 1.5 , '*b' ); hold on
    title('Stiffness Matrix')
     
    % Damping Matrix
    subplot( 2,1, 2)
    spy( assemb.D(VarDofs,VarDofs) , 1.5 , '*b' ); hold on
    title('Damping Matrix')

    % subplot( x,x, x)
    % spy( assemb.N(VarDofs,VarDofs) , 1.5 , '*y' ); hold on
    % spy( assemb.NF(VarDofs,VarDofs) , 1.5 , '*r' ); hold on
    % title('constraint ')
    % legend( 'constraint Jacobian' , 'constraint force Jacobian Matrix' )
end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off.                                                 20/03/2013 %     
% ---------------------------------------------------------------------------- %