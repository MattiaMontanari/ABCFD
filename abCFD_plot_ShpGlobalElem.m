% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	PLOTS ALL SHAPE FUNCTIONS OF A SINGLE ELEMENT  %
%  /----\ |  \|    |--  |   |   USING GLOBAL COORDINATES AND MATLAB'S          %
% /      \|__/ \__ |    |__/    INTERPOLATION FUNCTIONS.                       %
%                               WORKS FOR TRIANGUAL 2D ELEMENTS                %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_plot_ShpGlobalElem( mesh ,var, ElTag, MeshTag, VarTag, ElType, FiledInterOrd)
% INPUTS:
%   - ElTag     : Element tag among the whole mesh
%   - MeshTag   : mesh{ MeshTag } contatins the elements 'ElType'.
%                 Typically MeshTag = 2
%   - VarTag    : Choose the variable's tag corresponting to 'T','p','u','v'
%   - ElType    : It cannot be nothing but 'tri'
%   - FiledInterOrd: Grade of shape functions (from 1 to 5)

function abCFD_plot_ShpGlobalElem( mesh ,var, ElTag, MeshTag, VarTag, ElType, FiledInterOrd)

%% INITIALIZE

DofTagGlobal = mesh{ MeshTag }.ele{1}.dof{ VarTag}.TAG( ElTag ,:);
% Global coords
DofCordGlobal = var.dofCoord(:,DofTagGlobal);
% Look for element's verteces
vtx = abCFD_elem_FindVtx ( DofTagGlobal, DofCordGlobal, ElType, FiledInterOrd); 
% Initiazize unitare shape values
ZallNs = eye( numel( DofTagGlobal ) );

% Loop over each single shape function of the element
for id_N = 1 : numel( DofTagGlobal )
    
    % Comput polynomial
    PolyOrdN1 = FiledInterOrd ;
    PolyN1 = polyfitn( DofCordGlobal', ZallNs( : , id_N), PolyOrdN1  );
    
    %% PLOT SLICES
    % Initialize plot
    figure()
    hold on
    
    % The idea is to track three lines along the element's edges, compute
    % the polynomial of these lines ( y = mx + q ) and use these three
    % lines as limitis to the plot
    
    % Slope factors 'm'
    m12 = (vtx.GlobalCoord(2 , 1) - vtx.GlobalCoord(2 , 2)) /...
           (vtx.GlobalCoord(1 , 1) - vtx.GlobalCoord(1 , 2));
    m13 = (vtx.GlobalCoord(2 , 1) - vtx.GlobalCoord(2 , 3)) /...
           (vtx.GlobalCoord(1 , 1) - vtx.GlobalCoord(1 , 3));
    m23 = (vtx.GlobalCoord(2 , 2) - vtx.GlobalCoord(2 , 3)) /...
           (vtx.GlobalCoord(1 , 2) - vtx.GlobalCoord(1 , 3));   
    % Ordinate intercetta 'q'
    q12 = vtx.GlobalCoord(2 , 1) - m12 * vtx.GlobalCoord(1 , 1);
    q13 = vtx.GlobalCoord(2 , 1) - m13 * vtx.GlobalCoord(1 , 1);
    q23 = vtx.GlobalCoord(2 , 2) - m23 * vtx.GlobalCoord(1 , 2);
    
    % Intervals for accurate plot
    intervals = 50;
    % Define slices in at x-constant
    X = linspace( min( DofCordGlobal(1,:) ), max( DofCordGlobal(1,:) ), intervals-1);
    for id_slice = 1 : intervals - 2
       
        % Track #1 of a slice
        Xslice = X( id_slice );
        Xslice1 = ones( intervals , 1 ) .*  Xslice;
        
        Y12 = m12 * Xslice + q12;
        Y(1) = Y12;
        Y13 = m13 * Xslice + q13;
        Y(2) = Y13;
        Y23 = m23 * Xslice + q23;
        Y(3) = Y23;
        
        cond(1) = min(vtx.GlobalCoord(2,:)) <= Y12 && Y12 <= max(vtx.GlobalCoord(2,:));
        cond(2) = min(vtx.GlobalCoord(2,:)) <= Y13 && Y13 <= max(vtx.GlobalCoord(2,:));
        cond(3) = min(vtx.GlobalCoord(2,:)) <= Y23 && Y23 <= max(vtx.GlobalCoord(2,:));
         
        if sum( cond ) == 1
           % you are a one vertex
           Ymin = Y( find( cond == 1 ) ); %#ok<FNDSB>
           Ymax = Ymin;
           
        else 
            Ylimits = Y(find( cond == 1 )); %#ok<FNDSB>
            
        end
        
        Yslice1 = linspace( min( Ylimits), max( Ylimits ), intervals)';
        
        % Track #2 of a slice
        
        Xslice = X( id_slice +1 );
        Xslice2 = ones( intervals , 1 ) .*  Xslice;
        
        Y12 = m12 * Xslice + q12;
        Y(1) = Y12;
        Y13 = m13 * Xslice + q13;
        Y(2) = Y13;
        Y23 = m23 * Xslice + q23;
        Y(3) = Y23;
        
        cond(1) = min(vtx.GlobalCoord(2,:)) <= Y12 && Y12 <= max(vtx.GlobalCoord(2,:));
        cond(2) = min(vtx.GlobalCoord(2,:)) <= Y13 && Y13 <= max(vtx.GlobalCoord(2,:));
        cond(3) = min(vtx.GlobalCoord(2,:)) <= Y23 && Y23 <= max(vtx.GlobalCoord(2,:));
         
        if sum( cond ) == 1
           % you are a one vertex
           Ymin = Y( find( cond == 1 ) ); %#ok<FNDSB>
           Ymax = Ymin;
           
        else 
            Ylimits = Y(find( cond == 1 )); %#ok<FNDSB>
            
        end
        
        Yslice2 = linspace( min( Ylimits), max( Ylimits ), intervals)';
        
        
        % Plot
        XO = [Xslice1;Xslice2] ; YO = [ Yslice1 ; Yslice2];
        tri= delaunay( XO , YO );
        zg = polyvaln( PolyN1 , [XO(:),YO(:)] );
        patch('faces',tri,'Vertices',[XO,YO,zg],...
                    'CData', zg ,'LineWidth',1,...
                    'EdgeColor','interp','FaceColor','interp');
                
    end
        
        
end
     
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   0.1 - kick-off                                                  22/03/2013 %
% ---------------------------------------------------------------------------- %