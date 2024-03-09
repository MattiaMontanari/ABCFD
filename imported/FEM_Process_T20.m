function [T21,T10,T12] = FEM_Process_T20(T20,segments)
%T20: tableau de connectivité: Définition des éléments à partir des points
%segments: tableau donnant la définition des segments de l'élément en fonction des noeuds (dans la numérotation locale de l'élément parent)
%T21: tableau de connectivité donnant la définition des éléments en termes de segments orientés.
%T10: tableau de connectivité des segments en fonction des noeuds



%nombre de segments par élément
Nsegments_el = size(segments,1);
%nombre de noeuds par segment (2=linéaire, 3=quadratique,...)
Nnodes_seg = size(segments,2);
%nombre d'éléments
Nelements = size(T20,1);

%initialisaion des tables
T21 = zeros(Nelements,Nsegments_el);
T10 = zeros(Nelements*Nsegments_el,Nnodes_seg);

%remplissage, pour chaque élément on génère les segments sans se préoccuper des doublons
idx_seg = 1;
for i = 1:Nelements
    loc_elem = T20(i,:);
   T21(i,:) = idx_seg:(idx_seg+Nsegments_el-1);
   T10(idx_seg:(idx_seg+Nsegments_el-1),:) = loc_elem(segments);
   idx_seg = idx_seg+Nsegments_el;
end

%on inverse les segments pour que le premier noeud soit de plus petit indice
%NB: on inverse toujours en deux blocs: d'abord les deux premiers noeuds puis les suivants (cf. choix de numérotation locale)
idx_seg = 1;
for el = 1:Nelements
    for seg = 1:Nsegments_el
	%Si nécessaire on inverse l'ordre des noeuds et on change l'orientation du segment dans la table 21
        if (T10(idx_seg,1) > T10(idx_seg,2))
            if (size(T10,2)>2)
            T10(idx_seg,:) = [fliplr(T10(idx_seg,1:2)) fliplr(T10(idx_seg,3:end))];
            else
                T10(idx_seg,:) = fliplr(T10(idx_seg,1:2));
            end
            T21(el,seg) = -T21(el,seg);
        end
        idx_seg = idx_seg+1;
    end
end

%on supprime les doublons et on génère l'info (N) pour la renumérotation.
[T10,~,N] = unique(T10,'rows'); 
%renumérotation des segments de la table 21 après suppression des doublons
for el = 1:Nelements
    for seg = 1:Nsegments_el
        sens = sign(T21(el,seg));
        myseg = abs(T21(el,seg));
        T21(el,seg) = sens*N(myseg);
        idx_seg = idx_seg+1;
    end
end


T12 = zeros(size(T10,1),2);
for el = 1:Nelements
    for seg = 1:Nsegments_el
        sens = sign(T21(el,seg));
        pos = 1+(1-sens)/2;
        T12(abs(T21(el,seg)),pos) = sens*el;
    end
end
end