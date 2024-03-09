function [T21,T10,T12] = FEM_Process_T20(T20,segments)
%T20: tableau de connectivit�: D�finition des �l�ments � partir des points
%segments: tableau donnant la d�finition des segments de l'�l�ment en fonction des noeuds (dans la num�rotation locale de l'�l�ment parent)
%T21: tableau de connectivit� donnant la d�finition des �l�ments en termes de segments orient�s.
%T10: tableau de connectivit� des segments en fonction des noeuds



%nombre de segments par �l�ment
Nsegments_el = size(segments,1);
%nombre de noeuds par segment (2=lin�aire, 3=quadratique,...)
Nnodes_seg = size(segments,2);
%nombre d'�l�ments
Nelements = size(T20,1);

%initialisaion des tables
T21 = zeros(Nelements,Nsegments_el);
T10 = zeros(Nelements*Nsegments_el,Nnodes_seg);

%remplissage, pour chaque �l�ment on g�n�re les segments sans se pr�occuper des doublons
idx_seg = 1;
for i = 1:Nelements
    loc_elem = T20(i,:);
   T21(i,:) = idx_seg:(idx_seg+Nsegments_el-1);
   T10(idx_seg:(idx_seg+Nsegments_el-1),:) = loc_elem(segments);
   idx_seg = idx_seg+Nsegments_el;
end

%on inverse les segments pour que le premier noeud soit de plus petit indice
%NB: on inverse toujours en deux blocs: d'abord les deux premiers noeuds puis les suivants (cf. choix de num�rotation locale)
idx_seg = 1;
for el = 1:Nelements
    for seg = 1:Nsegments_el
	%Si n�cessaire on inverse l'ordre des noeuds et on change l'orientation du segment dans la table 21
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

%on supprime les doublons et on g�n�re l'info (N) pour la renum�rotation.
[T10,~,N] = unique(T10,'rows'); 
%renum�rotation des segments de la table 21 apr�s suppression des doublons
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