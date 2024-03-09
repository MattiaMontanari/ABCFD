function [Xe,W] = hammer_points(NP)
%position des points d'intégration de hammer
switch NP
    case 1
        Xe = [1;1]/3;
        W = 0.5;
    case 2
        Xe = [0.5 0.5 0 ; 0 0.5 0.5];
        W = [1 1 1]/6;
    case 3
        Xe = [1 0.6 1.8 0.6 ; 1 0.6 0.6 1.8]/6;
        W = [-0.2812500000 0.2604166667 0.2604166667 0.2604166667];
    otherwise
        error('TOTO');
end
end
