function [l1, l2] = principal_stretches(T1, T2)
% returns the principal stretches when deforming triangle T1 to become T2
% Input:
%       T1 and T2 are each 2-vectors with 3 rows, with columns corresponding to
%       x and y coordinates. Each row is a vertex in the triangle in the
%       sequence V1, V2 and V3
% The algorithm is based on SoftMatter article Khairy et al 2012, Supp note 7 and Figure S-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1 = T1(1,:);
V2 = T1(2,:);
V3 = T1(3,:);
v1 = 
Dm12 = 