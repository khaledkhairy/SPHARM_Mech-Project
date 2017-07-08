function [m_in, m_out] = get_inner_outer_surfaces(obj, d, flag)
% returns the inners surface removed from the obj surface by d/2, where d
% is the shell thickness
m_in = obj;
m_out = obj;
if sign(flag) == -1,
    [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props([obj.X(:,1) -obj.X(:,2) obj.X(:,3)], obj.F, 0);
else
    [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(obj.X, obj.F, 0);
end

m_in.X = obj.X + N .* d/2;  %%% now generate the inner surface
m_out.X = obj.X - N .* d/2;  %%% now generate the outer surface