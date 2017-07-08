function [F_areas] = qareas(m)
% return the areas of the quad mesh elements
F = m.F;
X = m.X
if size(F,2)~=4,error('This function requires a quad mesh');end
%% construct the triangles
persistent C
if size(C,1)~=2*size(F,1),
    disp('Setting up the associated triangulation one time');
    C = zeros(size(F,1) * 2, 3);
    C(1:size(F,1),:) = [F(:,1) F(:,2) F(:,3)];
    C(size(F,1)+1:end,:) = [F(:,1) F(:,3) F(:,4)];
end

u = X(:,1); v = X(:,2); w = X(:,3);
x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2);
twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
A = sum(twoA)/2;                    % this is the total area
F_areas = twoA/2;                   % this is the vector of face areas
F_areas = F_areas(1:size(F,1)) + F_areas(size(F,1)+1:end);