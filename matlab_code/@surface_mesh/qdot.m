function [D] = qdot(m)
% returns the four dot products. This can be used as a measure of quad
% element deformation
X = m.X;
F = m.F;
if size(F,2)~=4,error('This function requires a quad mesh');end
persistent A
if size(A,1)~= size(F,1)
    A = zeros(size(F,1), 4);
end
v1 = X(F(:,1),:)-X(F(:,2),:);v2 = X(F(:,1),:)-X(F(:,4),:);A(:,1) = dot(v1,v2,2);
v1 = X(F(:,2),:)-X(F(:,1),:);v2 = X(F(:,2),:)-X(F(:,3),:);A(:,2) = dot(v1,v2,2);
v1 = X(F(:,3),:)-X(F(:,2),:);v2 = X(F(:,3),:)-X(F(:,4),:);A(:,3) = dot(v1,v2,2);
v1 = X(F(:,4),:)-X(F(:,3),:);v2 = X(F(:,4),:)-X(F(:,1),:);A(:,4) = dot(v1,v2,2);
D = A;