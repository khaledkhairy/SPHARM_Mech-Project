function [t,p] = newton_steps_02(t,p,F,F_areas,newton_step, maxiter,verbose)
% Minimizes the areas of triangles (on the sphere) by calculating new
% positions t and p for the vertices.
% The algorithm uses Newton's method iteratively

F_A = double(F_areas./sum(F_areas(:)));
%%%%% Prepare initial quantities    
nvert = length(t);
p = mod(p,2*pi);
X = [t(:);p(:)];
u = zeros(length(t),1,'double');
v = zeros(length(t),1,'double');
w = zeros(length(t),1,'double');


%% prepare calculation of Jacobian sparsity pattern
JacPat = sparse(length(F),nvert);
for ix = 1:length(F),verts = F(ix,:);for vert = 1:length(verts),JacPat(ix,verts(vert)) = 1;end;end
JacPat = sparse([JacPat JacPat]);   % we have two of the same pattern
%indJ = find(JacPat);
[i,j] = find(JacPat);
indJ = sub2ind(size(JacPat),i,j);       % linear indices of non-zero entries into Jacobian
JacInd = JacPat;JacInd(indJ) = indJ;    % fill in the jacobian pattern with the indices
pos_vec = zeros(size(indJ));        % gives the position of the variable relative to the face definition
indcol = zeros(length(indJ), 1, 'uint64');
indrow = zeros(length(indJ), 1, 'uint64');

%%% prepare quantities for efficient area_plus calculation
disp('Initializing Jacobian');
indVertVal = zeros(length(indJ),6); % 6 columns corrspond to theta and phi values for the three verices of each triangle
for ix = 1:length(indJ),        % loop over the Jacobian Pattern non-zeros
    [row,col] = ind2sub(size(JacPat),indJ(ix));
    indvec = JacInd(row,:);
    indvec = indvec(indvec>0);
    pos_vec(ix) = find(indvec==indJ(ix));
    indVertVal(ix,:) = JacInd(indvec);
    indcol(ix)   = col;       %%
    indrow(ix)   = row;
end
disp('Done!');
Cv = zeros(size(indVertVal));   % the actual coordinates of the faces are stored in here
indCv = sub2ind(size(Cv),1:length(indJ),pos_vec');
J = sparse(size(JacPat,1), size(JacPat,2));
seps = sqrt(eps);



%% %%%%%%%%%%%%%%% Begin iterations %%%%%%%%%%%%%%%%%%%%%
%  we are taking Newton steps as follows:
%  assuming that our system of equations (constraints) is Areas(S), we need to find 
%  dS which makes Areas(S+dS) = 0;
%  According to the first order Taylor expansion:
%  Areas(S+dS) =  Areas(S) + JT(S).dS      where JT is the transpose of the Jacobian J
%  We can satisfy the system of equations by taking Newton steps:
%  We calculate dS by solving the linear system: JT(S).dS = -Areas(S), and then
%  incrementing S by dS.
%  The problem is that in general there are more unknowns than there are constraints (equations)
%  this means that we have many solutions. We restrict dS to the column space of A
%  to take the shortes dS possible. So dS can be written as dS = J.dV, and dV is obtianed
%  by solving the system JT(S).dS = JT(S).J(S).dV = -Areas(S), 
%  then we increment S by J.dV
%  Note the symbol S has been subsituted by X in the code
%  Also note that the calcualtion of the Jacobian elements does not
%  necessitate subtraction of the relative areas on the real object. This
%  is because we are only calcualting differences (i.e. changes).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%
warning off;
njvals_old = inf;
njvals_new = inf;
for iter = 1:maxiter,
    if verbose, disp([iter njvals_old]);end
    t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
    
    %%% we need the triangle areas (Areas) at the current configuration.
    %%% Areas is a column vector of length N-4 (where N is the # of vertices).
    u(:) = cos(pi/2-t(:)).*cos(p(:));
    v(:) = cos(pi/2-t(:)).*sin(p(:));
    w(:) = sin(pi/2-t(:));
    crossqpr = cross([u(F(:,2))-u(F(:,1)) v(F(:,2))-v(F(:,1)) w(F(:,2))-w(F(:,1))],[u(F(:,3))-u(F(:,1)) v(F(:,3))-v(F(:,1)) w(F(:,3))-w(F(:,1))],2);
    Areas = ((sqrt(sum(crossqpr.*crossqpr,2)))./2);

    %%% we need the Jacobian (J) at the current step. J is a (N-4) x N matrix
    CHG_vec = seps*X;                           %% fill in the finite difference values
    VAL = JacPat;VAL(indJ) = X(indcol);         %% Values of the current configuration
    Cv(:) =full(VAL(indVertVal(:)));            % fill in the values without change first
    Cv(indCv) = X(indcol)+CHG_vec(indcol);
    
    %%%% Calculate the Areas with the new values given in 6-column matrix Cv(t1, t2, t3,
    %%%% p1, p2, p3) where 1, 2 and 3 are the three verices of each
    %%%% triangle
    u1 = cos(pi/2-Cv(:,1)).*cos(Cv(:,4));   u2 = cos(pi/2-Cv(:,2)).*cos(Cv(:,5));   u3 = cos(pi/2-Cv(:,3)).*cos(Cv(:,6));
    v1 = cos(pi/2-Cv(:,1)).*sin(Cv(:,4));   v2 = cos(pi/2-Cv(:,2)).*sin(Cv(:,5));   v3 = cos(pi/2-Cv(:,3)).*sin(Cv(:,6));
    w1 = sin(pi/2-Cv(:,1));                 w2 = sin(pi/2-Cv(:,2));                 w3 = sin(pi/2-Cv(:,3));
    crossqpr = cross([u2-u1 v2-v1 w2-w1],[u3-u1 v3-v1 w3-w1],2);
    areas_plus = ((sqrt(sum(crossqpr.*crossqpr,2)))./2);
    
    %%% calculate the Jacobian elements
    Jvals = (areas_plus-Areas(indrow))./CHG_vec(indcol);
    J(indJ) = Jvals;
    
   %%% test for errors and convergence %%%%%%%%%%%%%%%%%%%
    if isnan(J), disp('isnan J detected');end
    if isinf(J), disp('isinf J detected');end
    J(isnan(J)) = 0;
    J(isinf(J)) = 100;
    njvals_new = norm(Jvals);
   % if njvals_old<njvals_new, disp('converged');break;end
    njvals_old = njvals_new;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%% calculate the solution
   [dv,~,~,~,~]  = gmres(J*J',-(Areas - F_A),1, 5e0, 10);
   if sum(dv(:))==0, error('zero increment');end
   %%% and take a step
   X = X + newton_step.*J'*dv;

    %%% plot current configuration if desired
    t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
    u(:) = cos(pi/2-t(:)).*cos(p(:));v(:) = cos(pi/2-t(:)).*sin(p(:));w(:) = sin(pi/2-t(:));
    if verbose==2, clf; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');axis square;view(-276,4);zoom(1);drawnow;end
    

end
%%% plot current configuration if desired
t = [ X(1:(nvert))];p = [X(nvert+1:end)];p = mod(p,2*pi);
u(:) = cos(pi/2-t(:)).*cos(p(:));v(:) = cos(pi/2-t(:)).*sin(p(:));w(:) = sin(pi/2-t(:));
if verbose, clf; patch('Vertices',[u v w],'Faces',F,'FaceVertexCData', hsv(size(F,1)),'FaceColor','flat');axis square;view(-123,-18);zoom(1);drawnow;end

warning on;

% function Jac = fdj(fun,t,p)
% %% calculate the forward difference Jacobian for a function that takes t,p as its argument
% %% and returns a vector
% f   = feval(fun,X);
% Jac = zeros(length(f),length(X));
% CHG = sqrt(eps)*sign(X).*(abs(X));
% Xplus = X;
% for ix = 1:length(X)
%     Xplus(ix) = X(ix) + CHG(ix);
%     fplus = feval(fun,Xplus);
%     Jac(:,ix) = (fplus-f)/CHG(ix);
%     Xplus(ix) = X(ix);
% end