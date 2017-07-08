function obj = r_inv(obj)
%%% Calculate the Euler angles a,b,g needed for a rotational invariant shape description of X_o shape. 
%%% Uses:
%%%         get_xyz_clks, rotate_x_shp, rotate_z_shp,
%%%         parametric_rotation_objective
%%% Author: Khaled Khairy
%%% To do: Compare this (slow) method with Brechbuehler et al. 1995 method
%%%% WORKS !!! please only modify if you will make a backup copy
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_o = shp_surface.bosh2nocs(obj.X_o);
%%%%
verbose = 1;
mc_iter_p = 400;
plot_flag_p = 0;

mc_iter_o = 400;
plot_flag_o = 0;

[xc yc zc] = shp_surface.get_xyz_clks(X_o);
xpos = xc(1);ypos = yc(1);zpos = zc(1);
if isempty(obj.sf)
[X_1 res_p] = parametrization_invariance(X_o, mc_iter_o, plot_flag_p); if verbose,disp('Canonical parameterization calculated.');end
else
    [X_1 res_p sf_r] = parametrization_invariance_sf(X_o,obj.sf, mc_iter_o, plot_flag_p); if verbose,disp('Canonical parameterization calculated.');end
    obj.sf = sf_r;
end
[X_2 ang res_o] = object_invariance(X_1, mc_iter_o, plot_flag_o);if verbose, disp('Canonical object orientation calculated.');end


%%% sosi -- in case of sf, the following operations will screw up the sf
%%% orientation/reflection/inversion --- please fix soon
X_o = X_2;
ang = [];
[X_o Ihist] = fix_signs(X_2);   % Ihist encodes the inversions needed to reconstruct the original configuration
ang = [ang(:)' Ihist(:)'];
X_o = fix_inversion(X_o);
pass = shp_surface.check_canonicity(X_o);
if pass == 0,
    X_o = fix_canonicity(X_o);
end
%%%%

%%% put the object at its original position
[xc yc zc] = shp_surface.get_xyz_clks(X_o);
xc(1) = xpos;
yc(1) = ypos;
zc(1) = zpos;

%% last rearrangements
X_o = ([yc(:)' zc(:)' xc(:)']);

%%% update the object
X_o = shp_surface.nocs2bosh(X_o);
X_1 = shp_surface.nocs2bosh(X_1);
X_2 = shp_surface.nocs2bosh(X_2);

obj.X_o = X_o;
obj.X_1 = X_1;
obj.X_2 = X_2;
obj.ang = ang;
obj.res_p = res_p;
obj.res_o = res_o;
obj = obj.update;
%[X_o X_1 X_2 ang res_p res_o]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [X_o res_small]= parametrization_invariance(X_o, mc_iter, plot_flag) %#ok<INUSD>
%% parameterization invariance
t = 0;p = 0;nY1 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];
t = pi/2;p = 0;nY2 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];
%t = 0;p = pi;nY2 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];

[xc yc zc] = shp_surface.get_xyz_clks(shp_surface.nocs2cs(shp_surface.tr(X_o,1)'));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
res_small = inf;
%%% let's do a quick Monte Carlo estimation to get a good starting value
%%% IMPORTANT: random angles need to be generated as in Recipes in the
%%% future
counter = 0;
while counter<mc_iter
    ang = [rand(1)*2*pi rand(1)*pi rand(1)*2*pi];       %%% this is not a good method to get random angles (use gaussian)
    [res] = shp_surface.parametric_rotation_objective(ang,xc,yc,zc, sqrtdD, nY1, nY2,0);
    if res < res_small,res_small = res; ang_min = ang;end
    counter = counter + 1;
end
ang = ang_min;



if res_small > 1e-10
    options =   optimset(...
        'Algorithm', 'levenberg-marquardt','LargeScale', 'off', ...
        'MaxFunEvals', 1e3,'MaxIter', 300,'DiffMaxChange', 0.2,'DiffMinChange', 1e-4,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    
    plot_flag = 0;
    [ang res_small] = fminunc(@shp_surface.parametric_rotation_objective,ang,options,xc,yc,zc, sqrtdD, nY1, nY2,plot_flag);
end
 [xc yc zc] = shp_surface.get_xyz_clks(shp_surface.nocs2cs(X_o));
%%dfig(1);clf;plot_shps((X_o),3,'red');view(3);lighting gouraud;camlight;


[xc] = shp_surface.sh_rot(xc, ang(1), ang(2), ang(3));
[yc] = shp_surface.sh_rot(yc, ang(1), ang(2), ang(3));
[zc] = shp_surface.sh_rot(zc, ang(1), ang(2), ang(3));


X_o = (shp_surface.cs2nocs([xc(:)' yc(:)' zc(:)'])); %  Note the normalization to be compatible with our current basis definition
%%
%%
function [X_o res_small sf_r]= parametrization_invariance_sf(X_o, sf,mc_iter, plot_flag) %#ok<INUSD>
%% parameterization invariance with sf --- still not working properly --- more tests needed XXXXX
% sf is the sh_surface object representing a scalar field that also needs
% to be rotated
t = 0;p = 0;nY1 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];
t = pi/2;p = 0;nY2 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];

[xc yc zc] = shp_surface.get_xyz_clks(shp_surface.nocs2cs(shp_surface.tr(X_o,1)'));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
res_small = inf;
%%% let's do a quick Monte Carlo estimation to get a good starting value
%%% IMPORTANT: random angles need to be generated as in Recipes in the
%%% future
counter = 0;
while counter<mc_iter
    ang = [rand(1)*2*pi rand(1)*pi rand(1)*2*pi];       %%% this is not a good method to get random angles (use gaussian)
    [res] = shp_surface.parametric_rotation_objective(ang,xc,yc,zc, sqrtdD, nY1, nY2,0);
    if res < res_small,res_small = res; ang_min = ang;end
    counter = counter + 1;
end
ang = ang_min;



if res_small > 1e-10
    options =   optimset(...
        'Algorithm', 'levenberg-marquardt', 'LargeScale', 'off', ...
        'MaxFunEvals', 1e3,'MaxIter', 300,'DiffMaxChange', 0.2,'DiffMinChange', 1e-4,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    
    plot_flag = 0;
    [ang res_small] = fminunc(@shp_surface.parametric_rotation_objective,ang,options,xc,yc,zc, sqrtdD, nY1, nY2,plot_flag);
end
 [xc yc zc] = shp_surface.get_xyz_clks(shp_surface.nocs2cs(X_o));
%%dfig(1);clf;plot_shps((X_o),3,'red');view(3);lighting gouraud;camlight;

[xc] = shp_surface.sh_rot(xc, ang(1), ang(2), ang(3));
[yc] = shp_surface.sh_rot(yc, ang(1), ang(2), ang(3));
[zc] = shp_surface.sh_rot(zc, ang(1), ang(2), ang(3));
sf_r = sf;
for ix = 1:numel(sf)
    disp(['transforming field: ' sf_r{ix}{1} ' ' num2str(ix) ' of ' num2str(numel(sf))]);
    sf_r{ix}{2} = sh_rot(sf{ix}{2}, ang(1), ang(2), ang(3));
end
X_o = (shp_surface.cs2nocs([xc(:)' yc(:)' zc(:)'])); %  Note the normalization to be compatible with our current basis definition

function [X_out ang res_small] = object_invariance(X_in, mc_iter, plot_flag)
t = 0;p = 0;nY1 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];
t = pi/2;p = 0;nY2 = [sh_basis.ylk_cos_sin(1,-1,p, t)/sh_basis.N_LK_nocs(1,-1) sh_basis.ylk_cos_sin(1,0,p, t)/sh_basis.N_LK_nocs(1,0) sh_basis.ylk_cos_sin(1,1,p, t)/sh_basis.N_LK_nocs(1,1)];

X_in = shp_surface.nocs2cs(X_in(:)');
[xc yc zc] = shp_surface.get_xyz_clks((X_in));    %   Note the normalization so that we can use our rotation routines
Ar = [xc(2:4) yc(2:4) zc(2:4)];
[V D] = eig(Ar'*Ar);%%% diagonalize the ellipsoid
sqrtdD = sqrt(diag(D));
res_small = inf;
counter = 0;
while counter<mc_iter   %%% let's do a Monte Carlo search to get a  good starting value
    %%% IMPORTANT: random angles need to be generated as in Recipes in the
%%% future
    ang = [rand(1)*2*pi rand(1)*pi rand(1)*2*pi];
    res = shp_surface.object_rotation_objective(ang,X_in, sqrtdD,nY1, nY2, plot_flag);
    if res < res_small,res_small = res; ang_min = ang;end
    counter = counter + 1;
end
ang = ang_min;

if res_small > 1e-10    %%% perform nonlinear optimization to refine Monte Carlo guess

    options =   optimset(...
        'Algorithm', 'levenberg-marquardt', 'LargeScale', 'off', ...
        'MaxFunEvals', 1e3,'MaxIter', 300,'DiffMaxChange', 0.2,'DiffMinChange', 1e-4,...
        'Display', 'off','Diagnostics', 'off',...
        'TolFun', 1e-6,'TolX', 1e-6);
    
    plot_flag = 0;
    [ang res_small] = fminsearch(@shp_surface.object_rotation_objective,ang,options,X_in, sqrtdD, nY1, nY2,plot_flag);
end
ang(1) = mod(ang(1), 2*pi);
ang(2) = mod(ang(2), 2*pi);
ang(3) = mod(ang(3), 2*pi);
X_out = shp_surface.cs2nocs([shp_surface.rotate_shp(X_in,ang)]');
%%
function [X_o Ihist]= fix_signs(X_o)
%%% fix the signs to become canonical (we want the largest values in X to
%%% be -ve, for Y +ve and for Z -ve , as by convention). It also returns
%%% the inversion history so that we can reconstruct the original
%%% orientation from the canonical shape
[xc yc zc] = shp_surface.get_xyz_clks(X_o);
Ihist = [1 1 1];
if sign(xc(4)) == 1, Ihist(1) = -1;xc = -xc;end
if sign(yc(2)) == -1, Ihist(2) = -1;yc = -yc;end
if sign(zc(3)) == 1, Ihist(3) = -1;zc = -zc;end

X_o = [xc(:)' yc(:)' zc(:)'];
%%
function X_o = fix_canonicity(X_o)
%%% if it fails canonicity test, try to fix it
disp('Trying to fix canonicity problem.');
[xc yc zc] = get_xyz_clks(X_o);
if  sign(zc(3))== 1,
    disp('Case 1');
    [xc yc zc] = get_xyz_clks(X_o); X_o = [xc(:)' yc(:)' -zc(:)'];
    pass = check_canonicity(X_o); 
    if pass == 0, 
        [xc yc zc] = get_xyz_clks((X_o)); disp([xc(2:4) yc(2:4) zc(2:4)]);
        warning('Canonical shape transformation inconsistent');
    else
        disp('fixed !');
    end
end
%%
function X_o = fix_inversion(X_o)
%%%% if max(|C20|) is not negative, we invert the shape to make it negative.
%%%% Inversion of SHP shapes occurs by changing the sign of all
%%%% coefficients with even L > 0 coefficients 

 [xc yc zc] = shp_surface.get_xyz_clks(X_o);

if length(xc)>4
 vec = ([xc(7) yc(7) zc(7)]);
 indx = find(abs(vec) == max(abs(vec)));    % get the index i (i = x, y, or z) with the largest absolute value
 
 % if this coefficient is -ve, then let us change the sign of all L even
 % coefficients
 if sign(vec(indx))==-1, X_o = fix_xc20(X_o);end
end
 %%
 function X_o = fix_xc20(X_o)


[xc yc zc] = shp_surface.get_xyz_clks(X_o);
L_max = get_L_max(X_o);
last = 0;
for ix = 0:L_max
    start = last+1; finish = (ix + 1).^2;
    vec = start:finish;
    
    last = finish;
    if kk_iseven(ix),    % i.e. if ix is odd
        xc(vec) = -xc(vec);
        yc(vec) = -yc(vec);
        zc(vec) = -zc(vec);
    end
end
X_o = [xc(:)' yc(:)' zc(:)'];















