function [obj, Rmx] = sh_rot(obj,g, b, a, verbose, Rmx)
%% rotate by Euler angles
if nargin<5, verbose = 0;end
%% first we need to convert from bosh to nocs to cs
c = sh_surface.bosh2nocs_xc(obj.xc);
c = sh_surface.nocs2cs_xc(c);
%% % rotate the spherical harmonic coefficients (normalized with CS
%%% included)
cp = c;
L_max = sqrt(length(c))-1;
[L,M] = shp_surface.indices_gen(c);
%%% construct cell array such that cell(L) has the CLKs of that L order
CL = cell(L_max+1,1);
for ix = 0:L_max,
    vec = c((ix+1)^2-2*ix -1+1:(ix+1)^2);
    CL{ix+1} = vec;
end

if nargin>5,    % then we are given the R matrix and only need to do summation and multiplication
    for ix = 1:length(L),   % loop over all indices
        if verbose, disp(['Rotating coefficient ' num2str(ix) ' of ' num2str(length(L))]);end
        l = L(ix);m = M(ix);
        cp(ix) = 0;
        counter = 0;
        clk_vec = CL{l+1};
        sq2 = sqrt(2);
        for mp = -l:l,  % loop over the summation indices
            counter = counter +1;
            clk = clk_vec(counter);
            cp(ix) = cp(ix) + clk * Rmx{ix}{counter};
        end
    end
else            % in case we are not given the R matrix
    if nargout>1, Rmx = cell(length(L),1);end
    for ix = 1:length(L),   % loop over all indices
        if nargout>1, Rmx{ix} = cell(2*L(ix)+1,1);end
        if verbose, disp(['Rotating coefficient ' num2str(ix) ' of ' num2str(length(L))]);end
        l = L(ix);m = M(ix);
        cp(ix) = 0;
        counter = 0;
        clk_vec = CL{l+1};
        sq2 = sqrt(2);
        for mp = -l:l,  % loop over the summation indices
            counter = counter +1;
            clk = clk_vec(counter);
            if (mp>0 && m>0),
                R = dlmpm(l,mp,m,b)*cos(m*g+mp*a)+ (-1)^mp*dlmpm(l,-mp,m,b)*cos(m*g-mp*a);
            elseif (mp==0 && m>0),
                R = dlmpm(l,0,m,b)*sq2*cos(m*g);
            elseif (mp<0 && m>0),
                R = (-1)^(mp+1)*dlmpm(l,mp,m,b)*sin(m*g+mp*a) +  dlmpm(l,-mp,m,b)*sin(m*g-mp*a);
            elseif (mp>0 && m==0),
                R = dlmpm(l,mp,0,b)*sq2*cos(mp*a);
            elseif (mp==0 &&m==0),
                R = dlmpm(l,0,0,b);
            elseif (mp<0 && m==0),
                R = (-1)^(mp+1) * dlmpm(l,mp,0,b)*sq2*sin(mp*a);
            elseif (mp>0 && m<0),
                R = (-1)^(m)*dlmpm(l,mp,m,b)*sin(m*g+mp*a) +  (-1)^(m+mp)*dlmpm(l,-mp,m,b)*sin(m*g-mp*a);
            elseif (mp==0 && m<0),
                R = (-1)^m*dlmpm(l,0,m,b)*sq2*sin(m*g);
            elseif (mp<0 && m<0),
                R = (-1)^(m+mp)*dlmpm(l,mp,m,b)*cos(m*g+mp*a) +  (-1)^(m+1)*dlmpm(l,-mp,m,b)*cos(m*g-mp*a);
            end
            cp(ix) = cp(ix) + clk * R;
            if nargout>1, Rmx{ix}{counter} = R;end 
        end
    end
end
%% convert back from cs to nocs to bosh
c = sh_surface.cs2nocs_xc(cp);
c = sh_surface.nocs2bosh_xc(c);
obj.xc = c;
