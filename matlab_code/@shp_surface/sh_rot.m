function [cp] = sh_rot(c,g, b, a)
%%% rotate the spherical harmonic coefficients (normalized with CS?)
cp = c;
L_max = sqrt(length(c))-1;
[L,M] = shp_surface.indices_gen(c);
%%% construct cell array such that cell(L) has the CLKs of that L order
CL = {};
for ix = 0:L_max,
    vec = c((ix+1)^2-2*ix -1+1:(ix+1)^2);
    CL{ix+1} = vec;
end
for ix = 1:length(L),   % loop over all indices
    l = L(ix);m = M(ix);
    cp(ix) = 0;
    counter = 0;
    clk_vec = CL{l+1};
    if sum(abs(clk_vec))>0,        %%%% todo : Check here whether clk_vec is all zero or not
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
        end
    end%%%%% todo : end the if-statement where we check clk_vec == 0 ?
end
%