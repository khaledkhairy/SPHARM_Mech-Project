function S_tr = nocs2bosh(S)
%%% convert the "new" CLKs to the old ones
%%% nocs2cs stands for conversion from coefficients that were calculated
%%% whithout the use of the Condon-Shortly phase factor to coefficients
%%% with. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_tr = nocs2cs_single(S);


%%
function Xnew = nocs2cs_single(X_o)


[xclks yclks zclks] = shp_surface.get_xyz_clks(X_o);
L_max = get_L_max(X_o);

counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        NLK_old = sh_basis.N_LK_nocs(L,K);
        NLK = sh_basis.N_LK_bosh(L,K);
        fac = NLK_old/NLK;
        fac = 1/fac;
        xclks(counter) = xclks(counter) * 1/fac;
        yclks(counter) = yclks(counter) * 1/fac;
        zclks(counter) = zclks(counter) * 1/fac;
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];