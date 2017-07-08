function xc = bosh2nocs_xc(xc)
%%% convert the "new" CLKs to the old ones
%%% nocs2cs stands for conversion from coefficients that were calculated
%%% whithout the use of the Condon-Shortly phase factor to coefficients
%%% with. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_max = sh_surface.get_L_max_xc(xc);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        NLK_old = sh_basis.N_LK_nocs(L,K);
        NLK = sh_basis.N_LK_bosh(L,K);
        fac = NLK_old/NLK;
        xc(counter) = xc(counter) * 1/fac;
    end
end
