function xc = cs2nocs_xc(xc)
%%% convert the "old" CLKs to the new ones
%%% cs2nocs stands for conversion from coefficients that were calculated
%%% whith the use of the Condon-Shortly phase factor to coefficients
%%% without. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_max = round(sqrt(length(xc))-1);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
        NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
        NLK = sh_basis.N_LK_nocs(L,K);
        xc(counter) = xc(counter) * CS*NLK_old/NLK;
    end
end
