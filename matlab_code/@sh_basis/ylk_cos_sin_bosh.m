function [Y, P]= ylk_cos_sin_bosh(phi, theta, L_max)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max... note this version does not negate the CS phase
%% factor. and the new NLK normalization which makes it compatible with the
%% functions that calculate the derivatives later on.

if min(size(phi))==1 && min(size(theta))==1,
    phi = phi(:);theta = theta(:);
    gdimp = length(phi);
    gdimt = length(theta);
    Y = zeros(gdimp, 1, (L_max+1)^2);
    P_dim = uint8(L_max^2/2+3*L_max/2 + 1);
    P = zeros(gdimp, 1, P_dim);
    counter = 0;%
    pcounter = 0;
    for L = 0:L_max
        plk_mat = legendre(L,cos(theta(:)'));
        for K = -L:L
            counter = counter + 1;
            plk = reshape(plk_mat(abs(K)+1,:),size(theta,1),size(theta,2));
            NLK = sh_basis.N_LK_bosh(L,K);
            if K ==0 && L ==0
                pcounter = pcounter +1;
                Y(:,counter) =  NLK*plk;
                P(:,pcounter) =  NLK*plk;
            elseif K >=0,
                pcounter = pcounter+1;
                P(:,pcounter) = NLK*(plk);
                Y(:,counter) = NLK.* plk.*cos(K.*phi);
            elseif K<0,
                Y(:,counter) = NLK*plk.* sin(abs(K).*phi);
            end
        end
    end
   
    
else
    %%%%% The following has been tested and works!
    gdimp = length(phi);
    gdimt = length(theta);
    
    if (min(size(phi))==1 && max(size(theta))==1),
        phi = phi(:)';
    end
    Y = zeros(gdimp, gdimt, (L_max+1)^2);
    P_dim = uint8(L_max^2/2+3*L_max/2 + 1);
    P = zeros(gdimp, gdimt, P_dim);
    counter = 0;%
    pcounter = 0;
    for L = 0:L_max
        plk_mat = legendre(L,cos(theta(:)'));
        for K = -L:L
            counter = counter + 1;
            plk = reshape(plk_mat(abs(K)+1,:),size(theta,1),size(theta,2));
            NLK = sh_basis.N_LK_bosh(L,K);
            if K ==0 && L ==0
                pcounter = pcounter +1;
                Y(:,:,counter) =  NLK*plk;
                P(:,:,pcounter) =  NLK*plk;
            elseif K >=0,
                pcounter = pcounter+1;
                P(:,:,pcounter) = NLK*(plk);
                Y(:,:,counter) = NLK.* plk.*cos(K.*phi);
            elseif K<0,
                Y(:,:,counter) = NLK*plk.* sin(abs(K).*phi);
            end
        end
    end
end
