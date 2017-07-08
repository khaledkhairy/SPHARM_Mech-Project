classdef sh_basis < handle
    properties
        basis = 'bosh';
        L_max;
        gdim;
        t;
        p;
        w;
        Y;
        Y_P;
        Y_T;
        Y_PP;
        Y_TT;
        Y_TP;
    end
    methods (Static)
        NLK     = N_LK_bosh(L,K);
        NLK     = N_LK_nocs(L,K);
        Y       = ylk_cos_sin(L,K,phi, theta);
        [t wt]  = gaussquad(a, b, c);
        [Y P]   = ylk_cos_sin_bosh(a, b, c);
        [Y P]   = ylk_cos_sin_dphi_bosh(a, b, c, d);
        [Y P]   = ylk_cos_sin_dphiphi_bosh(a, b, c, d);
        [Y P]   = ylk_cos_sin_dtheta_bosh(a, b, c, d);
        [Y P]   = ylk_cos_sin_dthetaphi_bosh(a, b, c, d);
        [Y P]   = ylk_cos_sin_dthetatheta_bosh(a, b, c, d);
        Y = ylk_bosh(L,K,phi, theta);
    end
    
    methods
        function obj = sh_basis(L_max, gdim)
             obj.L_max = L_max;
             obj.gdim  = gdim;
            [obj.t wt]                      = obj.gaussquad(gdim, 0, pi );
            [obj.p wp]                      = obj.gaussquad(gdim,0, 2*pi);
            [obj.p obj.t]                   = meshgrid(obj.p,obj.t);
            [wp wt]                         = meshgrid(wp, wt);
            [obj.Y P]                       = obj.ylk_cos_sin_bosh(obj.p, obj.t, L_max);
            obj.Y_P                         = obj.ylk_cos_sin_dphi_bosh(obj.p, obj.t, L_max, P);
            [obj.Y_T P_T]                   = obj.ylk_cos_sin_dtheta_bosh(obj.p, obj.t, L_max, P);
            obj.Y_PP                    	= obj.ylk_cos_sin_dphiphi_bosh(obj.p, obj.t, L_max, P);
            obj.Y_TT                        = obj.ylk_cos_sin_dthetatheta_bosh(obj.p, obj.t, L_max,P_T);
            obj.Y_TP                        = obj.ylk_cos_sin_dthetaphi_bosh(obj.p, obj.t, L_max, P_T);
            obj.w                           = wp(:).*wt(:);
        end
    end
end