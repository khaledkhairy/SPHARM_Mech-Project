classdef sh_surface
    properties
        basis;
        L_max;
        gdim = 60;          % default if no basis is given
        xc;
        x;y;z;
        needs_updating = 1;
        %%% display parameters
        use_camorbit = 0;
        edge_color   = 'none';
    end
    methods(Static)
        [l,m, flags] = indices_gen(c);
        [XF X C Y_LK t p] = get_mesh(obj, nico, Y_LK, C);
        xc = bosh2nocs_xc(xc);
        xc = nocs2bosh_xc(xc);
        xc = cs2nocs_xc(xc);
        xc = nocs2cs_xc(xc);
        L_max = get_L_max_xc(xc);
        xc = tr_xc(xc, L_max);
    end
    methods
        function obj = sh_surface(L_max, basis)
            if nargin == 0,         % then make a sphere
                L_max = 6;
                obj.basis = sh_basis(10,obj.gdim);
                obj.xc = zeros((L_max + 1)^2,1);
                obj.xc(1) = 1;
                obj.L_max = L_max;
                obj = obj.update;
            end
            if nargin ==1 && numel(L_max)>1,    % then we are initializing with xc
                %disp('Assuming argument to be SH shape vector');
                obj.xc = L_max;
                obj.L_max = sqrt(numel(L_max))-1;
                obj.basis = sh_basis(obj.L_max,obj.gdim);
                obj = obj.update;
            end
            if nargin ==1 && numel(L_max) == 1,
                obj.basis = sh_basis(L_max,obj.gdim);
                obj.xc = zeros((L_max + 1)^2,1);
                obj.xc(1) = 1;
                obj.L_max = L_max;
                obj = obj.update;
            end
            if nargin ==2
               obj.basis = basis;
               obj.L_max = L_max;
               obj.xc = zeros((L_max + 1)^2,1);
               obj.xc(1) = 1;
               obj = obj.update;
            end
        end
        function [obj km kp] = flip(obj)
           %%% flip the sh_surface (top/bottom) by negating all K<0 channels
           [l,m, flags] = indices_gen(obj.xc);
           km = zeros(size(l));
           kp = zeros(size(l));
           for ix = 1:numel(l)
               if m(ix)<0, 
                   obj.xc(ix) = -obj.xc(ix); 
                   km(ix)= 1;
               end
               if m(ix)>0, kp(ix) = 1;end
           end
        end
        function [obj km kp] = flop(obj)
            %%% flop the sh_surface (front/back) by negating all K<0
            %%% channels
            [l,m, flags] = indices_gen(obj.xc);
            km = zeros(size(l));
           kp = zeros(size(l));
           for ix = 1:numel(l)
               if m(ix)<0 && (mod(m(ix),2)==0), 
                   obj.xc(ix) = -obj.xc(ix);
                   km(ix) = 1;
               end
               if m(ix)>=0 && (mod(m(ix),2)~=0),
                   obj. xc(ix) = -obj.xc(ix);
                   kp(ix) = 1;
               end

           end
        end
   end
end















