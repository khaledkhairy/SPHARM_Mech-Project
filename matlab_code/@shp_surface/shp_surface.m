classdef shp_surface
    properties
        name = 'untitled';
        basis;
        L_max;
        gdim = 60;          % default if no basis is given
        xc;yc;zc;
        x;y;z;
        A;V;v;H;KG;T;h;S;Eb;
        X_o;X_1;X_2;ang;res_p;res_o;
        needs_updating = 1;
        sf = {};    % scalar field(s)
        %%% display parameters
        use_camorbit = 1;
        edge_color   = 'none';
        map_gen = 1;
        m_sphere = [];
        Rmx = {};   % stores rotation matrices for fast rotation once calculated
    end
    methods(Static)
        S_tr = nocs2cs(S);
        S_tr = nocs2bosh(S);
        S_tr = cs2nocs(S);
        S_tr = bosh2nocs(S);
        [l,m, flags] = indices_gen(c);
        [cp] = sh_rot(c,g, b, a);
        X_rot = rotate_shp(X_o,ang);
        pass = check_canonicity(X_o);
        res = parametric_rotation_objective(ang,xc,yc,zc,D, nY1, nY2, plot_flag);
        res = object_rotation_objective(ang,X_o,D,nY1, nY2, plot_flag);
        [xc yc zc]  = shp_sphere_gen(L_max, R,hfn);
        [X_o]       = tr(X_o, L_max);
        [xc yc zc]  = get_xyz_clks(X_o);
        [L_max]     = get_L_max(X_o);
    end
    methods
        function obj = shp_surface(arg1, arg2, arg3)
            if nargin == 0,         % then make a sphere
                L_max = 1;
                obj.basis = sh_basis(L_max,obj.gdim);
                [xc yc zc] =  obj.shp_sphere_gen(L_max, 1, @sh_basis.N_LK_bosh);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X_o);
                obj.X_o = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
                obj.L_max = L_max;
                obj = obj.update;
            end
%             if nargin ==1 && strcmp('char', class(arg1)),   %we import the surface from file
%                 obj = read_shp_surface_ascii(arg1, 60);
%             end
            if nargin ==1 && strcmp('surface_mesh', class(arg1)),  % then we are initializing by mapping a surface mesh
                % in this case input "arg1" is a surface mesh object
                % the following is inefficient and it is better to initialize the
                % the object with a predifined basis and L_max and then use
                % mesh2shp explicitly.
                disp('Using L_max 16 ! ');
                obj.L_max = 16;
                obj.basis = sh_basis(obj.L_max, obj.gdim);
                obj = mesh2shp(obj, arg1, obj.L_max);
                
            elseif nargin == 1 && strcmp('sh_basis', class(arg1)), % then we are initializing to a sphere with a basis
                obj.basis = arg1;
                L_max = arg1.L_max;
                obj.L_max = L_max;
                obj.gdim = arg1.gdim;
                [xc yc zc] =  obj.shp_sphere_gen(L_max, 1, @sh_basis.N_LK_bosh);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X_o); %#ok<PROP>
                obj.X_o = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
                obj = obj.update;
                
             elseif nargin == 1 && length(arg1)==1, % then the number is assumed to be L_max
                L_max = arg1;
                obj.basis = sh_basis(L_max,obj.gdim);
                [xc yc zc] =  obj.shp_sphere_gen(L_max, 1, @sh_basis.N_LK_bosh);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X_o); %#ok<PROP>
                obj.X_o = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
                obj.L_max = L_max;
                obj = obj.update;
                
            elseif (nargin ==1 && length(arg1)>1) && ~ischar(arg1),     % then we are initalizing an object with X_o directly
                obj = shp_surface(shp_surface.get_L_max(arg1));
                obj.X_o = arg1;    % arg1 here is really X_o
            elseif (nargin ==1 && ischar(arg1))
                if strcmp('discocyte',arg1)
                    load X_o_discocyte;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('sheet',arg1)
                    load X_o_sheet;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('plectrum',arg1)
                    load X_o_plectrum;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('stomatocyte_01',arg1)
                    load X_o_stomatocyte_01;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('stomatocyte_02',arg1)
                    load X_o_stomatocyte_02;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('stomatocyte_03',arg1)
                    load X_o_stomatocyte_03;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('bowling_pin',arg1)
                    load X_o_bowling_pin;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('disc',arg1)
                    load X_o_disc;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('bowling_pin',arg1)
                    load X_o_bowling_pin;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o;
                elseif strcmp('dros_embryo_01',arg1)
                    load X_o_dros_embryo_01.mat;
                    L = shp_surface.get_L_max(X_o);
                    basis = sh_basis(L, obj.gdim);
                    obj = shp_surface(L, basis);
                    obj.X_o = X_o; %#ok<*PROP>
                else
                    error('input was detected to be a character string: please supply a valid shape file name');
                end
                obj = update(obj);
                obj.needs_updating = 1;
            end
            
            if nargin ==2 && strcmp('surface_mesh', class(arg1)),  % then we are initializing by mapping a surface mesh
                % in this case input "arg1" is a surface mesh object and
                % arg2 is L_max
                % the following is inefficient and it is better to initialize the
                % the object with a predifined basis and L_max and then use
                % mesh2shp explicitly.
                
                obj.L_max = arg2;
                obj.gdim = 120;
                obj.basis = sh_basis(obj.L_max, obj.gdim);
                obj = mesh2shp(obj, arg1, obj.L_max);
            elseif nargin ==2 && strcmp('char', class(arg1)),   % we import the surface from file
                obj = import_shp3(obj, arg1, arg2);
                obj.name = arg1;
            elseif nargin == 2 && strcmp('sh_basis', class(arg2)),         % then make a sphere
                L_max = arg1;
                basis = arg2;
                obj.L_max = arg1;
                obj.basis = arg2;
                obj.gdim  = arg2.gdim;
                [xc yc zc] =  obj.shp_sphere_gen(L_max, 1, @sh_basis.N_LK_bosh);
                X_o = shp_surface.tr([xc(:)' yc(:)' zc(:)'], L_max);
                [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X_o);
                obj.X_o = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
                obj.needs_updating = 1;
                obj = obj.update;
            elseif nargin == 2 && strcmp('surface_mesh', class(arg2)), % map to sphere using m
                L_max = arg1;
                m = arg2;
                obj.L_max = L_max;
                obj.basis = sh_basis(obj.L_max, obj.gdim);
                obj = mesh2shp(obj, m, obj.L_max);  % note "m" is in a surface mesh object in this case
            elseif nargin == 2 && numel(arg1) ==1 && numel(arg2)==1,   % the we assume arg1 = L_max and arg2 = gdim
                L_max = arg1;
                gdim = arg2;
                basis = sh_basis(L_max, gdim);
                obj = shp_surface(basis);
            end
            if nargin==3    % then we are initializing with L, basis and a mesh
                obj.L_max = arg1;
                obj.basis = arg2;
                obj.gdim  = arg2.gdim;
                obj = mesh2shp(obj, arg3, obj.L_max);  % note "basis" is in reality a surface mesh object in this case
                obj.needs_updating = 1;
            end
            
        end
        function update_clks(obj)
            [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(obj.X_o);
        end
%         function disp_clks(obj)
%             [xc yc zc] = shp_surface.get_xyz_clks(obj.X_o);
%             cnt = 1;
%             for l = 0:obj.L_max
%                 for m = -l:l
%                     str = sprintf('%d  %d  %.2f %.2f %.2f', l, m, xc(cnt), yc(cnt), zc(cnt));
%                     disp(str);
%                     cnt = cnt + 1;
%                 end
%             end
%         end
        function obj = set.X_o(obj, X)
            obj.X_o = X;
            obj.L_max = shp_surface.get_L_max(X);
            if obj.basis.L_max ~= obj.L_max,
                obj.basis = sh_basis(obj.L_max, obj.gdim);
            end
            [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(X);
            obj.needs_updating = 1;
        end
        function obj = truncate(obj,L_max)
           obj.X_o = tr(obj.X_o,L_max);
           obj.L_max = L_max;
           [obj.xc obj.yc obj.zc] = get_xyz_clks(obj.X_o);
           new_basis = sh_basis(L_max,obj.basis.gdim);
           for ix = 1:numel(obj.sf)
               obj.sf{ix}{2}.xc = trunc_sh(obj.sf{ix}{2}.xc, L_max);
               obj.sf{ix}{2}.L_max = L_max;
               obj.sf{ix}{2}.basis = new_basis;
           end
        end
        function obj = rotate_fields(obj, g, b, a)
            % rotates all scalar fields relative to the shape 
            % by the Euler angles (in radians)
            % the convention is z-y-z
            disp(['Rotating field: ' obj.sf{1}{1} ' -- '  num2str(1) ' of ' num2str(numel(obj.sf))]);
            [obj.sf{1}{2} Rmx]= sh_rot(obj.sf{1}{2},g, b, a,0);
            if numel(obj.sf)>1
                for ix = 2:numel(obj.sf)
                    disp(['Rotating field: ' obj.sf{ix}{1} ' -- '  num2str(ix) ' of ' num2str(numel(obj.sf))]);
                    obj.sf{ix}{2}= sh_rot(obj.sf{ix}{2},g, b, a, 0, Rmx);
                end
            end
            obj.Rmx = Rmx;
        end
        function disp_sf(obj)
           for ix = 1:numel(obj.sf)
               disp([num2str(ix) ': ' obj.sf{ix}{1}]);
           end
        end
        function disp_clks(obj)
           [l m] = indices_gen(obj.xc);
           str = sprintf('\nL\tK\txc\t\tyc\t\tzc');disp(str);
           maxdisp = numel(l);
           maxN = 10;
           if maxdisp>maxN, maxdisp = maxN;disp(['showing only first' num2str(maxN)]);end
           for ix = 1:maxdisp
               str = sprintf('%d\t%d\t%.2f\t%.2f\t%.2f', l(ix), m(ix), obj.xc(ix), obj.yc(ix), obj.zc(ix));
               disp(str);
           end
        end
        function obj = set_center(obj, pos)
            %%% position is given as absolute xyz
            obj.xc(1)   = pos(1);
            obj.yc(1)   = pos(2);
            obj.zc(1)   = pos(3);
            obj.X_o     = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
        end
        function obj = scale(obj, s)
            if length(s) == 1, s = [s s s];end
            obj.xc(2:end) = obj.xc(2:end)*s(1);
            obj.yc(2:end) = obj.yc(2:end)*s(2);
            obj.zc(2:end) = obj.zc(2:end)*s(3);
            obj.X_o     = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
        end
        function obj = rotate_parameterization(obj,g, b, a, v)
                if nargin ==4, v = 0;end
                [xc yc zc] = get_xyz_clks(obj.X_o);
                disp('Rotating x');
                s = sh_surface(xc);[s, Rmx] = sh_rot(s,g, b, a,v);obj.xc = s.xc;
                disp('Rotating y');
                s = sh_surface(yc);s = sh_rot(s,g, b, a,v, Rmx);obj.yc = s.xc;
                disp('Rotating z');
                s = sh_surface(zc);s = sh_rot(s,g, b, a,v, Rmx);obj.zc = s.xc;
                obj.X_o = [obj.xc(:)' obj.yc(:)' obj.zc(:)'];
                for ix = 1:numel(obj.sf)
                    disp(['Also rotating scalar field ' obj.sf{ix}{1} ' ' num2str(ix) ' of ' num2str(numel(obj.sf))]);
                    obj.sf{ix}{2} = sh_rot(obj.sf{ix}{2},g, b, a,v, Rmx);
                end
                obj.Rmx = Rmx;
                
        end
        function obj = rotate_shps_around_self(obj,ang)
            %%% rotate(radians) the shapes described by the spherical harmonic coefficiencts S
            %%% by the Euler angles a b g (radians) around each shape's center
            %%% rotation conventions are y-z-y
            %c = c/sqrt(1/4/pi);
            %% generate rotation matrix
            a = ang(1);b = ang(2);g = ang(3);
            Rg = rot_mx(g,2);
            Rb = rot_mx(b,1);
            Ra = rot_mx(a,2);
            R = Ra*Rb*Rg;
            %% loop over the shapes and rotate them around the center
            c = [obj.xc(1) obj.yc(1) obj.zc(1)];
            obj.xc(1) = 0;obj.yc(1) = 0;obj.zc(1) = 0;
            C = [obj.xc(:)';obj.yc(:)';obj.zc(:)'];
            cr = R*C;
            cr = cr';
            tx = [cr(:,1)];tx(1) = tx(1) + c(1);
            ty = [cr(:,2)];ty(1) = ty(1) + c(2);
            tz = [cr(:,3)];tz(1) = tz(1) + c(3);
            obj.X_o = [tx; ty; tz ];
        end
        function obj = mesh2shp(obj, m, L_max)
            if nargin<3,L_max = obj.L_max;end
            disp(['Using L_max: ' num2str(L_max)]);
            if isempty(obj.m_sphere) || obj.map_gen == 1
                m = m.map2sphere;
                obj.m_sphere = m;
                obj.map_gen = 0;
            end
            if isempty(m.sf),     %i.e. no scalar field defined
                obj = obj.shp_analysis(m.X, m.t, m.p, L_max);
            else
                obj = obj.shp_analysis_with_field(m.X, m.t, m.p,m.sf, L_max);
            end
        end
        function obj = shp_analysis(obj, X, t, p, L_max)
            %%% The expansion of the three functions X = [x(t,p), y(t,p) z(t,p)] on the sphere.
            %%% [xyz] and the vectors are defined on the sphere at the spherical
            %%% coordinates t and p.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(isnan(p))==0 && sum(isnan(t))==0
                [L, K ] = shp_surface.indices_gen(1:(L_max + 1)^2);
                M = length(L); %% number of functions in expansion
                N = length(X(:,1)); %% number of data points
                A  = zeros(N, M, 'double');
                for S = 1:length(L),
                    A(:,S) = obj.basis.ylk_bosh(L(S),K(S),p(:)',t(:)')';
                end;
                [U, S, V] = svd(A, 'econ');invS = 1./(S);invS(invS==inf) = 0;
                [clks] = (V*invS) * (U'*X);
                obj.X_o = double(reshape(clks,1,M*3));
            else
                disp('nan found in phi or theta: aborting');
            end
        end
        function obj = shp_analysis_with_field(obj, X, t, p,sf, L_max)
            %%% The expansion of the three functions X = [x(t,p), y(t,p) z(t,p)] on the sphere.
            %%% [xyz] and the vectors are defined on the sphere at the spherical
            %%% coordinates t and p.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [L, K ] = shp_surface.indices_gen(1:(L_max + 1)^2);
            M = length(L); %% number of functions in expansion
            N = length(X(:,1)); %% number of data points
            A  = zeros(N, M, 'double');
            for S = 1:length(L),
                A(:,S) = obj.basis.ylk_bosh(L(S),K(S),p(:)',t(:)')';
            end;
            [U, S, V] = svd(A, 'econ');invS = 1./(S);invS(invS==inf) = 0;
            [clks] = (V*invS) * (U'*X);
            %%% fitting the scalar field(s)
            for ix = 1:length(sf)
                s = sf{ix};
                disp(['Fitting  ' s{1}]);
                sge = s{2};
                GE_new = (V*invS) * (U'*sge(:));
                sh = sh_surface(obj.L_max, obj.basis);
                sh.xc = GE_new;
                obj.sf{ix} = {s{1}, sh};
            end
            obj.X_o = double(reshape(clks,1,M*3));
        end
        function obj = rotate_true_around_y(obj, g)
            b = 0;a = 0;
            obj = rotate_parameterization(obj,g, b, a);
            [xc yc zc] = get_xyz_clks(obj.X_o);
            temp = xc;
            xc = zc;
            zc = -temp;
            obj.X_o = [xc(:)' yc(:)' zc(:)'];
        end
        function obj = rotate_true_around_z(obj, b)
            g = -pi/2;a = pi/2;
            obj = rotate_parameterization(obj,g, b, a);
            [xc yc zc] = get_xyz_clks(obj.X_o);
            temp = yc;
            yc = -xc;
            xc = temp;
            obj.X_o = [xc(:)' yc(:)' zc(:)'];
        end
        function obj = rotate_true_around_x(obj, b)
            g = 0;a = 0;
            obj = rotate_parameterization(obj,g, b, a);
            [xc yc zc] = get_xyz_clks(obj.X_o);
            temp = yc;
            yc = zc;
            zc = -temp;
            obj.X_o = [xc(:)' yc(:)' zc(:)'];
        end
        function obj = force_mirror_yz(obj)
            % enforce symmetry around y-z plane
            % works only when leading +ve Cx1-1, -ve Cy10, -ve Cz11
            [l,m, flags] = indices_gen(obj.xc);
            [xc yc zc] = get_xyz_clks(obj.X_o);
            counter = 1;
            for l = 0:obj.L_max%(int l = 0; l<this->b->L_max;l++)
                for m = -l:l%(int m = -l;m<=l;m++)
                    if (m >=0) xc(counter) = 0;end
                    if (m < 0)
                        
                        yc(counter) = 0;
                        zc(counter) = 0;
                    end
                    counter = counter + 1;
                end
            end
            obj.X_o = [xc(:)' yc(:)' zc(:)'];obj.update;
        end
        function obj = flush_field_L_max(obj, nsf, Lmax)
                nc = (Lmax + 1) * (Lmax + 1);
                sf = obj.sf{nsf}{2};
                sf.xc = tr_xc(sf.xc, L_max);
                if numel(sf.xc)<nc, sf.xc(nc) = 0;end
        end
        function obj = set_xc_field(obj, nsf, xc)
            obj.sf{nsf}{2}.xc = xc;
        end
    end
    methods
        function obj = set.A(obj,A)
            obj.A = A;
        end
        function obj = set.V(obj,V)
            obj.V = V;
        end
        function obj = set.H(obj,H)
            obj.H = H;
        end
        function obj = set.KG(obj, KG)
            obj.KG = KG;
        end
        function obj = set.T(obj, T)
            obj.T = T;
        end
        function obj = set.h(obj, h)
            obj.h = h;
        end
        function obj = set.S(obj, S)
            obj.S = S;
        end
        function obj = set.Eb(obj, Eb)
            obj.Eb = Eb;
        end
        
        function obj = set_L_max(obj, L)
            obj.L_max = L;
            obj.basis = sh_basis(L,obj.gdim);
            obj.X_o = shp_surface.tr(obj.X_o, L);
            [obj.xc obj.yc obj.zc] = shp_surface.get_xyz_clks(obj.X_o);
            if ~isempty(obj.sf)
                for ix = 1:length(obj.sf)
                    obj.sf{ix}{2}.L_max = L;
                    obj.sf{ix}{2}.gdim = obj.sf{ix}{2}.basis.gdim;
                    xc = obj.sf{ix}{2}.xc;
                    obj.sf{ix}{2}.xc = sh_surface.tr_xc(xc, L);
                end
            end
        end
        function obj = set_gdim(obj, gdim)
            obj.basis = sh_basis(obj.L_max,gdim);
            if ~isempty(obj.sf)
                for ix = 1:length(sf)
                    obj.sf{ix}{2}.gdim = obj.sf{ix}{2}.basis.gdim;
                end
            end
        end
    end
end















