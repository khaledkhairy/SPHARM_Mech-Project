classdef surface_mesh
    properties
        X;
        F;
        E;
        L;
        face_memb;
        isclosed_shape = 1;
        needs_updating = 1;
        needs_edge_info = 1;
        needs_map2sphere = 1;
        euler_violation = 0;
        A; V; v; F_areas; h; H; Eb; da;dA;quality;
        sf = {};    % scalar field
        %%% various configurations
        laplacian_smooth_iter = 10;
        laplacian_smooth_beta = .001;
        meshresample_keepratio= 0.7;
        %%% properties specific to the mapping to the sphere part
        ixN = 0;
        ixS = 0;
        bijective_plot_flag = 1;
        mapping_plot_flag = 2;
        newton_niter = 200;
        newton_step  = 0.05;
        optimization_method = 1;      % 1 will optimize based on area only, 2 will also calculate angles
        t;p;
        objfun = @qarea;
        %%% display parameters
        use_camorbit = 1;
    end
    methods(Static)
        %%% mesh manipulation functions
        [vertex1,face1] = subdivide_XF(vertex,face,nsub);
        %%% shape generators
        [X,F]       = sphere_mesh_gen(n);
        %%% map2sphere functions
        [t,p,dtline, W] = bijective_map_gen(X, F, L, plotflag, ixN, ixS);
        [t, Ap, b] = latitude_calc(L, ixN, ixS);
        [p, A, b, dtline,W] = longitude_calc(x, y, z,t, A, F, L, ixN, ixS);
        [t,p] = newton_steps_02( t,p,F,F_areas,newton_step,maxiter,verbose);
        mesh2raw(X,F, fn);
        raw2stl(fn1,fn2);
    end
    methods
        function obj = surface_mesh(X,F)
            if nargin==2,
                if size(X,1) ==3, X = X';end
                if size(F,1) ==3, F = F';end
                obj.X = X;
                obj.F = F;
            elseif nargin==0,
                [obj.X obj.F] = surface_mesh.sphere_mesh_gen(3);
            end
        end
        function obj = translate_to_center_of_mass(obj)
            obj.X(:,1) = obj.X(:,1) - mean(obj.X(:,1));
            obj.X(:,2) = obj.X(:,2) - mean(obj.X(:,2));
            obj.X(:,3) = obj.X(:,3) - mean(obj.X(:,3));
        end
        function obj = edge_info(obj)
            if obj.needs_edge_info
                nvert = length(obj.X);
                nfaces = length(obj.F);
                if size(obj.F,1)==3, obj.F = obj.F';end
                if size(obj.X,1)==3, obj.X = obj.X';end
                if (nvert*2-4~=nfaces), 
                    disp('Mesh does not represent a closed shape (for triangulations)');
                    obj.isclosed_shape = 0;
                end
                E = zeros(nvert,2);counter = 0;L = {};face_memb = {};
                % determine the links L and the member faces
                for ix = 1:length(obj.X),   % loop over the vertices
                    fmemb = ismember(obj.F, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
                    face_memb{ix} = find(sum(fmemb,2));
                    [ig, jg] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
                    links = [];
                    for ik = 1:length(ig),
                        links = [links obj.F(ig(ik),:)];
                    end
                    L{ix} = unique(links(links~=ix));   % only record the links that are not ix
                    % create the list of all edges
                    llinks = L{ix};
                    for jx = 1:length(llinks),
                        counter = counter + 1;
                        E(counter,:) = [ix llinks(jx)];
                    end;
                end;
                obj.E = E;
                obj.L = L;
                obj.face_memb = face_memb;
                if (length(obj.X)-length(E)/2+length(obj.F)~=2),
                    obj.euler_violation = 1;
                else
                    obj.euler_violation = 0;
                end
                obj.needs_edge_info = 0;
            end
        end
        function obj = map2sphere(obj, t, p)
            if obj.needs_map2sphere
                obj = obj.props;
                if nargin==1 && obj.optimization_method == 1,
                    if size(obj.F,2)==3
                        disp('Performing bijective mapping using solution of heat diffusion equations.');
                        %% Here we use Brechbuehler's idea of solving the heat diffusion equations for rough mapping
                        obj = obj.edge_info;
                        if (obj.ixN)
                            
                            [obj.t,obj.p,~,~] = surface_mesh.bijective_map_gen(obj.X, obj.F, obj.L, obj.bijective_plot_flag, obj.ixN, obj.ixS); % %%%% Generate Bijective Map
                        else
                            [obj.t,obj.p,~,~] = surface_mesh.bijective_map_gen(obj.X, obj.F, obj.L, obj.bijective_plot_flag); % %%%% Generate Bijective Map
                        end
                    elseif size(obj.F,2)==4        % assume quad mesh
                        disp('Performing quad mesh bijective mapping');
                        [obj.t,obj.p, E, L, face_memb, A,b,Aq, bq, dl]  =...
                            qbijective_map_gen(obj.X,obj.F, obj.ixN, obj.ixS);
                        dfig;[u, v, w] = kk_sph2cart(obj.t,obj.p,ones(size(obj.p)));plot_state(u,v,w,obj.F);drawnow
                    end
                else
                    %obj.t = t;
                    %obj.p = p;
                end
                if obj.optimization_method == 1
                    if size(obj.F,2)==3
                    %% using simple Newton iterations as in Khairy and Howard 2008 (MEDIA)
                    disp('Performing area preservation optimization only using gmres. Using planar approximation for triangles.');
%                    warning('under testing --- works best with high quality meshes');
                    [obj.t obj.p] = obj.newton_steps_02( obj.t, obj.p,obj.F,obj.F_areas,...
                        obj.newton_step,obj.newton_niter, obj.mapping_plot_flag);%#ok<PROP> %%%% Constrained optimization for uniform parametrization
                    obj.needs_map2sphere = 0;
                    elseif size(obj.F,2)==4
                        disp('Mapping optimization implementation for quad meshes is still in progress.');
                        obj.F_areas = qarea(obj.X, obj.F);
                        
                        [obj.t, obj.p] = ...
                            quad_optimize_02(obj.objfun,obj.t, obj.p,obj.F, ...
                            obj.F_areas, obj.newton_step, obj.newton_niter);
                    end
                    
% %                 elseif obj.optimization_method == 2
% %                     %% we need to also preserve angles (not only areas), i.e. resist shear deformation
% %                     disp('Not calculating starting value for bijective mappin. Performing Gu Wang 2004.');
% %                     [t p] = conformal_spherical_mapping_Gu(obj.X,obj.F);
% %                     obj.t = t;
% %                     obj.p = p;
% %                     obj.needs_map2sphere = 0;
% %                 elseif obj.optimization_method == 3
% %                     disp('Performing Gu Wang 2004, using bijective map as starting values');
% %                     [t p] = conformal_spherical_mapping_Gu(obj.X,obj.F, obj.t, obj.p);
% %                     obj.t = t;
% %                     obj.p = p;
% %                     obj.needs_map2sphere = 0;
% %                 elseif obj.optimization_method == 4
% %                     %% preserve angles and use large scale matlab optimization
% %                     disp('Performing both area and angle preservation optimization on the sphere using geodesic triangle calculation and matlabs large scale optimization');
% %                     obj = optimize_on_sphere(obj);
% %                     obj.needs_map2sphere = 0;
                else
                    warning('No optimization performed');
                end
                
            end
        end
        function obj = subdivide(obj, nsub)
            if nargin ==1, nsub = 1;end
            [obj.X obj.F] = surface_mesh.subdivide_XF(obj.X, obj.F, uint8(nsub));
            obj.needs_updating = 1;
            obj.needs_edge_info = 1;
            obj.needs_map2sphere = 1;
        end
        function obj = laplacian_smooth(obj)
            obj = obj.edge_info;
            for iter = 1:obj.laplacian_smooth_iter,
                for ix = 1:length(obj.X),
                    vec = obj.X(ix,:);
                    vix = obj.L{ix};
                    for ld = 1:30,
                        obj.X(ix,:) = obj.X(ix,:) + obj.laplacian_smooth_beta* sum(obj.X(vix,:)-vec(ones(length(vix),1),:))/length(vix);
                    end
                end
            end
        end
        function obj = optimize_mesh(obj)
            obj = obj.props;
            %%%%%%%%%%%% simplify using CGAL (computational geometry algorithms library), and then clean
            [obj.X, obj.F] = meshresample(obj.X,obj.F,obj.meshresample_keepratio);
            [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'duplicated');
            [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'isolated');
            [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'deep', '99');
            facecell= finddisconnsurf(obj.F);        % find disconnected surfaces
            if size(facecell,2)>1,
                disp(' ****** multiple surface fragments found, picking largest!');
                obj.F      = maxsurf(facecell);         % extract the largest surface
                %%%%%%%%%%%% simplify using CGAL (again) and clean mesh
                [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'duplicated');
                [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'isolated');
                if~isdeployed,[obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'deep', '99');end
            end
            if (length(obj.X)*2-4~=length(obj.F)),
                warning('KK: Mesh does not appear to represent a closed shape!');
                %%% use convhull to close the mesh
                obj.F = convhulln(X);
                [obj.X, obj.F] = meshcheckrepair(obj.X,obj.F,'deep', '99');
                [obj.X, obj.F] = subdivide(obj.X,obj.F,2);
                if (length(obj.X)*2-4~=length(obj.F)), warning('********* KK: Mesh does not appear to represent a closed shape!*******');end
            end
            obj.needs_updating = 1;
            obj.needs_edge_info = 1;
            obj.needs_map2sphere = 1;
        end
        function plot_fast(obj,c)
            patch('Vertices', obj.X, 'Faces', obj.F,'FaceColor', c,'EdgeColor','none','FaceAlpha',1);
        end
        function plot_quality(obj)
            obj = obj.props;
            patch('Vertices', obj.X, 'Faces', obj.F,'FaceColor', 'flat', 'FaceVertexCData',obj.quality,'EdgeColor','none','FaceAlpha',1);
            axis off;axis equal;
            colorbar;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
        end
        function [obj qF] = plot_map_quality(obj)
            obj = obj.props;
            [x y z] = kk_sph2cart(obj.t,obj.p, 1);
            X = [x(:) y(:) z(:)];
            C = obj.F;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            A = 0; V = 0; H = 0; h = 0; v = 0; F_areas = 0;
            u = X(:,1); v = X(:,2); w = X(:,3);
            x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
            q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
            crossqpr = cross(q,r,2);
            twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
            A = sum(twoA)/2;                    % this is the total area
            F_areas = twoA/2;                   % this is the vector of face areas
           
            %%% for triangle quality determination
            p = [x3-x2 y3-y2 z3-z2];
            d1 = sqrt(q(:,1).^2 + q(:,2).^2 + q(:,3).^2);
            d2 = sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2);
            d3 = sqrt(p(:,1).^2 + p(:,2).^2 + p(:,3).^2);
            qF = 4*F_areas*sqrt(3)./(d1.^2 + d2.^2 + d3.^2);
            patch('Vertices', X, 'Faces', C,'FaceColor', 'flat', 'FaceVertexCData',qF,'EdgeColor','none','FaceAlpha',1);
            axis off;axis equal;
            colorbar;
            if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
        end
        function obj = add_XF(obj, X2, F2)
            % assuming F2 to be n x 3 and X2 as well
            g  = size(obj.X,1);
            g = [g g g];
            g = g(ones(size(F2, 1),1),:);   
            obj.F = [obj.F;F2 + g];
            obj.X = [obj.X;X2];
            obj.needs_updating = 1;
            obj.needs_edge_info = 1;
            obj.needs_map2sphere = 1;
        end
        function obj = merge_XF(obj, X2, F2)
            
            disp('Merging of two meshes is not implemented yet!!');
            %%% determine the closest two points.
            
            %%% perform local 3D metamorphosis until the two points meet at
            %%% their center point
            
            
            
            %%% Delete the points and perform local remeshing
            
            
            %%% widen the tube by translating the local mesh outwards under
            %%% the effect of a force emanating from the line going
            %%% through the tube
            
            
            % merge the new sets of vertices and faces:assuming F2 to be n x 3 and X2 as well
% %             g  = size(obj.X,1);
% %             g = [g g g];
% %             g = g(ones(size(F2, 1),1),:);   
% %             obj.F = [obj.F;F2 + g];
% %             obj.X = [obj.X;X2];
% %             %%% now we need to find a proper hull that engulfs all the
% %             %%% vertices and produces a new mesh
% %             [ center, radii, evecs, v ] = ellipsoid_fit( obj.X);
% %             obj.F = convhulln(double(obj.X));
% %             obj = optimize_mesh(obj);
% %             %%
% %             obj.needs_updating = 1;
% %             obj.needs_edge_info = 1;
% %             obj.needs_map2sphere = 1;
        end
        function obj = baloon(obj)
            %%%% simulate blowing up of a baloon  
            niter =300;
            Afac = 2.0;
            maxdX = inf;%0.1;
            maxK = 0.6;
            
            pre = shear_stretch_precalc(obj.X, obj.F);
            [E,L,face_memb] = edge_info(obj.X,obj.F);
            [A, V, v, F_areas_o, h, H, wb, da, N, K, k_g, dA] = triangulated_props(obj.X, obj.F, 0);
            
            %Ko = 1/(A*Afac/4/pi);    % target Gaussian curvature
            X = obj.X;
            F = obj.F;
            
            vfac = zeros(size(X,1),1);
            
            for ix = 1:niter
                [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(obj.X, obj.F, 0);
                %K(K>maxK) = maxK;
                %K = abs(K-Ko);
                %H = H-min(H);H = H/norm(H);
                %H = H;%/norm(H);
                %[E_shear ai beta] = shear_calc(X,F,pre);
                %trifacb  = 1-full(beta);
                %trifaca  = 1-F_areas/F_areas_o;
                %clf;subplot(2,1,1);hist(ai);subplot(2,1,2);hist(beta);drawnow;
% %                 for(ix = 1:size(X,1))
% %                     vfac(ix) = mean(trifacb(face_memb{ix}));
% % %                     vfac(ix) = vfac(ix) + mean(trifaca(face_memb{ix}));
% %                 end
                fac  = 0.01;%(vfac)/2;
                mag = H.*fac;
                mag(mag>maxdX) = maxdX;
                X = X + [mag mag mag].*N;
                %dfig(1);cla;plot(vfac);drawnow;
                %dfig(2);cla;hist(H);drawnow;
                dfig(3);cla;
                patch('Vertices', X, 'Faces', F,'FaceColor', 'interp', 'FaceVertexCData', K, 'FaceAlpha',0.4);
                view(3);axis equal;colorbar;
                view(37.5, 90);
                drawnow;
                cameramenu;
                
%                 dfig(4);cla;
%                 quiver3(X(:,1), X(:,2), X(:,3), N(:,1), N(:,2), N(:,3));
%                 drawnow;
                
            end
        end
    end
end
























