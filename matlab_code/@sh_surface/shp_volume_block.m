function [block, in, outline] = shp_volume_block(obj, dim, nico, Y_LK, C)
%%%  generate the block information for an image volume with the object in it
%%%  to generate Y_LK call get_mesh beforehand once with:
%%%  [XF X C Y_LK] = get_mesh(obj, nico)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==1, 
    nico = 5;
    [XF X C Y_LK t p] = obj.get_mesh(nico);
    minx = floor(min(X(:,1)));maxx = ceil((max(X(:,1))));
    miny = floor(min(X(:,2)));maxy = ceil(max(X(:,2)));
    minz = floor(min(X(:,3)));maxz = ceil(max(X(:,3)));
    
    X(:,1) = X(:,1) - minx + 1;
    X(:,2) = X(:,2) - miny + 1;
    X(:,3) = X(:,3) - minz + 1;
    
    minx = floor(min(X(:,1)));maxx = ceil((max(X(:,1))));
    miny = floor(min(X(:,2)));maxy = ceil(max(X(:,2)));
    minz = floor(min(X(:,3)));maxz = ceil(max(X(:,3)));

    dim(1) = maxy+1;
    dim(2) = maxx+1;
    dim(3) = maxz+1;
end  

if nargin ==2, nico = 4; [XF X C] = obj.get_mesh(nico);end  % slow, needs to build basis
if nargin ==5, 
%     disp([size(Y_LK) size(C) nico]);
    [XF X C] = obj.get_mesh(nico, Y_LK, C);
end     % faster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = X(:,1);X(:,1) = X(:,2);X(:,2) = temp;

minx = floor(min(X(:,1)));maxx = ceil((max(X(:,1))));
miny = floor(min(X(:,2)));maxy = ceil(max(X(:,2)));
minz = floor(min(X(:,3)));maxz = ceil(max(X(:,3)));

% % if minx <1,minx = 1;end;if maxx>dim(1), maxx = dim(1);end
% % if miny <1,miny = 1;end;if maxy>dim(2), maxy = dim(2);end
% % if minz <1,minz = 1;end;if maxz>dim(3), maxz = dim(3);end

if minx < 1 || maxx>dim(1) || miny <1 || maxy>dim(2) || minz <1 || maxz>dim(3),
    block = [];
    in = [];
    outline = [];
else
    
    
    %%
    if nargout == 3,    % then calculate the object outlines as well
        x1 = round(X(:,1));
        x2 = round(X(:,2));
        x3 = round(X(:,3));
        ixval = find(x1<(dim(1)+1) & x1>0);x1 = x1(ixval);x2 = x2(ixval);x3 = x3(ixval);
        iyval = find(x2<(dim(2)+1) & x2>0);x1 = x1(iyval);x2 = x2(iyval);x3 = x3(iyval);
        izval = find(x3<(dim(3)+1) & x3>0);x1 = x1(izval);x2 = x2(izval);x3 = x3(izval);
        
        outline = sub2ind(dim,x1, x2, x3);
    end
    
    %% calculate the full volume indices
    [xi yi zi] = meshgrid(minx:maxx, miny:maxy, minz:maxz);
    eps = 0.5;
    in = zeros(length(minx:maxx), length(miny:maxy), length(minz:maxz), 'uint8');
    im = zeros(length(minx:maxx), length(miny:maxy), 'uint8');
    counter = 1;
    se = strel('disk',1);
    for zix = minz : maxz
        c = (X(:,3)<(zix+eps)) & (X(:,3)>(zix-eps)); % ones correspond to indices on the contour
        imix = sub2ind(size(im), round(X(c,1))-minx+1, round(X(c,2))-miny+1);
        im(:) = 0;
        im(imix) = 1;
        in(:,:,counter) = imfill(imclose(im, se));
        counter = counter + 1;
%         figure(1);cla;imshow(im*255);drawnow;
        
    end
    in = logical(in);
    
    block.minx = minx;
    block.miny = miny;
    block.minz = minz;
    
    block.maxx = maxx;
    block.maxy = maxy;
    block.maxz = maxz;
    if isempty(in),error('empty object');end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
