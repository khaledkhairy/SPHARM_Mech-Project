function obj = rotate(obj, angles)
%% rotate the shape around center of mass
%% angle are the Euler angles in the z-x-z convention in degrees
%  see:
% http://en.wikipedia.org/wiki/Euler_angles#Euler_angles_as_composition_of_Euler_rotations
% for definition
all_X = obj.X;
convert_to_cell = 0;
if ~iscell(all_X), convert_to_cell =1;X{1} = all_X; all_X = X;end
X = cell2mat(all_X');   % concatenate all X
%% center of mass
cm = sum(X,1)./length(X);cm = cm(:)';
% generate the rotation matrix
a = angles(1);
b = angles(2);
g = angles(3);
R = [(cosd(a)*cosd(g)-sind(a)*cosd(b)*sind(g)) ...
               (-cosd(a)*sind(g)-sind(a)*cosd(b)*cosd(g)) ...
                           sind(b)*sind(a);...
     (sind(a)*cosd(g)+cosd(a)*cosd(b)*sind(g)) ...
               (-sind(a)*sind(g)+cosd(a)*cosd(b)*cosd(g)) ...
                           -sind(b)*cosd(a);...
      sind(b)*sind(g)...
               sind(b)*cosd(g)...
                            cosd(b)];
% translate the center of mass to the origin
X = X-cm(ones(1,size(X,1)), :);
% apply the rotation
Xp = (X*R);
% translate back the center of mass to where is was
obj.X = Xp+cm(ones(1,size(X,1)), :);
