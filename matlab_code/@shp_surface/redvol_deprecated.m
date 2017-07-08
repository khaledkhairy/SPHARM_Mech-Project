function v = redvol(Area,Volume)
% input: area and volume
% output: the reduced volume, i.e. the ratio: Volume/V_sphere, where
% V_sphere is the volume of a sphere of area "Area".

% Area of a sphere = 4 * pi * r^2
r = sqrt(Area/4/pi);
% Volume of that sphere = 4/3 * pi * r^3
V_sphere = 4/3 * pi *r^3;
v = Volume/V_sphere;

