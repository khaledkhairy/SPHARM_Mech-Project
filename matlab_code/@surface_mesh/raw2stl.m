function raw2stl(rawFilename, stlFilename)
% RAW2STL
%   A simple utilility for converting 3D objects of RAW format to (ASCII) STL format
%   RAW2STL(RAWFILENAME, STLFILENAME) converts 
%
%   Babak Taati 
%   University of Toronto
%   Feb 2005
%   revised Jul 2010 (only style change)
%
%   Note: the information in a RAW file is not enough to resolve the normal direction (sign ambiguity)

try
    rawData = load(rawFilename); % should be an Nx9  matrix
catch
    error('invalid RAW file');
    return;
end

[triangleCount, shouldBeNine] = size(rawData);

if (shouldBeNine ~= 9) || (triangleCount== 0)
    error('invalid RAW file');
    return;
end

P1 = rawData(:,1:3); % P1(ii,:) will be the 1st vertex of the ii'th triangle
P2 = rawData(:,4:6); % P2(ii,:) will be the 2nd vertex of the ii'th triangle
P3 = rawData(:,7:9); % P3(ii,:) will be the 3rd vertex of the ii'th triangle

U = P2 - P1;    % a side of each triangle
V = P3 - P1;    % the other side of each triangle

triangleNormals = cross(U,V);   %   1- note the sign ambiguity (can't be resolved)
                                %   2- yet to be normalized ...

fid = fopen(stlFilename, 'wt'); % write to file (start)
                                
if fid == -1
    error('could not write to file');
    return;    
end

fprintf(fid, 'solid Object\n');

for ii = 1 : triangleCount
    normalizedTriangleNormal = triangleNormals(ii,:) / norm(triangleNormals(ii,:)); % ... normalize to unit length
    fprintf(fid, ' facet normal %f %f %f\n', normalizedTriangleNormal);
    fprintf(fid, '  outer loop\n');
    fprintf(fid, '   vertex %f %f %f\n', P1(ii,:) );
    fprintf(fid, '   vertex %f %f %f\n', P2(ii,:) );
    fprintf(fid, '   vertex %f %f %f\n', P3(ii,:) );    
    fprintf(fid, '  endloop\n');
    fprintf(fid, ' endfacet\n');
end

fprintf(fid, 'endsolid\n');
fclose(fid); % write to file (end)

return