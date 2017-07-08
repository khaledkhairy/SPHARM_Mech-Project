function [p, A, b, dtline,W] = longitude_calc4(x, y, z,t, A, F, L, ixN, ixS)
% Calculate longitude from diffusion as in Brechbuehler 1995
% Example:
%       p = longitude_calc(x, y, z, t, A, L, ixN, ixS)
% INPUT:
%       t: vector of latitude temperatures from lattitude_calc.m
%       A: the matrix A as obtained from lattitude_calc.m 
%       L: Cell array of link arrays
%       ixN and ixS: Identify the northpole and the south pole vertices
%       DL: Date Line represents the path from N to S
%       E and W: determine which vertices (that are neighbors of
%       DL-vertices) are on the east and which are on the west side
%       E and W are the same length, which is N_vert-2 only if every vertex
%       has a valid non-dl-neighbor. But generally their length is
%       different from DL.
% OUTPUT:
%       p: phi value associated with each vertex
%       A and b: matrix A and b as calculated in the reference above
% Author:
%       Khaled Khairy  September 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the matrix A as in the pseudo-code given on page 158
% For the north pole
links = L{ixN};     % obtain the array of indices that the North pole is liked to
for ix = 1:length(links),       % loop over links
    A(links(ix),links(ix)) = A(links(ix),links(ix))-1;      % cut the link with the North pole
end
% For the south pole
links = L{ixS};     % obtain the array of indices that the south pole is liked to
for ix = 1:length(links),       % loop over links
    A(links(ix),links(ix)) = A(links(ix),links(ix))-1;      % cut the link with the south pole
end
%%% Eliminate the ixN and ixS rows and columns in A
%A = A((1:end)~=ixN & (1:end)~=ixS,(1:end)~=ixN & (1:end)~=ixS);
%A(:,ixS) = [];A(ixS,:) = [];A(:,ixN) = [];A(ixN,:) = [];
%%%%%
A(1,1) = A(1,1) + 2;    % additional condition can be added to any row that is not a pole

% % Setup the vector b
dtline = [];
b = sparse(length(L),1);
previous = ixN;

%% determine date line based on steepest ascent in theta
nbrs = L{ixN};here = nbrs(end);       % any neighbor of the north pole
counter = 0;
maximum = 0;
%disp('Calculating dateline');
while(here~=ixS)
    %disp(counter);
    counter = counter + 1;if counter>length(L),error('could not determine p');end
    dtline = [dtline here];
    nbrs = L{here}; % get the direct neighbors of here (array of neighbor indices)
    for ix = 1:length(nbrs),
        if t(nbrs(ix)) > maximum,
            maximum = t(nbrs(ix));
            nextpos = ix;
        end
        if nbrs(ix) == previous
            prevpos = ix;
        end
    end
    previous = here;
    if exist('nextpos', 'var'), here = nbrs(nextpos(1));end
end



if isempty(dtline),
    %%   error('invalid date line !!');end
    nbrs = L{ixN};here = nbrs(1);       % any neighbor of the north pole
    counter = 0;
    maximum = 0;
    disp('Calculating dateline');
    while(here~=ixS)
        %disp(counter);
        counter = counter + 1;if counter>length(L),error('could not determine p');end
        dtline = [dtline here];
        nbrs = L{here}; % get the direct neighbors of here (array of neighbor indices)
        for ix = 1:length(nbrs),
            if t(nbrs(ix)) > maximum,
                maximum = t(nbrs(ix));
                nextpos = ix;
            end
            if nbrs(ix) == previous
                prevpos = ix;
            end
        end
        previous = here;
        here = nbrs(nextpos(1));
    end
end
%%%%%%%%%%% determine the western links
% disp('Assigning west and east labels');
%%%% now that we have the date line dl we need to determine the 
%%%% relative positions of the vertices that are neighbors of the 
%%%% points on the date line .
%%%% Let us loop over the vertices on the date line and exclude the
%%%% vertices that are on the date line itself from determination as far as
%%%% east and west is concerned. The info will be stored in 2 arrays E and
%%%% W that store the indices of the East and West vertices that are also
%%%% neighbors
S   = [x(ixS) y(ixS) z(ixS)];       % Coordinates of the South pole
N   = [x(ixN) y(ixN) z(ixN)];       % Coordinates of the North pole
E = [];W = [];DL = zeros(length(dtline),3);
here = dtline(1);
dl = dtline;
counter = 0;
%%%%%%%%%% try to walk along the dateline starting from one side only
for ix = 1:length(dtline) % as long as we are not at the south pole
    Wlinks = [];    W_count = 0;    
    if ix==1, prev = ixN;end
    if ix == length(dtline),next = ixS;else next = dtline(ix + 1);end
    % find the list of common vertices between here and prev
    
    here_links = L{here};   % get the array of links indices of here node.
    here_links = here_links(here_links ~= prev);  % the list of those links that are not prev
    here_links = here_links(here_links ~= next);
    
    prev_links = L{prev};   % get the array of links indices connected to previous node.
    prev_links = prev_links(prev_links ~= here);  % the list of those links that are not here
    
    next_links = L{next};   % get the array of links indices of next node.
    next_links = next_links(next_links ~= here);  % the list of those links that are not here
    if ix == 1
        common_links = intersect(here_links, prev_links);   % get a list of common links
%         Wlinks = [Wlinks;common_links(1)];% just assign the first one to be west.
        Wprev = common_links(1);
    end
    Wlinks = [Wlinks;Wprev];% just assign the first one to be west.
    W_count = W_count + 1;
    Wprev_links = L{Wprev(1)};
    Wprev_links = Wprev_links(Wprev_links~=here);
    Wprev_links = Wprev_links(Wprev_links~=prev);
    % to assign the next west link, it has to fulfill 2 conditions
    % 1] it must be linked to the previous west link.
    % 2] It must not be "next" or "prev".
    while ~isempty(intersect(here_links, Wprev_links))
        here_links = here_links(here_links ~=Wprev(1));
        Wprev = intersect(here_links,Wprev_links);
        %%if length(Wprev)>1,disp(Wprev);warning('Wprev should have only one element');end
        Wlinks = [Wlinks(:);Wprev(1)];    % register the index as Western
        W_count = W_count + 1;   %  increment the number of West links of this node
        %% update the valid links list
        Wprev_links = L{Wprev(1)};
        Wprev_links = Wprev_links(Wprev_links~=here);
        Wprev_links = Wprev_links(Wprev_links~=prev);
        Wprev_links = Wprev_links(Wprev_links~=next);
    end


    %%%%%
    
    W = [W;Wlinks(:)];
    dl_array = [here W_count 0]; %
    DL(ix,:) = dl_array;     % we need  this entry for modifying b
    prev = here;
    here = next;
end

%%%%%%%%%%% Modify b accordingly
%%%
% disp('Modifying b');
for ix = 1:length(dtline),      % loop over the date-line vertices
    dlarr = DL(ix,:);            % retrieve the array [dl_index nW nE]
    b(dtline(ix)) = -dlarr(2)* 2 *pi;   % dl(1) is the vertex index, dl(2) is the # of West linked vertices
end
for ix = 1:length(W),       % loop over the West vertices
    b(W(ix)) = b(W(ix)) + 2 * pi;   % increment the b vector for every occurence of a West vertex
end

% we need to delete the elements (ixN and ixS) before solving the system of
% equations
%b = b(2:end-1);
%b = b(((1:end)~= ixN) + ((1:end)~=ixS)==2);

%b(ixS) = [];b(ixN) = [];% eliminate the ixN and ixS rows 


% %  solve the linear system
% disp('Solving sparse linear system');
p = full(A\b);
% warning off;
% [L1,U1] = luinc(A,1e-3); warning on;
% [p, flag1, relres1, iter1, resvec1] = bicg(A,b, 1e-6, 1000, L1, U1);

%p = ([0 p' 0]');
