function [rowcoords, colcoords] = criticalpointsbilinear(vx, vy, edgeDistance)
% CRITICALPOINTSBILINEAR finds the critical points in the 2D vector field
% defined by matrices VX and VY and returns their row- and column-
% coordinates. Critical point coordinates are found through bilinear
% interpolation.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

% Find piecewise linear level curves at 0 for vx and vy using MATLAB's
% built in contour function
czerox = contourc(vx, [0 0]);
czeroy = contourc(vy, [0 0]);

% Identify possible critical point sites by finding cells which contain
% level curves in both x and y components
xcrossings = processlevelcurves(czerox, size(vx), edgeDistance);
ycrossings = processlevelcurves(czeroy, size(vy), edgeDistance);
critcells = find(xcrossings & ycrossings)';

% Initialize solution vectors
coords = zeros(2*length(critcells), 2);
istore = 0;

% Iterate over all cells that may contain a critical point and detect
% precise coordinates using bilinear interpolation
for icell = critcells
    
    % Define the corners of the cell
    [irow, icol] = ind2sub(size(vx)-1, icell);
    corners = sub2ind(size(vx), [irow irow; irow+1 irow+1], ...
        [icol icol+1; icol icol+1]);
    
    if ~any(isnan(vx(corners(:))) | isnan(vy(corners(:))))
        % Find crossings through bilinear interpolation
        icoords = bilinearIntersection(vx(corners), vy(corners), 0);
        numcrossings = size(icoords,2);
        
        if ~isempty(icoords)
            
            % Scale coordinates from the unit 2x2 matrix input to
            % bilinearIntersection back to their position in the full matrix
            icoords(1,:) = icoords(1,:) + irow - 1;
            icoords(2,:) = icoords(2,:) + icol - 1;
            
            % Store results
            coords((1:numcrossings)+istore, :) = icoords';
            istore = istore + numcrossings;
            
        end
    end
    
end

% Trim zeros off result vectors
coords = coords(1:istore, :);

% Remove repeated coordinates
coords = unique(coords, 'rows');

rowcoords = coords(:,1);
colcoords = coords(:,2);


function crossingMatrix = processlevelcurves(contourc_output, matSize, edgeDist)
% Takes the output of MATLAB's contourc and uses it to find the corners of
% all cells which are crossed by the 0-level curve.

crossingMatrix = zeros(matSize-1);

% Find delimiters that mark the start of each contour line and remove them
delimiters = contourc_output(1,:) == 0;
edgePoints = contourc_output(:,~delimiters);

% Find the coordinates of every cell that any level curve crosses
colcoords = [floor(edgePoints(1,:)-10*eps), floor(edgePoints(1,:)+10*eps)];
rowcoords = [floor(edgePoints(2,:)-10*eps), floor(edgePoints(2,:)+10*eps)];

% Remove points outside the matrix
badrow = rowcoords<1+edgeDist | rowcoords>=matSize(1)-edgeDist;
badcol = colcoords<1+edgeDist | colcoords>=matSize(2)-edgeDist;
bad = badrow | badcol;

crossingMatrix(sub2ind(matSize-1, rowcoords(~bad), colcoords(~bad))) = 1;

