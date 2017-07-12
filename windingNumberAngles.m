function index = windingNumberAngles(vx, vy, loc, radius)
% Estimate the winding number (Poincare index) around a circle of radius
% RADIUS centred on the the row and column coordinates in LOC in the vector
% field defined by VX and VY, by computing the angular differences between
% all consecutive vectors.

% First check that calculation is not off the edge of the array
if any(floor(loc(1:2))-radius < 0) || floor(loc(1))+radius>size(vx,1) ||...
        floor(loc(2))+radius>size(vx,2)
    index = nan;
    return
end

% Convert vectors to complex values
vf = vx + 1i*vy;

% Compute the indices of the counterclockwise circuit
circRow = floor(loc(1)) + [-radius+1:radius, radius:-1:-radius+1];
circCol = floor(loc(2)) + [0:-1:-radius+1, -radius+1:radius, radius:-1:1];
circt = sub2ind(size(vf), circRow, circCol);

% Calculate the total angle difference
vf = vf(circt);
totAngDiff = sum(anglesubtract(angle([vf(2:end) vf(1)]), angle(vf)));
index = round(totAngDiff/(2*pi));

