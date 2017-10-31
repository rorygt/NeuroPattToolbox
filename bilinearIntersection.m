function coords = bilinearIntersection(cell1, cell2, v)
% BILINERINTERSECTION Finds the coordinates of the intersections of the
% V-level contour lines of unit squares CELLA and CELLB, calculated using
% bilinear interpolation.
%
% Inputs: 
% CELLA must be a matrix of size [2 2] giving the values at the 4 vertices
%       of one of the data set to be interpolated.
% CELLB must have the same form as CELLA for the other data set.
% V must be a real number that defines the value of the curves in A and B
%
% Outputs:
% COORDS is a 2 x (0, 1 or 2) matrix, with the first row containing the
%   x-coordinates of all detected intersections and the second row
%   containing the y-coordinates.
%
% Contour lines are defined through bilinear interpolation by the equation
%   a*x + b*y + c*x*y + d = v. The coefficients are defined as follows:
% If the cell is the unit square with corners [f01 f11; f00 f10], then:
%  a = f10 - f00
%  b = f01 - f00
%  c = f00 + f11 - f01 - f10
%  d = f00

% Calculate bilinear interpolation coefficients for both cells
% These interpolants are based on a unit square on the cartesian plane:
%     [ (0,1)  (1,1) ]      
%     |              |      However, the MATLAB indexing of this square 
%     [ (0,0)  (1,0) ]      is [(0,1),   (0,0),   (1,1),   (1,0)]
%                           =  [cell(1), cell(2), cell(3), cell(4)]
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

a1 = cell1(4) - cell1(2);
b1 = cell1(1) - cell1(2);
c1 = cell1(2) + cell1(3) - cell1(1) - cell1(4);
d1 = cell1(2) - v;

a2 = cell2(4) - cell2(2);
b2 = cell2(1) - cell2(2);
c2 = cell2(2) + cell2(3) - cell2(1) - cell2(4);
d2 = cell2(2) - v;

coords = [];

% Catch cases where cell values are NaN
if any(isnan([a1 a2 b1 b2 c1 c2 d1 d2]))
    return
end

% Catch linear case where quadratic does not apply (divide by 0)
if abs(b2) < eps && abs(c2) < eps
    if abs(a2) < eps || (abs(b1) < eps && (abs(c1) < eps || abs(d2) < eps))
        return
    else
        x = -d2 / a2;
        y = -(a1 * x + d1) / (b1 + c1 * x);
        % Check that coordinates are in the unit square
        if y >= 0 && y < 1 && x >= 0 && x < 1
            coords = [1-y; x] + 1;
        end
        return
    end
end

% Set up quadratic equation for x
x2coeff = a1 * c2 - a2 * c1;
x1coeff = a1 * b2 - a2 * b1 - c1 * d2 + c2 * d1;
x0coeff = b2 * d1 - b1 * d2;

% Solve the equation
if abs(x2coeff) < eps
    if x1coeff == 0
        % No solutions
        return
    else
        % Equation is linear
        x = -x0coeff / x1coeff;
    end
else
    % Solve with quadratic formula
    discrim = x1coeff^2 - 4*x2coeff*x0coeff;
    if discrim<0
        % No solutions
        return
    elseif abs(discrim)<eps
        % One solution
        x = -x1coeff / (2*x2coeff);
    else
        % Two solutions
        x = (-x1coeff + [-1 1]*sqrt(discrim) ) / (2*x2coeff);
    end
end

% Solve with polynomial solver, then filter results
%x = roots([x2coeff x1coeff x0coeff])';
%x = x(imag(x) == 0);
%x = unique(x);

% Remove points where y is not defined
x((b2 + c2 * x) == 0) = [];
x(x < 0) = [];
x(x >= 1) = [];

% Find corresponding y-values by substitution
y = -(a2 * x + d2) ./ (b2 + c2 * x);

% Remove values outside the unit square
x(y < 0 | y >= 1) = [];
y(y < 0 | y >= 1) = [];

% Convert Cartesian coordinates to MATLAB row/column indexing coordinates
%     [ (0,1)  (1,1) ]         [ (1,1)  (1,2) ]
%     |              |  ---->  |              |
%     [ (0,0)  (1,0) ]         [ (2,1)  (2,2) ]
coords = [1-y; x] + 1;

end