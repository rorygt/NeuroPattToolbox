function [u, v, convergenceLoop] = opticalFlowStep(im1, im2, ...
    nanIndices, surroundLocs, alpha, beta,...
    showFlag, u0, v0, im0, im3, angleFlag)
% Find the optical flow U, V between two images IM1 and IM2 from a video
% sequence.
%
% Rory Townsend, Aug 2017
% rory.townsend@sydney.edu.au

%% Fixed point iteration parameters
% Maximum fractional change between iterations to be counted as a fixed
% point
maxChange = 0.01;
% Maximum number of iterations
maxIter = 1000;
% Starting relaxation parameter for fixed point iteration
relaxParam = 1.1;
% Step to decrease relaxation parameter every iteration to ensure convergence
relaxDecStep = 0.02;
% Minimum relaxation parameter
relaxParamMin = 0.2;

%% Default inputs
if ~exist('nanIndices', 'var')
    nanIndices = find(isnan(im1));
end
if nargin<6
    alpha=0.5; % Smoothing factor
end
if nargin<7
    beta = 0.1; % Charbonnier penalty weighting
end
if nargin < 8
    showFlag = 0; % Do not show output
end
if nargin < 11
    u0 = zeros(size(im1(:,:,1))); % Use empty matrices as starting guesses
    v0 = zeros(size(im2(:,:,1)));
end

u = u0;
v = v0;

%% Find derivatives

% Spatial derivatives
[Ex1, Ey1, newNan1] = phasegradient(im1, nanIndices, angleFlag, surroundLocs);
[Ex2, Ey2, newNan2] = phasegradient(im2, nanIndices, angleFlag, surroundLocs);
Ex = (Ex1+Ex2)/2;
Ey = (Ey1+Ey2)/2;

nanIndices = unique(cat(1, nanIndices, newNan1, newNan2));
Ex(nanIndices) = 0;
Ey(nanIndices) = 0;

% Temporal derivative
if ~isreal(im1)
    Et = anglesubtract(angle(im1), angle(im2), angleFlag);
else
    if nargin < 10 || isempty(im0) || isempty(im3)
        % Take centred difference
        Et = anglesubtract(im2, im1, angleFlag);
    else
        % Use 5 point stencil
        Et = anglesubtract(1/12 * anglesubtract(im0, im3),...
            2/3 * anglesubtract(im1, im2), angleFlag);
    end
end

%fixedPointFlag = 0;
dataE = inf(size(im1));
smoothE = dataE;

% Temporarily turn off warnings to supress singular matrix warning
%warning off all

% Loop over different non-linear penalty functions until a fixed point is
% reached
for convergenceLoop = 1:maxIter
    
    lastDataE = dataE;
    lastSmoothE = smoothE;
    
    % Compute the first order error in data and smoothness
    dataE = Ex.*u + Ey.*v + Et;
    [upx, upy] = phasegradient(u, nanIndices, 0, surroundLocs);
    [vpx, vpy] = phasegradient(v, nanIndices, 0, surroundLocs);
    smoothE = upx.^2 + upy.^2 + vpx.^2 + vpy.^2;
    
    % Compute nonlinear penalty functions
    dataP = 0.5/beta*(beta^2 + dataE.^2) .^ (-1/2);
    smoothP = 0.5/beta*(beta^2 + smoothE) .^ (-1/2);
    
    % Check if data and smoothing errors have reached a fixed point
    dataEChange = abs(dataE-lastDataE) ./ abs(dataE);
    smoothEChange = abs(smoothE-lastSmoothE) ./ abs(smoothE);
    
%    % TESTING ONLY: Show loop number and convergence properties at each step
%      totError = sum(abs(dataE(:)) + alpha * abs(sqrt(smoothE(:))));
%      disp([convergenceLoop, max(max(dataPChange)), max(max(smoothPChange)), ...
%          totError])
    
    % Exit loop if fixed point has been reached
    if max(max(dataEChange)) < maxChange &&...
            max(max(smoothEChange)) < maxChange
        break
    end
    
    % Organize the discretized optical flow equations into a system of
    % linear equations in the form Ax=b. 
    [nrow, ncol] = size(im1);
    N = nrow * ncol;
    
    linear = false;
    if linear
        % Use original Horn-Schunk equations
        %gamma = 1 / alpha;
        gamma = dataP / alpha;
        delta = 4*smoothP;
        surroundTerms = surroundLocs.laplacian .* repmat(smoothP(:), 1, N);
    else
        % Use non-linear penalty function for more robust results (but
        % calculation may take more time)
        gamma = dataP / alpha;
        delta = 0;
        
        % Surrounding terms are a combination of laplacian and first
        % spatial derivative terms
        [psx, psy] = phasegradient(smoothP, [], 0, surroundLocs);        
        surroundTerms = surroundLocs.dx .* repmat(psx(:), 1, N) + ...
            surroundLocs.dy .* repmat(psy(:), 1, N) + ...
            surroundLocs.laplacian .* repmat(smoothP(:), 1, N);
    end
    
    % Calculate b vector
    b = [gamma(:) .* Et(:) .* Ex(:); gamma(:) .* Et(:) .* Ey(:)];
    % Add diagonal terms
    A = sparse(diag([-delta(:) - Ex(:).^2 .* gamma(:) ; ...
        -delta(:) - Ey(:).^2 .* gamma(:)]));
    % Add off-diagonal terms for ui-vi dependence
    uvDiag = -Ex(:) .* Ey(:) .* gamma(:);
    A = A + diag(uvDiag, N) + diag(uvDiag, -N);
    % Add other terms for surrounding locations
    A = A + [surroundTerms, sparse(N,N); sparse(N,N), surroundTerms];
    
    % Solve this system of linear equations, adding a small value along the
    % diagonal to avoid potentially having a singular matrix
    xexact = sparse(A+1e-10*eye(2*N))\b;
    
    % Reshape back to grids
    u = (1-relaxParam)*u + relaxParam*reshape(xexact(1:N), nrow, ncol);
    v = (1-relaxParam)*v + relaxParam*reshape(xexact((1:N)+N), nrow, ncol);
    
    % Gradually reduce the relaxation parameter to ensure the fixed point
    % iteration converges
    if relaxParam > relaxParamMin
        relaxParam = relaxParam - relaxDecStep;
    end
    
%    % TESTING ONLY: Show optical flow field at each step
%    quiver(u,v)
%    drawnow
end



% Turn warnings back on
%warning on all

% Plot results
if showFlag == 1
    subplot(1,2,1)
    imagesc(im1)
    hold on
    quiver(u,v)
    hold off
    subplot(1,2,2)
    imagesc(im2)
    colorbar
end

end


