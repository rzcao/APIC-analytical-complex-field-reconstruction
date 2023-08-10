function [unwrappedPhase] = phase_unwrapCG(wrappedPhase,weights)
%phase_unwrapCG Performs 2D phase unwrapping using conjugate gradient
%
% This implementation is based on the work by Ghiglia and Romero (JOSA A, 1994)
% DOI: 10.1364/JOSAA.11.000107
% Input: 
%   1. warppedPhase: the phase to be unwrapped
%   2. weights: weights to be used in the unwrapping. Ideally, small weight
%               weights are assigned to places where the phase measurements
%               might be corrupted.
% Output:
%   unwrappedPhase: the unwrapped phase.

if nargin<2
    weights = ones(size(wrappedPhase));
end

% normalize the weights
weights = weights/max(weights(:));
sumW = sum(weights(:));
if any(weights(:)<0)
    error('Weights cannot be negative. Please check your input');
end

% prepare for calculating the phase
[xsize,ysize] = size(wrappedPhase);
[Y,X] = meshgrid(1:ysize,1:xsize);
denominator2use = 2*( cos(pi*(X-1)/xsize) + cos(pi*(Y-1)/ysize) -2);

% prepare for the modified discrete Laplacian operator
squaredW = weights.^2;
weightsLapX = min(squaredW(2:end,:),squaredW(1:end-1,:));
weightsLapX(end+1,:) = 0;
weightsLapY = min(squaredW(:,2:end),squaredW(:,1:end-1));
weightsLapY(:,end+1) = 0;
[dx,dy] = dGrad(wrappedPhase);
dx = wrapToPi(dx);
dy = wrapToPi(dy);

% initialization
unwrappedPhase = zeros(size(wrappedPhase));
r = dLap(dx,dy,weightsLapX,weightsLapY);
n0 = norm(r(:));
threshold = 1e-8;
clear p;

while any(r(:) ~= 0)
    dctR = dct2(r); % cosine transformed signal
    temp = dctR./denominator2use;
    temp(1,1) = dctR(1,1); % use the same bias
    z = idct2(temp);
    
    if ~exist('p','var') % same as k=1
        p = z;
        rzNow = sum(sum(r.*z));
    else
        rzNow = sum(sum(r.*z));
        beta = rzNow/rzPrev;
        p = z + beta*p;
    end
    rzPrev = rzNow;
    [dx,dy] = dGrad(p);
    dx = wrapToPi(dx);
    dy = wrapToPi(dy);
    Qpk = dLap(dx,dy,weightsLapX,weightsLapY);
    alpha = rzNow./sum(sum(p.*Qpk));
    unwrappedPhase = unwrappedPhase + alpha*p;
    r = r - alpha.*Qpk;
    
    if norm(r(:)) < threshold*n0 || sumW == numel(wrappedPhase)
        break;
    end
end


end


%% local functions
function [dx,dy] = dGrad(X) % calculate the discrete gradient
    dx = X(2:end,:) - X(1:end-1,:);
    dx(end+1,:) = 0;
    dy = X(:,2:end) - X(:,1:end-1);
    dy(:,end+1) = 0;
end

function L = dLap(dx,dy,wx,wy) % calculate the discrete Laplacian
    L = zeros(size(dx));
    if nargin == 4
        dx = wx.*dx;
        dy = wy.*dy;
    end
    Lx = dx(2:end,:) - dx(1:end-1,:);
    Ly = dy(:,2:end) - dy(:,1:end-1);
    L(1,:) = dx(1,:);
    L(2:end,:) = Lx;
    L(:,2:end) = L(:,2:end) + Ly;
    L(:,1) = L(:,1) + dy(:,1);
end