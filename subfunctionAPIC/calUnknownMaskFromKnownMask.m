function [unknownMask,linearArea] = calUnknownMaskFromKnownMask(knownMask,fullMask)
%The function calculate the unknown part of the entire image/matrix using the known part, the second input can be the indicator of what is measured
%   Detailed explanation goes here

if nargin == 1
    fullMask = ones(size(knownMask));
end

unknownMask = fullMask - knownMask;
unknownMask = (unknownMask == 1);

if nargout == 2
    unknownCorrTemp = fftshift(fft2(unknownMask));
    unknownCorr = fftshift(ifft2(ifftshift(unknownCorrTemp.*conj(unknownCorrTemp))));

    crossCorrTemp = fftshift(fft2(knownMask));
    crossCorr = fftshift(ifft2(ifftshift(unknownCorrTemp.*conj(crossCorrTemp))));

    linearArea = (crossCorr>10^-9) - (unknownCorr>10^-9);
    linearArea = (linearArea == 1);
end

end

