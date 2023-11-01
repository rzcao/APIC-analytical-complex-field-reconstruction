function [CTF_abe,zernikeCoeff] = findAbeFromOverlap(recFTframe,mycoord,CTF,varargin)
%findAbeFromOverlap Determine aberration from multiple KK reconstructions with different illumination angles
%   This function estimate the aberration from the image-level KK
%   reconstruction. Evenly spaced LEDs along the ring, together with more
%   images, give a better result.
%   Options
%       1. zernike mode: zernike mode to optimize for abeeration estimation
%       2. use image: select which image to use in the algorithm.
%       3. weighted: use weighted matrix in the algorithm, in which case
%                    the algorithm focuses more on larger signals.
%       4. max CTF: define the maximal spatial frequency of the CTF
%                   manually. Use this for a more accurate (subpixel)
%                   cutoff frequency estimation.
%       5. closest n pairs: for a given spectrum, pair it with the closest
%                           n other spectrums and extract their overlaps
%                           for aberration correction. By default, the
%                           program finds the closest 1 pair.
%
% The source code is licensed under GPL-3. 
%
% By Ruizhi Cao, Dec 5, 2022, Biophotonics Lab, Caltech
% Modified on Oct 5, 2023

[xsize,ysize,nBrightField] = size(recFTframe);
useimg   = 1:nBrightField;
nPairs2use_eachSpectrum = 1;
zernikeMode2Opt = 3:35;
marginPix = 2;
useWeights = true;
desiredMed2MaxWeightRatio = 2; % ratio of weight of the median signal to that of the maximal signal 
mycoord = round(mycoord);
CTFThreshold = 1*10^-2;
idx = 1;
while idx<=length(varargin)
    switch lower(varargin{idx})
        case {'zernikemode','zernike mode'}
            zernikeMode2Opt = varargin{idx+1};
            idx = idx + 2;
        case {'useimage','use image'}
            useimg = varargin{idx+1};
            idx = idx + 2;
        case {'use weight','weighted','weight'}
            useWeights = varargin{idx+1};
            idx = idx + 2;
        case {'max ctf','cutoff','max freq'}
            maxCTFUser = varargin{idx+1};
            idx = idx + 2;
        case {'find overlap','find pairs','closest n pairs','n pairs'}
            if isnumeric(varargin{idx+1}) && varargin{idx+1}>0
                nPairs2use_eachSpectrum = varargin{idx+1};
                idx = idx + 2;
            else
                error('Please specify the argument for ''closest n pairs'' with a natural number.');
            end
        otherwise
            error('Supported options are ''use image'', ''weighted'', ''cutoff'', ''n pairs'', and ''zernike mode''.');
    end
end

if any(zernikeMode2Opt<3)
    warning('The first 3 zernike mode is dropped, as those produce a worse estimation.');
    zernikeMode2Opt = zernikeMode2Opt(zernikeMode2Opt>=3);
end

if any(useimg<1) || any(useimg>nBrightField)
    error('The index of image must be positive, and should not be greater than the number of images.');
end

if nPairs2use_eachSpectrum > round(nBrightField/2)
    warning('The program is asked to pair one spectrum with more than half of the acquried spectrums. This is likely to be wrong');
end

% generate a library for the spectrum pairs
useimg = sort(useimg);
idxNextLib = zeros(length(useimg)*nPairs2use_eachSpectrum,1);
idxTemp = 1:length(useimg);
for idxPair = 1:nPairs2use_eachSpectrum
    idxNextLib(idxTemp*nPairs2use_eachSpectrum - nPairs2use_eachSpectrum + idxPair) = mod(useimg+idxPair-1,nBrightField)+1;
end
useimg = repmat(useimg(:).',[nPairs2use_eachSpectrum,1]);
useimg = useimg(:).';

xc = floor(xsize/2+1); 
yc = floor(ysize/2+1);
nZernike = length(zernikeMode2Opt);
maxCTF = find(abs(CTF(xc,yc+1:end))<CTFThreshold,1);

dcAmp = zeros(1,nBrightField);
for idx = 1:nBrightField
    dcAmp(idx) = abs(recFTframe(xc+mycoord(idx,1),yc+mycoord(idx,2),idx));
end
maxDC = max(dcAmp);

% Least-square fit of the aberration
numPhaseMeas = []; % phase measuerement in each set
resPix = min([xsize,ysize]) - 1 - 2*maxCTF;
if marginPix*2>resPix
    marginPix = floor(resPix/2);
end

if exist('maxCTFUser','var')
    if maxCTF < maxCTFUser
        warning('The thresholding-based cutoff frequency estimate is smaller than the given value. The code will overwrite this value. Please check the given cutoff freq.');
        maxCTFUser = maxCTF;
    end
    temp = linspace(-maxCTF/maxCTFUser,maxCTF/maxCTFUser,2*maxCTF+1);
else
    temp = linspace(-1,1,2*maxCTF+1);
end
[Yz,Xz] = meshgrid(temp,temp);
[theta,r] = cart2pol(Yz,Xz);
idx2use = (r<=1);
zernikeTemp = zernfun2(zernikeMode2Opt,r(idx2use),theta(idx2use));
Hz = zeros((2*maxCTF+1)^2,nZernike); % zernike operator that generates aberration
for idx = 1:nZernike
    Hz(idx2use,:) = zernikeTemp;
end
clear zernikeTemp theta r;

bd = calBoundary([xc,yc],[4*maxCTF+1,4*maxCTF+1]); 
% boundary of the 'OTF' (nonzero signal in FT of the intensity image)
[Y,X] = meshgrid(1:4*maxCTF+1,1:4*maxCTF+1);
if isa('CTF','logical')
    tempCTF = CTF(bd(1):bd(2),bd(3):bd(4));
else
    tempCTF = CTF(bd(1):bd(2),bd(3):bd(4))>=CTFThreshold;
end
bdsmall = calBoundary([2*maxCTF+1,2*maxCTF+1],[2*maxCTF+1,2*maxCTF+1]);
% boundary of the CTF (nonzero signal in the spectrum)

Hdiff = []; % operator that calculates the aberration differences
% note here Hdiff is basically D_{il} in our derivation.

phaseMeas = [];
ncolHdiff = (2*maxCTF+1)^2;
weightsVct = [];
offsetIdx = [];
countsIdxNext = 1;

for idx = useimg
    idxNext = idxNextLib(countsIdxNext);
    countsIdxNext = countsIdxNext+1;
    relShift = mycoord(idx,:) - mycoord(idxNext,:);
    maskOneside1 = (-(X-2*maxCTF-1 - mycoord(idx,1))*mycoord(idx,1) ...
                   - (Y-2*maxCTF-1 - mycoord(idx,2))*mycoord(idx,2) > -0.5*norm(mycoord(idx,:)));
    maskOneside2 = (-(X-2*maxCTF-1 - mycoord(idxNext,1))*mycoord(idxNext,1) ...
                   - (Y-2*maxCTF-1 - mycoord(idxNext,2))*mycoord(idxNext,2) > -0.5*norm(mycoord(idxNext,:)));
    overlapCTF = (tempCTF & maskOneside1) & circshift(tempCTF & maskOneside2,relShift);
    overlapCTF2 = circshift(overlapCTF,-relShift);
    
    % index of the zero-freq in the vectorized spectrum
    originSpectrum1 = (maxCTF+mycoord(idx,2))*(2*maxCTF+1) + maxCTF+mycoord(idx,1)+1;
    originSpectrum2 = (maxCTF+mycoord(idxNext,2))*(2*maxCTF+1) + maxCTF+mycoord(idxNext,1)+1;
    
    
    posIdx = find(overlapCTF(bdsmall(1):bdsmall(2),bdsmall(3):bdsmall(4)));
    negIdx = find(overlapCTF2(bdsmall(1):bdsmall(2),bdsmall(3):bdsmall(4)));
    nMeasure = length(posIdx);
    
    phaseTemp = recFTframe(bd(1):bd(2),bd(3):bd(4),idx).*tempCTF.*...
                conj(circshift(recFTframe(bd(1):bd(2),bd(3):bd(4),idxNext).*tempCTF,relShift));
    if useWeights
        weightsTemp = abs(phaseTemp(overlapCTF(:)));
        
        % reduce the phase disturbance due to non-continuous edge.
        temp = false(size(overlapCTF));
        temp(2*maxCTF+1+mycoord(idx,1),:) = true;
        temp(:,2*maxCTF+1+mycoord(idx,2)) = true;
        temp(2*maxCTF+1+mycoord(idx,1),2*maxCTF+1+mycoord(idx,2)) = false;
        weightsTemp(temp(overlapCTF(:))) = 0;
        
        weightsVct = [weightsVct;weightsTemp.*(weightsTemp>0)];
    end
    
    tempXCTF = find(sum(overlapCTF,2) >= 0.1);
    tempYCTF = find(sum(overlapCTF,1) >= 0.1);
    cropCoord = [min(tempXCTF),max(tempXCTF),min(tempYCTF),max(tempYCTF)]; % find the smallest rectangle (matrix) that contians the overlapped region
    tempWeights = log10(abs(phaseTemp(cropCoord(1):cropCoord(2),cropCoord(3):cropCoord(4)))+1); % weights for phase unwrapping
    phaseUnwrapTemp = phase_unwrapCG(angle(phaseTemp(cropCoord(1):cropCoord(2),cropCoord(3):cropCoord(4))),tempWeights./max(tempWeights(:))); % unwrap the phase before aberration extraction
    phaseRaw = angle(phaseTemp(cropCoord(1):cropCoord(2),cropCoord(3):cropCoord(4)));
    
    % force phase unwrap using multiplies of 2pi
    wrappedPhase = phaseUnwrapTemp - phaseRaw;
    x2use = -(2*pi):pi/8:(2*pi);
    N = histcounts(wrappedPhase(overlapCTF(cropCoord(1):cropCoord(2),cropCoord(3):cropCoord(4))),x2use);
    idx2use = (x2use >= -(pi+pi/4)) & (x2use <= pi+pi/4);
    N(~idx2use) = 0;
    idxPk = find(N == max(N));idxPk = idxPk(1); % use the first peak when multiple maxima are found.
    offsetPk = mean2(wrappedPhase(wrappedPhase>=(x2use(idxPk)-pi/4) & wrappedPhase<=(x2use(idxPk)+pi/4)));
    
    phaseTemp = phaseRaw + round((wrappedPhase - offsetPk)/(2*pi))*2*pi + offsetPk;
    
    % force the phase difference at the zero-freq (unaltered part of the light) to be zero
    phaseTemp = phaseTemp - phaseTemp((2*maxCTF+1)+mycoord(idx,1)-cropCoord(1)+1,...
                                      (2*maxCTF+1)+mycoord(idx,2)-cropCoord(3)+1);
    
    phaseMeas = [phaseMeas;phaseTemp(overlapCTF(cropCoord(1):cropCoord(2),cropCoord(3):cropCoord(4)))];
    Hdiff = [Hdiff;sparse(repmat(1:nMeasure,[1,2]),[posIdx.',negIdx.'],...
               [ones(1,nMeasure),-ones(1,nMeasure)],nMeasure,ncolHdiff)];
    numPhaseMeas = [numPhaseMeas,nMeasure];
    offsetIdx = [offsetIdx;originSpectrum1,originSpectrum2];
end

Hoffset = zeros(sum(numPhaseMeas),nZernike); % operator that calculates the offset
% Note: Hoffset is basically D^0_{il} in our derivation
for idx = 1:length(useimg)
    idxSt = 1+sum(numPhaseMeas(1:idx-1));
    Hoffset(idxSt:idxSt+numPhaseMeas(idx)-1,:) = repmat(Hz(offsetIdx(idx,1),:) - Hz(offsetIdx(idx,2),:),[numPhaseMeas(idx),1]);
end

if useWeights
    medWeights = median(weightsVct);
    ratioMax2Med = (maxDC.^2)/medWeights;
    factor = ceil(desiredMed2MaxWeightRatio*log10(ratioMax2Med)/(desiredMed2MaxWeightRatio-1));
    weightsVct = log10(weightsVct./max(weightsVct(:))*10^factor +1);
    
    Hoverall = weightsVct.*(Hdiff*Hz - Hoffset); 
    % Here, Hoverall is our final operator encoded with weights W*D
else
    Hoverall = Hdiff*Hz - Hoffset;
end
clear Hdiff Hoffset;

%% without considering noise: (W*) (D*Z*x) = (W*)y;
% if useWeights
%     zernikeCoeff_new = (Hoverall'*Hoverall + eye(nZernike))\(Hoverall'*(weightsVct.*phaseMeas));
% else
%     zernikeCoeff_new = (Hoverall'*Hoverall + eye(nZernike))\(Hoverall'*phaseMeas);
% end

%% compensate for noise related phase difference offset: (W*) (D*Z*x + H_e*e) = (W*)y; Herror = H_e'*H_e
sigOffset = zeros(length(useimg),1);
blockLL = zeros(length(useimg),nZernike);
weightsVctSq = weightsVct.^2; % squared weights
weightsSqSum = ones(length(useimg),1);
for idx = 1:length(useimg)
    idxSt = 1+sum(numPhaseMeas(1:idx-1));
    
    if useWeights
        blockLL(idx,:) = sum(weightsVct(idxSt:idxSt+numPhaseMeas(idx)-1).*Hoverall(idxSt:idxSt+numPhaseMeas(idx)-1,:),1);
        weightsSqSum(idx) = sum(weightsVctSq(idxSt:idxSt+numPhaseMeas(idx)-1));
        sigOffset(idx) = sum(weightsVctSq(idxSt:idxSt+numPhaseMeas(idx)-1).*phaseMeas(idxSt:idxSt+numPhaseMeas(idx)-1));
    else
        blockLL(idx,:) = sum(Hoverall(idxSt:idxSt+numPhaseMeas(idx)-1,:),1);
        sigOffset(idx) = sum(phaseMeas(idxSt:idxSt+numPhaseMeas(idx)-1));
    end
end

% % the following assumes phase offset error for each pair is independent from each other
% Herror = diag(weightsSqSum);
% Hnew = [Hoverall'*Hoverall,blockLL';blockLL,Herror];
% if useWeights
%     zernikeCoeff_temp = (Hnew + eye(nZernike+length(useimg)))\[(Hoverall'*(weightsVct.*phaseMeas));sigOffset];
% else
%     zernikeCoeff_temp = (Hnew + eye(nZernike+length(useimg)))\[(Hoverall'*phaseMeas);sigOffset];
% end
% zernikeCoeff_new = zernikeCoeff_temp(1:nZernike);

% the following assumes the phase offset error is added onto each measurement
Herror          = zeros(nBrightField,nBrightField);
sigOffset2use   = zeros(nBrightField,1);
blockLL2use     = zeros(nBrightField,nZernike);

tempIdxblockLL  = reshape(1:length(useimg),nPairs2use_eachSpectrum,[]); % group pairs whose index-wise separation is fixed. Each
                                                                        % row contains groups that share the same separation.
tempSum         = reshape(weightsSqSum,nPairs2use_eachSpectrum,[]);     % intermediate sum of contributions from each fixed distance group.
for idxPair = 1:nPairs2use_eachSpectrum
vecWeightsSum = tempSum(idxPair,:);
HerrorEach = diag((vecWeightsSum) + circshift((vecWeightsSum),[0,idxPair]));
HerrorEach = HerrorEach - circshift(diag((vecWeightsSum)),[0,idxPair])-circshift(diag((vecWeightsSum)),[0,idxPair]).';
Herror = Herror + HerrorEach;

blockLLEach = blockLL(tempIdxblockLL(idxPair,:),:);
blockLLEach = blockLLEach - circshift(blockLLEach,[idxPair,0]);
blockLL2use = blockLL2use + blockLLEach;

sigOffsetEach = sigOffset(tempIdxblockLL(idxPair,:));
sigOffsetEach = sigOffsetEach - circshift(sigOffsetEach,[idxPair,0]);
sigOffset2use = sigOffset2use + sigOffsetEach;
end

Hnew = [Hoverall'*Hoverall,blockLL2use';blockLL2use,Herror];
H2use = (Hnew(:,1:end-1)); % Due to the fact that a global phase offset added onto
                           % all measurements does not change their subtractions, we 
                           % assume that the last measurement has phase offset error equals zero.
if useWeights
    sig2use = [(Hoverall'*(weightsVct.*phaseMeas));sigOffset2use];
else
    sig2use = [(Hoverall'*phaseMeas);sigOffset2use];
end
zernikeCoeff_temp = (H2use.'*H2use + eye(nZernike + nBrightField-1))\(H2use.'*sig2use);
zernikeCoeff_new = zernikeCoeff_temp(1:nZernike);

%% generate the extracted aberration
temp = reshape(Hz*zernikeCoeff_new,[2*maxCTF+1,2*maxCTF+1]);
CTF_abe = double(CTF);
bdC = calBoundary([xc,yc],[2*maxCTF+1,2*maxCTF+1]);
CTF_abe(bdC(1):bdC(2),bdC(3):bdC(4)) = CTF_abe(bdC(1):bdC(2),bdC(3):bdC(4)).*exp(1i*temp);

%% extract the coefficients
zernikeCoeff = zeros(1,max(zernikeMode2Opt)+1);
zernikeCoeff(zernikeMode2Opt+1) = zernikeCoeff_new;

if marginPix>0
    expandFactor = (marginPix + maxCTF)/maxCTF;
    fitAbe = imresize(phase_unwrapCG(angle(CTF_abe)),[ceil(xsize*expandFactor),ceil(ysize*expandFactor)],'bilinear');
    bd = calBoundary(floor(size(fitAbe)/2+1),[xsize,ysize]);
    fitAbe = fitAbe(bd(1):bd(2),bd(3):bd(4));
    temp = CTF_abe.*(abs(CTF_abe)>0.01) + (abs(CTF_abe).*(abs(CTF_abe)<=0.01) + 10^-5).*exp(1i*fitAbe);
    CTF_abe =  CTF_abe.*(abs(CTF_abe)>0.01) + (abs(CTF_abe).*(abs(CTF_abe)<=0.01) + 10^-5).*exp(1i*medfilt2(phase_unwrapCG(angle(temp))));
end

end

