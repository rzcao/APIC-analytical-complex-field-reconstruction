function [recFTframe,varargout] = recFieldKK(imStack,mycoord,varargin)
% use kramers-kronig to do field reconstruction
% Requried inputs:
%   1. imStack: The image acquired under critial-angle illuminations
%   2. mycoord: Coordinate of the illumination vector in terms of pixel in
%       the spatial frequency space. i.e. the pixel-wise shift of the
%       zero-frequency after Fourier transform (FT).
% Optional arguments:
%   1. CTF: enable CTF in the reconstruction. When enabled, the
%       reconstructed field whose spatial frequency lies outside the CTF
%       support will be set to zero.
%   2. Support: the support of the image's FT. When used, we only keep the
%       image's FT within the support (To reduce noise).
%   3. Normalizaiton: whether to normalize the acquired images such that
%       they have the same effective intensity (mean value). Disabled by
%       default.
%   4. Padding: Choose the zero-padding factor of the FT of the images.
%       Default padding factor is 3. This should integer that is larger
%       than 1.
%   5. Wiener: whether to use Wiener filter to mitigate noise. Disabled by
%       default.
%   6. noise floor: set the noise floor of the acquired image. It helps to
%       prevent image from showing weired phases in the dark region where
%       the SNR is too low for the algorithm to extract the right phases.
%       The noise floor is set automatically when no number is specified.
%   7. use data intensity: treat the square root of measured intensity as
%       the ground truth for amplitude. This is enabled by default.
%
% By Ruizhi Cao, Dec 7, 2022, Biophotonics Lab, Caltech
% Modified on July 20, 2023


padfactor       = 3; % default padding factor
useCTF          = false;
useSupport      = false;
normIntensity   = false;
useNoiseFilter  = false;
autoNoiseFloor  = true;
useDataIntensity= true;
nOpt = length(varargin);
idx = 1;
while idx <= nOpt
    switch lower(varargin{idx})
        case 'ctf'
            CTF = varargin{idx+1};
            useCTF = true;
            idx = idx+2;
        case {'support','otf','mask'}
            maskFilt = varargin{idx+1};
            useSupport = true;
            idx = idx+2;
        case {'normalization','norm'}
            normIntensity = varargin{idx+1};
            idx = idx+2;
        case {'pad','padding','padfactor','factor'}
            padfactor = varargin{idx+1};
            idx = idx+2;
        case {'wiener','reg','regularization','wiener filer','wienerfilter'}
            useNoiseFilter = varargin{idx+1};
            idx = idx+2;
        case {'noise floor','noisefloor','noise'}
            regImage = varargin{idx+1};
            autoNoiseFloor = false;
            idx = idx+2;
            if regImage<=0
                autoNoiseFloor  = true;
                warning('Noise level must be positive. Program uses its default noise level.');
            end
        case {'data intensity','replace intensity','use data intensity'}
            useDataIntensity = (varargin{idx+1}~=0);
            idx = idx+2;
        otherwise
            error('Supported options are ''pad'', ''normalization'', ''wiener'', ''noise floor'', ''CTF'' and ''support''.');
    end
end

if rem(padfactor,1)~= 0
    warning('Padding factor must be an integer. It is rounded up.');
    padfactor = ceil(padfactor);
end

if padfactor<1
    error('Pad factor must be larger than 1.');
end

if useNoiseFilter && ~useCTF
    error('To suppress noise, CTF must be given to the program to estimate noise''s amplitude.');
end

[xsize,ysize,nBrightField] = size(imStack);
xc = floor(xsize/2+1);
yc = floor(ysize/2+1);

% normalize the effective intensity of each image
if normIntensity
    meanIFrame = mean(mean(imStack,1),2);
    meanIFrame = meanIFrame./mean(meanIFrame);
    for idx = 1:nBrightField
        imStack(:,:,idx) = imStack(:,:,idx)/meanIFrame(idx);
    end
end

% assign spaces for the padded images and FT array
paddedFT = zeros(xsize*padfactor,ysize*padfactor);
imStackPad = zeros(xsize*padfactor,ysize*padfactor,nBrightField);

xcpad = floor(xsize*padfactor/2+1);
ycpad = floor(ysize*padfactor/2+1);
bdpad = calBoundary([xcpad,ycpad],[xsize,ysize]);

for idx = 1:nBrightField
    ftIm = fftshift(fft2(imStack(:,:,idx)));
    
    if useSupport
        paddedFT(bdpad(1):bdpad(2),bdpad(3):bdpad(4)) = maskFilt.*ftIm;
    else
        paddedFT(bdpad(1):bdpad(2),bdpad(3):bdpad(4)) = ftIm;
    end
    imStackPad(:,:,idx) = real(ifft2(ifftshift(paddedFT)));
    imStackPad(:,:,idx) = imStackPad(:,:,idx).*(imStackPad(:,:,idx)>=0);
end

recFTframe = zeros(xsize,ysize,nBrightField); 
% recFTframe: store the reconstructed FT, used to calculate the aberration

imSum = mean(imStackPad(:,:,1:nBrightField),3);
maxSig = max(imSum(:));
if autoNoiseFloor
    regImage = min(1/padfactor, maxSig/10000); %65535: 16 bit camera
end
[Y,X] = meshgrid(1:ysize*padfactor,1:xsize*padfactor);

% build wiener filters
if useNoiseFilter
    wienerPixelTol = 3;
    wienerRingWidth = 5;
    R = abs(X(bdpad(1):bdpad(2),bdpad(3):bdpad(4))-xcpad + 1i*(Y(bdpad(1):bdpad(2),bdpad(3):bdpad(4)) - ycpad));
    CTFmax = sqrt(sum(CTF(:))/pi);
    wienerMask = (R>(CTFmax + wienerPixelTol)) & (R<(CTFmax + wienerPixelTol + wienerRingWidth));
    if isa('CTF','logical')
        wienerMask = wienerMask & ~CTF;
    else
        wienerMask = wienerMask & (abs(CTF)<1e-2);
    end
end
mask2use = double(CTF>0.9);

%% reconstruction
for idx = 1:nBrightField
    realPart = log(imStackPad(:,:,idx)+regImage);
    tempFT = fftshift(fft2(realPart));
    maskOneside = (-(X-xcpad)*mycoord(idx,1) - (Y-ycpad)*mycoord(idx,2) > 0.5*norm(mycoord(idx,:)));
    temp = tempFT.*maskOneside;
    recField = exp(0.5*realPart + 1i*imag(ifft2(ifftshift(temp))));
    
    mycoord(idx,:) = round(mycoord(idx,:)); % rounded up as it's going to be used as indices
    temp = circshift(fftshift(fft2(recField))/padfactor,mycoord(idx,:));
    
    if useCTF
        zeroFreqWeight = temp(xcpad + mycoord(idx,1),ycpad + mycoord(idx,2));
        if double(CTF(xc + mycoord(idx,1),yc + mycoord(idx,2))) < 10^-2
            error('Peak corresponds to the zero-frequency cropped out using current CTF. Please increase the CTF.');
        end
        
        if useNoiseFilter
            wienerFilterTemp = maskOneside(bdpad(1)-mycoord(idx,1):bdpad(2)-mycoord(idx,1),...
                                           bdpad(3)-mycoord(idx,2):bdpad(4)-mycoord(idx,2)) & wienerMask;
            wienerFilterTemp = wienerFilterTemp & circshift(R,mycoord(idx,:))<2*CTFmax;
            noiseReg = 0.2*sum(sum(abs(temp(bdpad(1):bdpad(2),bdpad(3):bdpad(4))).*wienerFilterTemp))/sum(wienerFilterTemp(:));
            recFTframe(:,:,idx) = temp(bdpad(1):bdpad(2),bdpad(3):bdpad(4)).*mask2use./(double(CTF)+10^-3);
            recFTframe(:,:,idx) = recFTframe(:,:,idx).*abs(recFTframe(:,:,idx))./(abs(recFTframe(:,:,idx)) + noiseReg);
        else
            recFTframe(:,:,idx) = temp(bdpad(1):bdpad(2),bdpad(3):bdpad(4)).*mask2use./(double(CTF)+10^-3);
        end
        recFTframe(xc + mycoord(idx,1),yc + mycoord(idx,2),idx) = zeroFreqWeight;
    else
        recFTframe(:,:,idx) = temp(bdpad(1):bdpad(2),bdpad(3):bdpad(4));
    end
    
    if useDataIntensity
        temp = ifft2(ifftshift(recFTframe(:,:,idx)));
        recFTframe(:,:,idx) = fftshift(fft2(sqrt(imStack(:,:,idx)).*exp(1i*angle(temp))));
    end

end

if nargout == 2
    varargout{1} = mask2use;
end

end