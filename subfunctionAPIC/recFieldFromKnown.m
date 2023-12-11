function [ftRecons,maskRecons] = recFieldFromKnown(imStack,kIllu,ftRecons,maskRecons,CTF_abe,varargin)
%recFieldFromKnown reconstruct field based known FT of the measurement
% Input: 1. imStack: meaurements from the experiment
%        2. kIllu: illumination vector, which corresponds to the shift (in 
%                  pixel) in spatial frequency domain.
%        3. ftRecons: known FT part
%        4. maskRecons: mask with 1 denotes the known FT, 0 the unknown
%        Options:
%           1. drift: when enabled (set to true), the program takes
%                     possible calibration error of the illumination vector
%                     into consideration
%           2. regularization: specify the weight of the L2 regularizer
%           3. unknown ratio: only measurements with unknown spectrum that
%                             are larger than (this ratio)*CTF will be used.
%           4. high freq threshold (threshold): If the unkonwn spectrum is 
%                       smaller than this ratio, the newly measured 
%                       spectrum will be averaged, and will be added to the
%                       final spectrum in a later time. The final spectrum 
%                       is expanded immediately when the unkonwn spectrum 
%                       is larger than this ratio.
%           5. conserve energy: treat the amplitude of the measured data as
%                               the ground truth.
%           6. timer on: calculate the remaining time to finish.
%
% The source code is licensed under GPL-3. 
%
% By Ruizhi Cao, Nov 11, 2022, Biophotonics Lab, Caltech
% Modified on Oct 31, 2023

    correctDrift    = false; % whether to take angle calibration error in to consideration
    useTimer        = false; % whether to calculate the runtime estimate
    pixelTol        = 2;     % expand the calculated unknown mask by (roughly) pixelTol pixels when drift is enabled
    highFreqTHLD    = 0.3;   % can be used for an improved reconstruction quality when working with dataset with large overlap ratio.
    valueOnHold     = false; % works with nonzero highFreqTHLD
    unknownRatio    = 0;     % control which measurement to use in the reconstruction. 
                             % If the ratio unknown part to the (reconstructed) known
                             % part is below this threshold, the corresponding measurement
                             % is not used in the reconstruction.
    userDefinedReg  = false;
    userBrightness  = false;
    autoAmpMatch    = false; % whether to correct illumination intensity differences
    marginPixel     = 5;
    replaceMag      = true;  % whether to match the energy (intensity) of the complex 
                             % field reconstruction with the actual measurement.
    
    idx = 1;
    while idx <= length(varargin)
        switch lower(varargin{idx})
            case {'drift','driftcorrection','drift correction'}
                correctDrift = varargin{idx+1};
                if idx+2 <= length(varargin) && isnumeric(varargin{idx+2})
                    pixelTol = varargin{idx+2};
                    idx = idx+3;
                else
                    idx = idx+2;
                end
            case {'reg','regularization','regularizer'}
                myreg = varargin{idx+1};
                userDefinedReg = true;
                idx = idx+2;
            case {'unknown ratio','unknownratio','minratio','min ratio'}
                unknownRatio = varargin{idx+1};
                idx = idx+2;
            case {'high freq threshold','threshold','thres'}
                highFreqTHLD = varargin{idx+1};
                if highFreqTHLD > 0.5
                    warning('The threshold should not exceed 0.5. Set to 0.5 instread.');
                    highFreqTHLD = 0.5;
                end
                idx = idx+2;
            case {'timeron','timer on','time'}
                useTimer = true;
                idx = idx+1;
            case {'brightness','intensity'}
                myBrightness = varargin{idx+1};
                userBrightness = true;
                idx = idx + 2;
            case {'intensitycorrection','intensity correction','intensity match'}
                autoAmpMatch = varargin{idx+1};
                idx = idx + 2;
            case {'conserve energy','conserveenergy','usedataintensity','dataintensity','data intensity','use data intensity'}
                replaceMag = varargin{idx+1};
                idx = idx + 2;
            otherwise
                error('Supported arguments are ''drift'', ''reg'', ''unknown ratio'', ''threshold'', ''brightness'', ''intensity correction'', ''conserve energy'', and ''timer on''.');
        end
    end
    
    if ~userDefinedReg
        % when there is no regularization factor specified, the L2 regularizer's
        % weight is chosen based on whether drift correction is enabled
        if correctDrift
            myreg  = 1; % factor for L2 norm regularizer
        else
            myreg  = 0.01;
        end
    end
    
    maxAmpCTF = max(abs(CTF_abe(:)));
    if  maxAmpCTF < 0.99
        warning('CTF is unnormalized, it will be normalized by default of the program.');
        CTF_abe = CTF_abe./maxAmpCTF;
    end
    CTF = abs(CTF_abe) > 5*10^-3;    
    
    % center and size of the measurement
    [xsize,ysize,numIm] = size(imStack);
    if numIm ~= length(kIllu(:,1))
        error('Number of image and number of illumination angles mismatch.');
    end
    
    if userBrightness && length(myBrightness) ~= numIm
        error('Number of intensity calibration points disagrees with the number of images.');
    end
    
    [Y,X] = meshgrid(1:ysize,1:xsize);
    xc = floor(xsize/2+1);
    yc = floor(ysize/2+1);
    R = abs(X-xc + 1i*(Y-yc));
    
    if marginPixel > 0.1*min(xsize,ysize)
        marginPixel = ceil(max(0.08*min(xsize,ysize),2));
    end
    bdCrop = calBoundary([xc,yc],[xsize-2*marginPixel,ysize-2*marginPixel]);
    
    % center and size of the reconstruction
    [tempX,tempY] = size(ftRecons);
    if tempX ~= tempY
        error('The known reconstruction must be a square image.');
    end
    imsizeRecons = tempX;
    xcR = floor(imsizeRecons/2+1);
    ycR = floor(imsizeRecons/2+1);
    
    if ~exist('maxCTF','var')
        maxCTF = round((sum(CTF(xc,:))-1)/2);
        if sum(CTF(xc,:)) ~= sum(CTF(:,yc))
            error('The input image is not a square image. Please use zero padding.');
        end
    end
    
    if correctDrift
       CTF_larger = (R < maxCTF+1);
    end
    
    areaCTF = sum(CTF(:));
    
    if highFreqTHLD ~= 0
        unknownTHLD = highFreqTHLD*areaCTF; % threshold in terms of the area of the unknown mask
        ftExpanded = zeros(imsizeRecons,imsizeRecons); % store temporary reconstructed spectrum (with repeats)
        maskExpanded = zeros(imsizeRecons,imsizeRecons); % number of repeats for the expanded spectrum
    end

    bd = calBoundary([xcR,ycR],[xsize,ysize]);
    
    if useTimer
        myReconstructionTimer = CalTimeRemain(numIm);
    end
    
    %% perform reconstruction
    for idx = 1:numIm
        flagSubPixel = false;
        bd2use = bd - repmat(kIllu(idx,:),[2,1]);
        if any(mod(bd2use(:),1)~= 0)
            temp = bd2use - round(bd2use);
            subpixelShift = -[temp(1),temp(3)];
            bd2use = round(bd2use);
            flagSubPixel = true;
        end
        % introduce aberration to the known part to match up with real measurement
        if flagSubPixel
            ampCTF = exact_shift(abs(CTF_abe),-subpixelShift,'real');
            aglCTF = exact_shift(phase_unwrapCG(angle(CTF_abe)),-subpixelShift,'real');
            ampCTF = ampCTF.*(ampCTF>0.12);
            ampCTF = ampCTF.*(ampCTF<=1) + (ampCTF>1);
            CTF2use = ampCTF.*exp(1i*aglCTF);
        else
            if idx == 1
                CTF2use = CTF_abe.*(abs(CTF_abe)>10^-3);
                ampCTF = abs(CTF2use);
            end
        end
        lowFT=ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)).*CTF2use;

        knownMask = maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)).*CTF;
        [unknownMask,linearArea] = calUnknownMaskFromKnownMask(knownMask,CTF);
        
        if highFreqTHLD ~= 0  &&  sum(unknownMask(:)) > unknownTHLD
            % expand the spectrum of the reconstruction
            ftRecons = ftRecons + ftExpanded./(maskExpanded + eps);
            maskRecons = maskRecons + (maskExpanded ~= 0);
            
            % reset temporary spectrum and its weight mask
            maskExpanded(:) = 0;
            ftExpanded(:) = 0;
            valueOnHold = false;
            
            % regenerate the known field and the known and unknown masks
            lowFT=ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)).*CTF2use;
            knownMask = maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)).*CTF;
            [unknownMask,linearArea] = calUnknownMaskFromKnownMask(knownMask,CTF);
        end

        if sum(unknownMask(:)) > unknownRatio*areaCTF
            fieldKnown = ifft2(ifftshift(lowFT));
            imKnown = fieldKnown.*conj(fieldKnown);

            % calculate the correlation to correct the intensity
            if autoAmpMatch
                corrF = sum(sum((imKnown(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4)) - mean2(imKnown(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4))))...
                                .*(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx) - mean2(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx)))))./ ...
                        sum(sum((imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx) - mean2(imStack(bdCrop(1):bdCrop(2),bdCrop(3):bdCrop(4),idx))).^2));
            else
                corrF = 1;
            end
            
            if userBrightness
                corrF = corrF*myBrightness(idx);
            end
            
            imReal = imStack(:,:,idx)*corrF - imKnown;
            ftImSub = fftshift(fft2(imReal))*xsize*ysize;
            confinedCTFBoxSize = maxCTF + 1; % 1 more point as residual

            % vertices of a smaller square (smaller than the full image) that contains the CTF
            boxVert = [xc - confinedCTFBoxSize,xc + confinedCTFBoxSize,...
                       yc - confinedCTFBoxSize,yc + confinedCTFBoxSize]; 

            if any(boxVert<1) || boxVert(2)>xsize || boxVert(4)>ysize
                error('There is aliasing in the captured image, please reduce the CTF size.');
            end

            if correctDrift
                kernelTol = ones(2*pixelTol+1,2*pixelTol+1);
                unknownMaskOriginal = unknownMask;
                temp = imfilter(unknownMask,kernelTol);
                temp = (temp>0.05*max(temp(:))).*CTF_larger;
                unknownMask(boxVert(1):boxVert(2),boxVert(3):boxVert(4)) = ...
                                temp(boxVert(1):boxVert(2),boxVert(3):boxVert(4));
            end

            Hreduced = calConvMtx(rot90(conj(lowFT(boxVert(1):boxVert(2),boxVert(3):boxVert(4))),2),...
                                  linearArea(xc-confinedCTFBoxSize*2:xc+confinedCTFBoxSize*2,...
                                             yc-confinedCTFBoxSize*2:yc+confinedCTFBoxSize*2),...
                                  unknownMask(boxVert(1):boxVert(2),boxVert(3):boxVert(4)));
            Htemp = Hreduced'*Hreduced;
            absMean = mean2(abs(Htemp));
            
            if correctDrift
                vecWeight = (diag(unknownMaskOriginal(unknownMask)) + 0.00001)/1.00001;
                recFTvct = (Htemp + absMean*myreg*vecWeight)\(Hreduced'*ftImSub(linearArea(:) == 1));
            else
                recFTvct = (Htemp + absMean*myreg*eye(size(Htemp,1)))\(Hreduced'*ftImSub(linearArea(:) == 1));
            end

            ftTrue = zeros(xsize,ysize);
            ftTrue(unknownMask) = recFTvct;
            if correctDrift
                unknownMask = unknownMaskOriginal;
                ftTrue = ftTrue.*unknownMaskOriginal;
            end
                        
            if replaceMag % whether to maintain the (pixel-wise) energy
                fieldTemp = exp(1i*angle(ifft2(ifftshift(ftTrue + lowFT)))).*sqrt(imStack(:,:,idx)*corrF);
                ftTemp = fftshift(fft2(fieldTemp));
                ftTrue(unknownMask) = ftTemp(unknownMask);
            end
            
            ftTrue = ftTrue.*conj(CTF2use)./(abs(CTF2use) + 0.005).*(ampCTF>0.05)*1.005; % correct aberration
            
            % stitch the high-freq spectrum into the reconstructed spectrum when necessary
            if highFreqTHLD == 0
                ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
                    ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + ftTrue;

                maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
                    maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + unknownMask.*(ampCTF>0.05);
            else
                ftExpanded(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
                    ftExpanded(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + ftTrue;

                maskExpanded(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
                    maskExpanded(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + unknownMask.*(ampCTF>0.05);
                valueOnHold = true;
            end
        end
        
        if useTimer
            myReconstructionTimer.timeRemain(idx); % calculate the estimted runtime
        end
        
    end
    
    if useTimer
        myReconstructionTimer.delete;
    end
    
    if valueOnHold
        ftRecons = ftRecons + ftExpanded./(maskExpanded + eps);
        maskRecons = maskRecons + (maskExpanded ~= 0);
    end
    
end

