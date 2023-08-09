function [shiftVct_new] = calLEDPosPerPatch(shiftVct,imsize,varargin)
%calLEDPosPerPatch calculate the k-vector of a tilted illumination in a new patch based on the known k-vector and original patch's position
% Arguments:
%   1. shiftVct: the original illumination vector
%   2. imsize: size of the image which is used when calculating the
%              illumination vector.
%   3. position: position of the patch. Can use absolute position or
%                relative position. See "8. absolute" for the details when
%                using the absolute coordinate.
%   4. wavelength: wavelength of the illumination, eg '550nm'
%   5. magnification: magnification of the system, eg 10
%   6. pixel size: pixel size of the camera,       eg '5um'
%   7. height: height of the LED array,            eg '10cm'
%   8. absolute: tell the program that the position of the patch is in the
%                absolute (pixel-wise) coordinate (i.e. The center of the
%                patch's in the image in terms of pixel). In this case,
%                one must specify "origin", which is the absolute position
%                over which the k-illumination vectors are optimized. By
%                defualt, this is set to false.
%   9. origin: the coordinate of the initial patch that is used in
%              optimizing the k-vector of the tilted illumination.
%   10. relative: tell the program that the position of the patch is the
%                 offset w.r.t to the original patch. It is chosen by
%                 defualt.
%
%   WARNING: Please note that the positive direction of the patch shift 
%            (direction of the physical shift in real world) and the LED
%            position must be aligned.
%
% By Ruizhi Cao, Biophotonics Lab, Caltech, Feb 24, 2023

nIllunimation         = numel(shiftVct)/2;

if isnumeric(imsize)
     if numel(imsize)>2
         imsize = size(imsize);
     elseif numel(imsize) == 1
         imsize = [imsize,imsize];
     end
else
    error('Please set the second input as the size of the image which is used to calculate the spatial freq.');
end
absCoord    = false; % whether to use absolute patch coordinate

nFacPix    = 1;     % defualt unit for pixel size:  um/micron
nFacH      = 10^3;  % defualt unit for height:      mm
nFacLambda = 10^-3; % defualt unit for wavelength:  nm
shiftVct_new  = zeros(nIllunimation,2);

mag    = 1;  % magnification of the system
orgPos = []; % position of the original patch
pos    = []; % position of the new patch
lambda = []; % wavelength of the illumination light

idx = 1;
while idx <= length(varargin)
    switch lower(varargin{idx})
        case {'pos','position','center'}
            pos = varargin{idx+1};
            idx = idx+2;
        case {'original','origin'}
            orgPos = varargin{idx+1};
            idx = idx+2;
        case {'wavelength','lambda'}
            temp = varargin{idx+1};
            if isa(temp,'string')
                temp = convertStringsToChars(temp);
            end
            unitFlag = (temp - '.') > 11; %% ASCII order:  . / 0 1 2 3
            unit2use = temp(unitFlag);
            if ~isempty(unit2use)
                switch lower(unit2use)
                    case {'um','micron'}
                        nFacLambda = 1;
                    case 'mm'
                        nFacLambda = 10^3;
                    case 'nm'
                        nFacLambda = 10^-3;
                    otherwise
                        error('Unit of the pixel size must be ''mm'', ''um'' or ''nm''');
                end
            end
            
            lambda = str2double(temp(~unitFlag));
            if isnan(lambda)
                error('please check the input of wavelength. Except for the unit, this should only contain digits and ''.''.');
            end
            idx = idx+2;
        case {'offset','relative','rel'}
            absCoord = false;
            idx = idx+1;
        case {'absolute','abs'}
            absCoord = true;
            idx = idx+1;
        case {'mag','magnification'}
            mag = varargin{idx+1};
            idx = idx+2;
        case {'pixelsize','pixel size','psize'}
            temp = varargin{idx+1};
            if isa(temp,'string')
                temp = convertStringsToChars(temp);
            end
            unitFlag = (temp - '.') > 11; %% ASCII order:  . / 0 1 2 3
            unit2use = temp(unitFlag);
            if ~isempty(unit2use)
                switch lower(unit2use)
                    case {'um','micron'}
                        nFacPix = 1;
                    case 'mm'
                        nFacPix = 10^3;
                    case 'nm'
                        nFacPix = 10^-3;
                    otherwise
                        error('Unit of the pixel size must be ''mm'', ''um'' or ''nm''');
                end
            end
            
            pixelsize = str2double(temp(~unitFlag));
            if isnan(pixelsize)
                error('please check the input of pixelsize. Except for the unit, this should only contain digits and ''.''.');
            end
            idx = idx+2;
        case {'height','h'}
            temp = varargin{idx+1};
            if isa(temp,'string')
                temp = convertStringsToChars(temp);
            end
            unitFlag = (temp - '.') > 11; %% ASCII order:  . / 0 1 2 3
            unit2use = temp(unitFlag);
            if ~isempty(unit2use)
                switch lower(unit2use)
                    case {'um','micron'}
                        nFacH = 1;
                    case 'mm'
                        nFacH = 10^3;
                    case 'nm'
                        nFacH = 10^-3;
                    case 'cm'
                        nFacH = 10^4;
                    otherwise
                        error('Unit of the pixel size must be ''cm'', ''mm'', ''um'' or ''nm''');
                end
            end
            
            h = str2double(temp(~unitFlag));
            if isnan(h)
                error('please check the input of height. Except for the unit, this should only contain digits and ''.''.');
            end
            idx = idx+2;
            
        otherwise
            error('Supported arguments are ''position'', ''pixel size'', ''wavelength'', ''magnification'', ''height'', and ''absolute'' (''relative'').');
    end
end

%% check if all required arguments exist
if ~exist('h','var') || ~exist('pixelsize','var')
    error('To calculate the shift with respect to the location of the path, the pixel size and height of the LED array are needed.');
end

if isempty(pos) || (isempty(orgPos) && absCoord)
    error(['The position is not well defined. Check if the new position is given to the program.',...
        char(13),'If you use ''absolute'' coordinate, check the original position as well.']);
end

if isempty(lambda)
    error('Please specify the wavelength using option ''wavelength''.');
end

%% convert the position difference of the patch to relative pixel shift
if absCoord
    pos = pos - orgPos;
end

%% calculate angle shift
pixelsize = pixelsize*nFacPix/mag; % convert pixelsize to um
h         = h*nFacH;           % convert height to um
lambda    = lambda*nFacLambda; % convert height to um

if lambda > 2 || lambda < 0.3
    error('The wavelength is either smaller than 300 nm or larger than 2 micron.')
end

%% calculate the angle of the LEDs
xc = floor(imsize(1)/2+1);
yc = floor(imsize(2)/2+1);

resFreqx = 1/(2*pixelsize)/(xc-1);
resFreqy = 1/(2*pixelsize)/(yc-1);

freqMax = 1/lambda;

unitVctXY = shiftVct(:,1)*resFreqx + 1i*shiftVct(:,2)*resFreqy;
freqAbs = abs(unitVctXY);
unitVctXY = unitVctXY./freqAbs;

if max(freqAbs) > freqMax
    error('The calculated maximal spatial frequency is larger than the frequency of the light. Please check the input.');
end

compDis = tan(asin(freqAbs/freqMax)).*unitVctXY*h;
disX = real(compDis);
disY = imag(compDis);

disX_new = disX - pixelsize*pos(1);
disY_new = disY - pixelsize*pos(2);

%% convert back to illumination k-vector
unitVctXY_new = disX_new + 1i*disY_new;
disAbs_new = abs(unitVctXY_new);
unitVctXY_new = unitVctXY_new./disAbs_new;

freq_new = sin(atan(disAbs_new/h))*freqMax;
compFreq = freq_new.*unitVctXY_new;

shiftVct_new(:,1) = real(compFreq)/resFreqx;
shiftVct_new(:,2) = imag(compFreq)/resFreqy;
end

