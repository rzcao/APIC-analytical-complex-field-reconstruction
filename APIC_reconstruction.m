% This is the main reconstrution code for APIC, the Angular Ptychograhic
% Imaging with Closed-form method. APIC uses NA-matching and darkfield
% measurements and reconstructs a high-resolution, aberration-free complex
% sample field.
% 
% A tutorial of this code and the dataset can be found at:
% https://github.com/rzcao/APIC-analytical-complex-field-reconstruction
%
% The source code is licensed under GPL-3. (Due to conflict of interest of
% our provisional patent, we will release the editable source code of the
% aberration correction and the darkfield related complex field
% reconstruction algorithm after Oct 30, 2023.)
% 
% By Ruizhi Cao, Biophotonics Lab, Caltech
% Modified on Aug 7, 2023

%% check dependancies 
if ~any(any(contains(struct2cell(ver), 'Image Processing Toolbox')))
    error('Image processing toolbox is required.');
end
addpath(genpath('subfunctionAPIC')); % add subfunctions to matlab's search path

%% options for the reconstruction
saveResult          = false;    % whether to save the reconstruction results.
matchingTol         = 0.03;     % maximal tolerence in the deviation between the illumination NA of the NA-matching measurement and the objective NA
useAbeCorrection    = true;     % whether to use aberration correction
enableROI           = true;     % whether to use ROI in the reconstruction, this is the reconstruction size used in the code
ROILength           = 256;      % define the ROI of the reconstruction
ROIcenter           = 'auto';   % define the center of ROI. Example: ROIcenter = [256,256]; ROIcenter = 'auto';
paddingHighRes      = 3;        % define the upsampling ratio for the final high-res image
visualizeNAmatchingMeas = false; % whether to visualize the result using only NA-matching meaurements
%% load data
folderName      = 'reducedData'; % note this is case sensitive
fileNameKeyword = 'Siemens';     % e.g.: Siemens HEpath Thyroid
additionalKeyword = 'defocus0';  % this is only needed when there are multiple files contain the keyword in the above line

temp = dir(folderName);
if isempty(temp)
    error(['No folder with name ''',folderName,''' under current directory.']);
end
dataFolderName = temp.name;
temp = dir([folderName,filesep,dataFolderName,filesep,'*',fileNameKeyword,'*.mat']);
if isempty(temp)
    error(['No .mat file with name ''', fileNameKeyword, ''' found in folder ',dataFolderName]);
elseif length(temp)>1
    temp = dir([folderName,filesep,dataFolderName,filesep,'*',fileNameKeyword,'*',additionalKeyword,'*.mat']);
    if length(temp) > 1
        error('Multiple raw data files found in the folder. Consider to specify the full name or check the additional keyword.');
    end
end
fileName = [folderName,filesep,dataFolderName,filesep,temp.name]; 
load(fileName);

%% Select measurement whose illumination NA matches up with the objective NA
NA_diff = abs(sqrt(na_calib(:,1).^2+na_calib(:,2).^2)-na_cal); % calculate the absolute difference between illumination NA and objective NA 
slt_idx = find(NA_diff<matchingTol)'; % select the NA matching measurements based on the tolerence threshold
nNAmatching = length(slt_idx);
[xsize,ysize] = size(I_low(:,:,1));
xc = floor(xsize/2+1); 
yc = floor(ysize/2+1);

if ROILength>xsize || ROILength>ysize
    error(['ROI length cannot exceed ',num2str(min(xsize,ysize))]);
end

% Get the calibrated illumination angles for NA-matching measurements
x_illumination = freqXY_calib(slt_idx,2); % note the second column is assigned to x. 
y_illunimation = freqXY_calib(slt_idx,1); % note the first column is assigned to y.
NA_pixel = na_rp_cal; % calibrated maximum spatial freq in FT space
disp(['Number of NA-matching measurments found: ',num2str(nNAmatching)]);
%% select dark field measurement
NA_diff = abs(na_calib(:,1)+1i*na_calib(:,2)) - 0.008 - na_cal;
slt_idxDF = find(NA_diff>0)';

% LED illumination angle, darkfield measurements
x_illumination = [x_illumination;freqXY_calib(slt_idxDF,2)]; % note the second column is assigned to x. 
y_illunimation = [y_illunimation;freqXY_calib(slt_idxDF,1)]; % note the first column is assigned to y.

% change center to where the zero frequency is
x_illumination = x_illumination-xc;
y_illunimation = y_illunimation-yc;
if enableROI
    x_illumination = x_illumination*ROILength/xsize;
    y_illunimation = y_illunimation*ROILength/ysize;
    NA_pixel = NA_pixel*ROILength/xsize;
    na_rp_cal = na_rp_cal*ROILength/xsize;
    
    if isnumeric(ROIcenter)
        bdROI = calBoundary(ROIcenter,ROILength); 
    elseif strcmpi(ROIcenter,'auto')
        bdROI = calBoundary([xc,yc],ROILength); % ROI locates in the center of the image
    else
        error('ROIcenter should be a 1-by-2 vector or ''auto''.');
    end
    if any(bdROI(:)<1) || bdROI(2)>xsize || bdROI(4)>ysize
        error('ROI exceeds the boundary. Please check ROI''s center and length');
    end
    xsize = ROILength;
    ysize = ROILength;
else
    bdROI = [1,1;xsize,ysize]; % by default, use the maximum ROI
end
I_kohler = sum(I_low(bdROI(1):bdROI(2),bdROI(3):bdROI(4),:),3);
I = I_low(bdROI(1):bdROI(2),bdROI(3):bdROI(4),[slt_idx,slt_idxDF]);

clear I_low
%% preparing for reconstruction
imStack = zeros(size(I)); % allocate space for the filtered measurement.
                          % It is ordered based on the associated angles of
                          % each measurements. The first nNAmatching images
                          % are the NA-matching measurements.

% order measurement under NA-matching angle illumination
[theta,pupilR] = cart2pol(x_illumination(1:nNAmatching),y_illunimation(1:nNAmatching));
[~,idxMap] = sort(theta);

% order dark field measurements
[~,pupilR_DF] = cart2pol(x_illumination(nNAmatching+1:end),y_illunimation(nNAmatching+1:end));
[~,idxMapDF] = sort(pupilR_DF);
idxMapDF = idxMapDF + nNAmatching;
idxAll = [idxMap;idxMapDF];

% load the optimized CTF size and position
enlargeF = 4;
[Y,X] = meshgrid(1:ysize*enlargeF,1:xsize*enlargeF);
xc = floor(xsize*enlargeF/2+1);
yc = floor(ysize*enlargeF/2+1);
R_enlarge = abs(X-xc + 1i*(Y-yc));

k_illu = [x_illumination(idxAll),y_illunimation(idxAll)]; % ordered illumination k-vector
[Y,X] = meshgrid(1:ysize,1:xsize);
xc = floor(xsize/2+1);
yc = floor(ysize/2+1);
R = abs(X-xc + 1i*(Y-yc));

% calculate maximum spatial frequency based on the measurement.
pupilRadius = max([na_rp_cal, max(pupilR), vecnorm(fix([x_illumination(1:nNAmatching),y_illunimation(1:nNAmatching)]).')]);

CTF = imresize(double(R_enlarge<pupilRadius*enlargeF),[xsize,ysize],'bilinear');
CTF = max(circshift(rot90(CTF,2),[mod(xsize,2),mod(ysize,2)]),CTF);
% generate a binary mask for the (Fourier) support of acquired intensity images
binaryMask = (R <= 2*pupilRadius);


% taper the edge to avoid ringing effect
edgeMask = zeros(xsize,ysize);
pixelEdge = 3;
edgeMask(1:pixelEdge,:) = 1;
edgeMask(end-pixelEdge+1:end,:) = 1;
edgeMask(:,1:pixelEdge) = 1;
edgeMask(:,end-pixelEdge+1:end) = 1;
edgeMask = imgaussfilt(edgeMask,5);
maxEdge = max(edgeMask(:));
edgeMask = (maxEdge - edgeMask)/maxEdge;

noiseLevel = zeros(1,length(idxAll));
for idx = 1:length(idxAll)
    ftTemp = fftshift(fft2(I(:,:,idxAll(idx))));
    noiseLevel(idx) =  max([eps mean2(abs(ftTemp(~binaryMask)))]);
    ftTemp = ftTemp.*abs(ftTemp)./(abs(ftTemp) + noiseLevel(idx));
    if idx>nNAmatching
        imStack(:,:,idx) = ifft2(ifftshift(ftTemp.*binaryMask)).*edgeMask;
    else
        imStack(:,:,idx) = ifft2(ifftshift(ftTemp.*binaryMask));
    end
end
clear I;

% allocate space for the final high-resolution image
imsizeRecons = paddingHighRes*xsize;
ftRecons = zeros(imsizeRecons,imsizeRecons);
maskRecons = zeros(imsizeRecons,imsizeRecons);

% center and grid of the final high-res image
xcR = floor(imsizeRecons/2+1);
ycR = floor(imsizeRecons/2+1);
[YR,XR] = meshgrid(1:imsizeRecons,1:imsizeRecons);
R_recons = abs((XR-xcR) + 1i*(YR-ycR));
clear YR XR;

%% field reconstruction of NA-matching angle measurements and aberration extraction
timerAPIC = tic;
[recFTframe,mask2use] = recFieldKK(imStack(:,:,1:nNAmatching),k_illu(1:nNAmatching,:),...
    'CTF',CTF,'pad',4,'norm',true,'wiener',false); % recFTframe: reconstructed field of NA-matching measurements


[CTF_abe,zernikeCoeff] = findAbeFromOverlap(recFTframe,k_illu(1:nNAmatching,:),CTF,'weighted',true);%  find aberration 

% visualization
cmapAberration = cNeoAlbedo; % load colormap for aberration
figure;imagesc(angle(CTF_abe.*(abs(CTF_abe)>10^-3)));axis off;axis image;
title('Reconstructed pupil, APIC');caxis([-pi pi]);colorbar;colormap(cmapAberration);


% correct and stitch the reconstructed spectrums using NA-matching angle measurements
bd = calBoundary([xcR,ycR],[xsize,ysize]);
normMask = zeros(size(maskRecons));
maskRecons(xcR,ycR) = 1;
for idx = 1:nNAmatching
    bd2use = bd - repmat(round(k_illu(idx,:)),[2,1]);
    maskOneside = (-(X-xc - k_illu(idx,1))*k_illu(idx,1)...
                   - (Y-yc - k_illu(idx,2))*k_illu(idx,2) > -0.5*norm(k_illu(idx,:)));
    mask2useNew = mask2use & maskOneside; % mask to filter the reconsturcted spectrums
    
    unknownMask = (1-maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4))).*mask2use;
    maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) =...
        maskRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + unknownMask;
    offsetPhase = angle(CTF_abe(xc+round(k_illu(idx,1)),yc+round(k_illu(idx,2))));

    normMask(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) =...
        normMask(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + mask2useNew;
    if useAbeCorrection
        ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
            ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + ...
            recFTframe(:,:,idx).*conj(CTF_abe)*exp(1i*offsetPhase)./(abs(CTF_abe) + 10^-3);
    else
        ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) = ...
            ftRecons(bd2use(1):bd2use(2),bd2use(3):bd2use(4)) + ...
            recFTframe(:,:,idx).*mask2use./(mask2use + 10^-3);
    end
end
normMask(xcR,ycR) = nNAmatching; 
ftRecons = ftRecons.*(normMask>0.5)./(normMask + 10^-5).*maskRecons;
tempKKmask = maskRecons;
himMatching = ifft2(ifftshift(ftRecons)); % reconstructed complex field using NA-matching measurements

tempMask = imresize(edgeMask,size(ftRecons),'bilinear');
himMatching = himMatching.*sqrt(tempMask);
ftRecons = fftshift(fft2(himMatching));

%% visualization of reconstructed field with NA-matching angle illumination
if visualizeNAmatchingMeas
    cmapPhase = cDivVlag; % load colormap for phase
    figure('position',[35,5,1280,560]);
    subplot(121),imagesc(abs(himMatching));caxis([0 inf]);colorbar;colormap(gca,'gray');
    axis image;axis off;
    title('Amplitude, using NA-matching angle measurements, aberration corrected');
    subplot(122),imagesc(angle(himMatching));axis image;axis off;
    title('Phase, using NA-matching angle measurements, aberration corrected');
    colormap(gca,cmapPhase);colorbar;caxis([-pi pi]);
end
%% reconstruction using dark field measurements
if ~useAbeCorrection
    CTF_abe = CTF;
end

[ftRecons,maskRecons] = recFieldFromKnown(imStack(:,:,nNAmatching+1:end),k_illu(nNAmatching+1:end,:),...
                                          ftRecons,maskRecons,CTF_abe,'drift',true,'threshold',0.3,...
                                          'reg',0.01,'intensity correction', true...
                                          ,'use data intensity',true,'timer on');

himAPIC = ifft2(ifftshift(ftRecons)); % high-resolution image of APIC
runtimeAPIC = toc(timerAPIC); % reconstruction time using APIC

%% visualization of full reconstruction of APIC
imsizeRecons = size(himAPIC,1); % size of the final reconstruction
[xsize,ysize] = size(edgeMask);
[~,edgePixel] = max(edgeMask(floor(xsize/2+1),:));

edgePixeltemp = round(edgePixel/2);
bdDisp = [edgePixeltemp*imsizeRecons/xsize,imsizeRecons - edgePixeltemp*imsizeRecons/xsize+1]; % coordinates for display purpose

cmapPhase = cDivVlag; % load colormap for phase
figure('position',[35,5,1280,560]);
subplot(121),imagesc(abs(himAPIC(bdDisp(1):bdDisp(2),bdDisp(1):bdDisp(2))));colorbar;
caxis([0 inf]);axis image;axis off;colormap(gca,'gray');
title('Amplitude, APIC, w/ aberration correction');
subplot(122),imagesc(angle(himAPIC(bdDisp(1):bdDisp(2),bdDisp(1):bdDisp(2))));
axis image;axis off;colormap(gca,cmapPhase);colorbar;caxis([-pi pi]);
title('Phase, APIC, w/ aberration correction');

if exist('I_kohler','var') && any(size(I_kohler) ~= size(himAPIC))
    I_kohler = imresize(I_kohler,size(himAPIC));
end

amp_kohler = sqrt(I_kohler);
figure;imagesc(amp_kohler(bdDisp(1):bdDisp(2),bdDisp(1):bdDisp(2)));colorbar;
axis image;axis off;
caxis([0 inf]);colormap('gray');title('Amplitude, Kohler illumination');

%% save result
if saveResult
    CTF_abe = CTF_abe.*(abs(CTF_abe)>10^-3);
    [~,name2use,~] = fileparts(fileName);
    saveDir = 'Results';
    if ~exist(saveDir, 'dir')
       mkdir(saveDir)
    end
    save([saveDir filesep name2use '_APIC_',num2str(ROILength),'.mat'],'himMatching','himAPIC','CTF_abe','edgeMask','edgePixel','zernikeCoeff');
end
