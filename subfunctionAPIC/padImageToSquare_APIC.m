function [I_lowPad,xsize,ysize,xc,yc,freqXY_calib,xCropPad,yCropPad] = padImageToSquare_APIC(I_low,xsize,ysize,xc,yc,freqXY_calib,na_calib)
% pad measurements to square if the input images are not
    warning('For a rectangular image, we scale the pixel related k-vector. Also, to improve the efficiency, consider to use a squre patch instead.');
    maxsize = max([xsize,ysize]);
    numIm = size(I_low,3);
    
    ratioXY = [freqXY_calib(:,1) - yc, freqXY_calib(:,2) - xc]./na_calib;
    
    I_lowPad = zeros(maxsize,maxsize,numIm);
    for idx = 1:numIm
        I_lowPad(:,:,idx) = padarray(I_low(:,:,idx),[maxsize-xsize,maxsize-ysize],...
                                     mean2(I_low(:,:,idx)),'post');
    end
    xCropPad = xsize;
    yCropPad = ysize;
    
    
    if abs(mean(ratioXY(:,1))./mean(ratioXY(:,2)) - yCropPad/xCropPad)>0.01
        error(['For a rectangle image, the illumination k-vector (in pixels) does not match its illumination NA.',...
               ' Please make sure the k-vectors are calibrated with images that are of the same dimension of the reconstruction.',...
               ' If your k-vectors are matched with the dimension of some sqaure image, please perform padding manually and then ',...
               'feed into APIC. In padding, pleae pad an image using the mean intensity of that particular image.']);
    end
    freqXY_calib(:,1) = freqXY_calib(:,1)*maxsize/yCropPad; % note the yCropPad is used in the rescaling
    freqXY_calib(:,2) = freqXY_calib(:,2)*maxsize/xCropPad;
    xsize = maxsize;
    ysize = maxsize;
    xc = floor(xsize/2+1); 
    yc = floor(ysize/2+1);
end

