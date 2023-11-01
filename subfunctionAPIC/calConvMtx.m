function [convMatrix] = calConvMtx(ftKnown,maskMeas,maskUnknown)
%Genrate the convolution matrix (in a faster way, and works with sparse entires/measurements)
% The first argument should be the Fourier transform of the known part,
% the second (logic array) should indicate measurement in the same domain
% as the first input, and the third one should indicate the unknown entries
% and, again, in the same domain as the first input.
%
% The source code is licensed under GPL-3. 
%
% By Ruizhi Cao, Sep 23, 2022, Biophotonics Lab, Caltech
% Modified on June 1, 2023

    [xsize,ysize] = size(ftKnown);
    if xsize ~= ysize
        error('Function ''calConvMtx'' only works with square images/matrices.');
    end
    n = xsize;
    nConv = 2*n-1;
    
    temp = zeros(nConv,n);
    temp(1:n,1:n) = ftKnown;
    temp2 = zeros(1,n);
    cellMat2 = zeros([nConv,n,n+1]);
    for idx = 1:n
        temp2(1) = temp(1,idx);
        cellMat2(:,:,idx) = toeplitz(temp(:,idx),temp2);
    end
    nonzeroMtx = (cellMat2 ~= 0);
    nonzeroIdx = any(any(nonzeroMtx,1),2);
    
    convMatrix = zeros(sum(maskMeas(:)),sum(maskUnknown(:)));
    idxMat = toeplitz([1:n,ones(1,n-1)*n+1],[1,ones(1,n-1)*n+1]);
    
    if isvector(maskMeas)
        maskMeas = reshape(maskMeas(:),[nConv,nConv]);
        % its columns indicates the entries that ultimately generate the measurement
    end
    
    if isvector(maskUnknown)
        try
            maskUnknown = reshape(maskUnknown(:),[n,n]);
            % its columns indicates the entries that is unknown
        catch
            error('Function ''calConvMtx'' assumes the two images in the convolution are of the same size.');
        end
    end
    
    num2fillLeft = sum(maskMeas,1);
    num2fillRight = sum(maskUnknown,1);
    nonEmptyLeft = any(maskMeas,1);
    nonEmptyRight = any(maskUnknown,1);
    
    idxPrevLeft = 0;
    for idx1 = 1:nConv
        idxPrevRight = 0;
        if nonEmptyLeft(idx1)
            for idx2 = 1:n
                if nonEmptyRight(idx2) && nonzeroIdx(idxMat(idx1,idx2))
                    convMatrix(idxPrevLeft +1 : idxPrevLeft + num2fillLeft(idx1),...
                               idxPrevRight+1 : idxPrevRight+ num2fillRight(idx2)) = cellMat2(maskMeas(:,idx1),maskUnknown(:,idx2),idxMat(idx1,idx2));
                end
                idxPrevRight = idxPrevRight+num2fillRight(idx2);
            end
        end
        idxPrevLeft = idxPrevLeft + num2fillLeft(idx1);
    end
    
end

