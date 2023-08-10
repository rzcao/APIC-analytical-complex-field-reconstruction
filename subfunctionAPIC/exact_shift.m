function [shift_im] = exact_shift(im,relative_pixel,varargin)
%UNTITLED6 size_flage=1: output_matrix is of the same size of the input matrix
[xsize,ysize,~]=size(im);
sameSize = true;
isRealSpace = true;

idx = 1;
while idx <= length(varargin)
    switch varargin{idx}
        case {'FT','ft','fourier','Fourier'}
            isRealSpace = false;
            idx = idx+1;
        case {'real','Real','image','Image'}
            isRealSpace = true;
            idx = idx+1;
        case {'size','Size'}
            switch varargin{idx+1}
                case 'same'
                    sameSize = true;
                case 'padded'
                    sameSize = false;
            end
            idx = idx+2;
        otherwise
            error('Option is not supported.');
    end
            
end


if sameSize
    o_xsize=xsize;
    o_ysize=ysize;
end

if mod(xsize,2)==0
    im(xsize+1,:,:)=0;
end

if mod(ysize,2)==0
    im(:,ysize+1,:)=0;
end
relative_pixel(:,1)=rem(relative_pixel(:,1),xsize);
relative_pixel(:,2)=rem(relative_pixel(:,2),ysize);
[xsize,ysize,num]=size(im);
[Y,X]=meshgrid(1:ysize,1:xsize);
xc=floor(xsize/2+1);
yc=floor(ysize/2+1);
if isRealSpace
    yr=Y-yc;
    xr=X-xc;
else
    yr=Y-1;
    xr=X-1;
end
r_shift=sqrt(relative_pixel(:,1).^2+relative_pixel(:,2).^2);
shift_angle=zeros(num);
f_r=zeros(num,2);
for ii=1:num
    if relative_pixel(ii,2)==0
        shift_angle(ii)=pi/2.*sign(relative_pixel(ii,1));
    else
        shift_angle(ii)=atan(relative_pixel(ii,1)./relative_pixel(ii,2));
        if relative_pixel(ii,2)<0
            shift_angle(ii)=shift_angle(ii)-pi;
        end
    end
    
    if r_shift(ii)~=0
        f_r(ii,1)=xsize./r_shift(ii);
        f_r(ii,2)=ysize./r_shift(ii);
    end
    
end

if sameSize
    final_xsize=o_xsize;
    final_ysize=o_ysize;
else
    final_xsize=xsize;
    final_ysize=ysize;
end
shift_im=zeros(final_xsize,final_ysize,num);
for ii=1:num
    fr_temp=f_r(ii,:);
    if fr_temp~=0
        my_angle=shift_angle(ii);
        if isRealSpace
            ft=fftshift(fft2(im(:,:,ii)));
            ft=ft.*exp(-1i*2*pi*(xr.*sin(my_angle)/fr_temp(1)+yr.*cos(my_angle)/fr_temp(2)));
            temp=ifft2(ifftshift(ft));
        else
            ft=ifft2(ifftshift(im(:,:,ii)));
            ft=ft.*exp(1i*2*pi*(xr.*sin(my_angle)/fr_temp(1)+yr.*cos(my_angle)/fr_temp(2)));
            temp=fftshift(fft2(ft));
        end
        shift_im(:,:,ii)=temp(1:final_xsize,1:final_ysize);
    else
        shift_im(:,:,ii)=im(1:final_xsize,1:final_ysize,ii);
    end
end


end

