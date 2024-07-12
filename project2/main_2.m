clear 
clc
close all

%FWT
image=double(imread('harbour512x512.tif'));
fwt_coeff=fwt_M_scale(image,4);
imshow(uint8(fwt_coeff))

%BITRATE
PSNR=zeros(1,10);
bit_rate=zeros(1,10);
dd=zeros(1,10);

for l=0:9
    coeff_qtzd=midtread(fwt_coeff,2^l);
    % figure
    % imshow(coeff_qtzd);
    image2=ifwt(coeff_qtzd,4);
    figure
    imshow(uint8(image2))


    d=sum((image-image2).^2,"all")/(512*512);
    dd(l+1)=d;
    PSNR(l+1)=10*log10(255^2/d);

  
    bit_rate(l+1)=entr(coeff_qtzd);
    %bit_rate(l+1)=entropy(coeff_qtzd);
end
figure
plot(bit_rate,PSNR,'ob-');
grid on

coeff_qtzd=midtread(fwt_coeff,2^0);
image2=ifwt(coeff_qtzd,4);
e=uint8(image2-image);
figure
imshow(e)

% image2=ifwt(image,1);
% figure
% imshow(uint8(image2))
