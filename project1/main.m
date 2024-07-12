%% part 1 Histogram equalization
clear 
clc
close all

f=imread("images\lena512.bmp");
subplot(3,2,1);
imshow("images\lena512.bmp");
title('Original Image');
subplot(3,2,2)
histogram(f,'BinWidth',1);
title('Original Image Histogram');
axis([0 255 0 2500])

%low-contrast image
a=0.2;
b=50;
for x=1:512
    for y=1:512
g(x,y)=min(max(round(a*f(x,y))+b,0),255);
    end
end

subplot(3,2,3)
imshow(g);
title('The Low-Contrast Image');
subplot(3,2,4)
histogram(g,'BinWidth',1)
title('The Low-Contrast Image Histogram');
axis([0 255 0 10500])

% 直方图均衡化
prpixel=zeros(1,256);
for i=0:255  
   prpixel(i+1)=length(find(g==i))/(512*512);  %find(g==i) find the gray value is i in the image
end 

cdf=cumsum(prpixel);
cdf_pixel=round(255*cdf);
equalizedImage=uint8(arrayfun(@(x) cdf_pixel(x),g));

subplot(3,2,5)
imshow(equalizedImage);
title('The Equalized Image');
subplot(3,2,6)
histogram(equalizedImage,'BinWidth',1);
title('The Equalized Image Histogram');

% hold on
% figure
% j=histeq(g,255);
% subplot(221)
% imshow(j)
% subplot(222)
% histogram(j,'BinWidth',1)
% 
% subplot(223)
% imshow(equalizedImage);
% subplot(224)
% histogram(equalizedImage,'BinWidth',1)




%% part 2 Image Denoising

clear; clc; close all;

im = imread('images\lena512.bmp');

% Generate Gaussian noise
gaussian_noise = mynoisegen('gaussian', 512, 512, 0, 64);
im_gaussian = im + uint8(gaussian_noise);

% Generate salt & pepper noise
salt_pepper_noise = mynoisegen('saltpepper', 512, 512, 0.05, 0.05);
im_saltp = im;
im_saltp(salt_pepper_noise == 0) = 0;
im_saltp(salt_pepper_noise == 1) = 255;

% Plot the histograms
figure;
subplot(3, 2, 1);
imshow(im);
title('Original Image');
subplot(3, 2, 3);
imshow(im_gaussian);
title('Gaussian Noisy Image');
subplot(3, 2, 5);
imshow(im_saltp);
title('Salt & Pepper Noisy Image');

subplot(3, 2, 2);
% imhist(im);
histogram(im,'BinWidth',1);
% xlabel(r);
title('Original Image Histogram');
subplot(3, 2, 4);
% imhist(im_gaussian);
histogram(im_gaussian,'BinWidth',1);
title('Gaussian Noisy Image Histogram');
subplot(3, 2, 6);
% imhist(im_saltp);
histogram(im_saltp,'BinWidth',1);
title('Salt & Pepper Noisy Image Histogram');

% Apply 3x3 mean filter to Gaussian noisy image
h_mean = ones(3, 3) / 9;
im_gaussian_mean = conv2(double(im_gaussian), h_mean, 'same');
% imshow(uint8(im_gaussian_mean))

% Apply 3x3 mean filter to salt & pepper noisy image
im_saltp_mean = conv2(double(im_saltp), h_mean, 'same');

% Plot histograms
figure;
subplot(221);
imshow(uint8(im_gaussian_mean));
title('Gaussian Noisy Image after Mean Filter');
subplot(222);
histogram(uint8(im_gaussian_mean),'BinWidth',1);
% imhist(uint8(im_gaussian_mean));
title('Gaussian Noisy Image Histogram after Mean Filter');
subplot(223);
imshow(uint8(im_saltp_mean));
title('Salt & Pepper Noisy Image after Mean Filter');
subplot(224);
histogram(uint8(im_saltp_mean),'BinWidth',1);
% imhist(uint8(im_saltp_mean));
title('Salt & Pepper Noisy Image Histogram after Mean Filter');

im_gaussian_median = median_filter(im_gaussian,3);
im_saltp_median = median_filter(im_saltp,3);

figure;
subplot(221);
imshow(uint8(im_gaussian_median));
title('Gaussian Noisy Image after Median Filter');
subplot(222);
histogram(uint8(im_gaussian_median),'BinWidth',1);
% imhist(uint8(im_gaussian_median));
title('Gaussian Noisy Image Histogram after Median Filter');
subplot(223);
imshow(uint8(im_saltp_median));
title('Salt & Pepper Noisy Image after Median Filter');
subplot(224);
histogram(uint8(im_saltp_median),'BinWidth',1);
% imhist(uint8(im_saltp_median));
title('Salt & Pepper Noisy Image Histogram after Median Filter');



%% part 3 Frequency Domain Filtering
clear
clc
close all

% read the image and implement the blurring
lena=double(imread('images\lena512.bmp'));
h=myblurgen('gaussian',8);
blur=min(max(conv2(lena,h,"same"),0),255);
figure;
subplot(121)
imshow(uint8(lena))
title('original')
subplot(122)
imshow(uint8(blur));
title('blurred')

%sketching the spectra
LENA=log(abs(fft2(lena))+1);
BLUR=log(abs(fft2(blur))+1);
H=log(abs(fft2(h,512,512))+1);
figure;
subplot(131)
imshow(fftshift(LENA),[]);
title('spectrum of original image')
subplot(132)
imshow(fftshift(H),[]);
title('spectrum of blurring function')
subplot(133)
imshow(fftshift(BLUR),[]);
title('spectrum of out-of-focus image')

% compute the noise variance
mu=blur-double(uint8(blur));
figure;
imshow(mu,[]);
MU=mu(:);
n=var(MU);

% deblur the image and show
f_hat=deblur(double(uint8(blur)),h,n);
imshow(uint8(f_hat))