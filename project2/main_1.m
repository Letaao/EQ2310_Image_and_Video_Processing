%% bit-rate estimation

clear 
clc
close all

M=8;

A=zeros(M, M);
for i=0:M-1
    for k=0:M-1
        if i==0
            alpha=sqrt(1/M);
        else
            alpha=sqrt(2/M);
        end
        A(i+1,k+1)=alpha*cos(((2*k+1)*i*pi)/(2*M));
    end
end

% Display the DCT matrix A
disp('DCT Matrix A:');
disp(A);


image=double(imread("peppers512x512.tif"));
Blocks=zeros(8,8,64*64);
for i=1:64
    for k=1:64
        Blocks(:,:,64*(k-1)+i)=image((k-1)*8+1:(k-1)*8+8,(i-1)*8+1:(i-1)*8+8);
    end
end
DCT_co=zeros(8,8,64*64);
for i=1:4096
    DCT_co(:,:,i)=A*Blocks(:,:,i)*A';
end

% DCT_co_qtzd=midtread(DCT_co,1);
% dd=sum(sum(sum((DCT_co-DCT_co_qtzd).^2)))/(512*512)
% 
% 
% 
% Blocks_recstr=zeros(8,8,64*64);
% for i=1:4096
%     Blocks_recstr(:,:,i)=A'*DCT_co_qtzd(:,:,i)*A;
% end
% 
% d=sum(sum(sum((Blocks-Blocks_recstr).^2)))/(512*512)
% 
% 
% PSNR=10*log10(255^2/d)

% for i=1:4096
%     en(i)=entr(DCT_co_qtzd(:,:,i));
% end
% bit_rate=mean(en)

Blocks_recstr=zeros(8,8,64*64);
en=zeros(1,4096);
PSNR=zeros(1,10);
bit_rate=zeros(1,10);
for l=0:9
    DCT_co_qtzd=midtread(DCT_co,2^l);
    for i=1:4096
        Blocks_recstr(:,:,i)=A'*DCT_co_qtzd(:,:,i)*A;
        en(i)=entr(DCT_co_qtzd(:,:,i));
    end
    d=sum(sum(sum((Blocks-Blocks_recstr).^2)))/(512*512);
    PSNR(l+1)=10*log10(255^2/d);
    bit_rate(l+1)=mean(en);
end
plot(bit_rate,PSNR,'ob-');
grid on





