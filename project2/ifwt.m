function image = ifwt(fwt_coeff,I)
% do a M-scale iFWT 
J=I;
[M,N]=size(fwt_coeff);
while J>0
    LL=fwt_coeff(1:M/(2^J),1:N/(2^J));
    LH=fwt_coeff(M/(2^J)+1:M/(2^(J-1)),1:N/(2^J));%-128;
    HL=fwt_coeff(1:M/(2^J),N/(2^J)+1:N/(2^(J-1)));%-128;
    HH=fwt_coeff(M/(2^J)+1:M/(2^(J-1)),N/(2^J)+1:N/(2^(J-1)));%-128;
    recon_image=IFWT(LL,HL,LH,HH);
    fwt_coeff(1:M/(2^(J-1)),1:N/(2^(J-1)))=recon_image;
    J=J-1;
end
image=fwt_coeff;

end



function[recon_image]=IFWT(LL,HL,LH,HH)

% x2n_L=sqrt(0.5)*LL-0.25*sqrt(2)*LH-0.25*sqrt(2)*circshift(LH,1,1);
% x2nplus1_L=sqrt(2)*LH+0.5*x2n_L+0.5*circshift(x2n_L,-1,1);
% nColumns = size(x2n_L,2);
% x1 = [x2n_L,x2nplus1_L]'; 
% x1 = reshape(x1(:),nColumns,[])'; 
% 
% 
% x2n_R=sqrt(0.5)*HL-0.25*sqrt(2)*HH-0.25*sqrt(2)*circshift(HH,1,1);
% x2nplus1_R=sqrt(2)*HH+0.5*x2n_R+0.5*circshift(x2n_R,-1,1);
% nColumns = size(x2n_R,2);
% x2 = [x2n_R,x2nplus1_R]'; 
% x2 = reshape(x2(:),nColumns,[])'; 
% 
% 
% x2n=sqrt(0.5)*x1-0.25*sqrt(2)*x2-0.25*sqrt(2)*circshift(x2,1,2);
% x2nplus1=sqrt(2)*x2+0.5*x2n+0.5*circshift(x2n,-1,2);
x2n_L=1*LL-0.25*2*LH-0.25*2*circshift(LH,1,1);
x2nplus1_L=2*LH+0.5*x2n_L+0.5*circshift(x2n_L,-1,1);
nColumns = size(x2n_L,2);
x1 = [x2n_L,x2nplus1_L]'; 
x1 = reshape(x1(:),nColumns,[])'; 


x2n_R=1*HL-0.25*2*HH-0.25*2*circshift(HH,1,1);
x2nplus1_R=2*HH+0.5*x2n_R+0.5*circshift(x2n_R,-1,1);
nColumns = size(x2n_R,2);
x2 = [x2n_R,x2nplus1_R]'; 
x2 = reshape(x2(:),nColumns,[])'; 


x2n=1*x1-0.25*2*x2-0.25*2*circshift(x2,1,2);
x2nplus1=2*x2+0.5*x2n+0.5*circshift(x2n,-1,2);

[nRows,nClosB] = size(x2n);
nClosA = size(x2nplus1,2);
x = zeros(nRows,nClosA+nClosB);
x(:,1:2:end) = x2n;
x(:,2:2:end) = x2nplus1;
%x1=uint8(x);
% figure
% imshow(x1)
recon_image=x;
% error_image=f-x;
% figure
% imshow(error_image)
end