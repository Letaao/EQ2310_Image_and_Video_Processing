function fwt_coeff = fwt_M_scale(image,I)
% do a M-scale FWT 
J=I;
[M,N]=size(image);
while J>0
    A=image(1:M/(2^(I-J)),1:N/(2^(I-J)));
    A=FWT(A);
    image(1:M/(2^(I-J)),1:N/(2^(I-J)))=A;
    J=J-1;
end
fwt_coeff=image;
end
% function [output_image,low11,low12,high21,high22]=FWT(f)
% %列
% x2n_column=f(:,1:2:end);
% x2nplus1_column=f(:,2:2:end);
% 
% high_col=1/sqrt(2)*(x2nplus1_column - 0.5*x2n_column - 0.5*circshift(x2n_column,-1,2));
% low_col=sqrt(2)*(x2n_column + 0.25*high_col + 0.25*circshift(high_col,1,2));
% 
% %行
% x2n_row1=low_col(1:2:end,:);
% x2nplus1_row1=low_col(2:2:end,:);
% 
% high21=1/sqrt(2)*(x2nplus1_row1 - 0.5*x2n_row1 - 0.5*circshift(x2n_row1,-1,1));
% low11=sqrt(2)*(x2n_row1 + 0.25*high21 + 0.25*circshift(high21,1,1));
% 
% x2n_row2=high_col(1:2:end,:);
% x2nplus1_row2=high_col(2:2:end,:);
% high22=1/sqrt(2)*(x2nplus1_row2 - 0.5*x2n_row2 - 0.5*circshift(x2n_row2,-1,1));
% low12=sqrt(2)*(x2n_row2 + 0.25*high22 + 0.25*circshift(high22,1,1));
% 
% 
% 
% out1=[low11,low12];
% out2=[high21,high22];
% output_image=[out1;out2];
% end
function [output_image]=FWT(f)
%列
x2n_column=f(:,1:2:end);
x2nplus1_column=f(:,2:2:end);

high_col=1/2*(x2nplus1_column - 0.5*x2n_column - 0.5*circshift(x2n_column,-1,2));
low_col=1*(x2n_column + 0.25*2*high_col + 0.25*2*circshift(high_col,1,2));

%行
x2n_row1=low_col(1:2:end,:);
x2nplus1_row1=low_col(2:2:end,:);

high21=1/2*(x2nplus1_row1 - 0.5*x2n_row1 - 0.5*circshift(x2n_row1,-1,1));
low11=1*(x2n_row1 + 0.25*2*high21 + 0.25*2*circshift(high21,1,1));

x2n_row2=high_col(1:2:end,:);
x2nplus1_row2=high_col(2:2:end,:);
high22=1/2*(x2nplus1_row2 - 0.5*x2n_row2 - 0.5*circshift(x2n_row2,-1,1));
low12=1*(x2n_row2 + 0.25*2*high22 + 0.25*2*circshift(high22,1,1));



out1=[low11,low12];
out2=[high21,high22];
output_image=[out1;out2];
end