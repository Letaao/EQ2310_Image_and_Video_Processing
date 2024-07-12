clc
clear
%mother-daughter_qcif/mother-daughter_qcif.yuv /
% foreman_qcif\foreman_qcif.yuv
height=144;
width=176;
numfram=50;
size=16;
Y = yuv_import_y('mother-daughter_qcif/mother-daughter_qcif.yuv',[width height],numfram);
% mov = yuv2mov(Y);
% movie(mov);
numblock=(width/16)*(height/16);
stepsize=[8,16,32,64];

f=zeros(height,width,numfram);
for i=1:50
    f(:,:,i)=Y{i,1};
end

RR=zeros(1,4);
DD=zeros(1,4);

for q=3:6
rr=zeros(height/16,width/16,numfram);
dd=zeros(height/16,width/16,numfram);
m=1:16;
for i= 1:numfram  
    for row=1:height/16
        for column=1:width/16    
T=dct16(f(16*(row-1)+m,16*(column-1)+m,i));
Q=midtread(T,2^q);

rr(row,column,i)=entr(Q);
dd(row,column,i)=immse(T,Q);
        end
    end
end

RR(1,q-2)=sum(rr,"all");
DD(1,q-2)=sum(dd,"all");

end
r=RR/50/1000*30*256;
P=10*log10(255*255./(DD/50/99));
plot(r,P,'go-');
xlabel('rate(kbps)');
ylabel('PSNR(dB)')
grid on