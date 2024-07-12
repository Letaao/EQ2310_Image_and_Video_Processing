clc
clear

height=144;
width=176;
numfram=50;
size=16;
Y = yuv_import_y('mother-daughter_qcif/mother-daughter_qcif.yuv',[width height],numfram);
% numblock=(width/16)*(height/16);
% stepsize=[8,16,32,64];

f=zeros(height,width,numfram);
for i=1:50
    f(:,:,i)=Y{i,1};
end

R=zeros(height/16,width/16,numfram);
D=zeros(height/16,width/16,numfram);
R_plot=zeros(1,4);
D_plot=zeros(1,4);
m=1:16;
num=zeros(2,4);

for q=3:6
%the first frame
r1=RATE(f(:,:,1),2^q);
d1=intra_mode_d(f(:,:,1),2^q);
R(:,:,1)=r1;
D(:,:,1)=d1;

for j=2:50
    last_frame=f(:,:,j-1);
    current_frame=f(:,:,j);

    d_copy=copy_mode_d(last_frame,current_frame,2^q);
    d_intra=intra_mode_d(current_frame,2^q);

    r_copy=zeros(height/16,width/16)+1/256;
    r_intra=RATE(current_frame,2^q);

    r_final=zeros(height/16,width/16);
    d_final=zeros(height/16,width/16);

    J_copy=Lagrangian_cost(d_copy,r_copy,2^q);
    J_intra=Lagrangian_cost(d_intra,r_intra,2^q);

    selection=J_intra-J_copy;

    for row=1:height/16
        for column=1:width/16
            if(selection(row,column)>=0)
                r_final(row,column)=r_copy(row,column);
                d_final(row,column)=d_copy(row,column);

                f_dct=dct16(last_frame(16*(row-1)+m,16*(column-1)+m));
                f_q=midtread(f_dct,2^q);
                f_idct=idct16(f_q);
                current_frame(16*(row-1)+m,16*(column-1)+m)=f_idct;
                f(:,:,j)=current_frame;

            else
                r_final(row,column)=r_intra(row,column);
                d_final(row,column)=d_intra(row,column);
            end
        end
        R(:,:,j)=r_final;
        D(:,:,j)=d_final;
    end
end
num(1,q-2)=length(find(R==0.003906250000000));
num(2,q-2)=99*50-num(1,q-2);

R_plot(1,q-2)=sum(R,"all");
D_plot(1,q-2)=sum(D,"all");

end

R_plot=R_plot/50/1000*30*256;
PSNR=10*log10(255*255./(D_plot/50/99));
plot(R_plot,PSNR,'y*-');
grid on
% figure
% 
% bar(num');
% xticklabels({'2^3', '2^4', '2^5', '2^6'});
% xlabel('quantization stepsize');
% ylabel('num of blocks')
% legend('copy','intra');

