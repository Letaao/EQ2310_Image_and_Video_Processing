function [d_intra]=intra_mode_d(input,stepsize)
[height,width]=size(input);
d_intra=zeros(height/16,width/16);
m=1:16;
%intra mode
for row=1:height/16
    for column=1:width/16

qian2=dct16(input(16*(row-1)+m,16*(column-1)+m));
hou2=midtread(qian2,stepsize);
d_intra(row,column)=immse(hou2,qian2);
    end
end
end


