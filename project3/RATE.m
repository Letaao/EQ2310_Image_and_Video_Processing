function [r_intra]=RATE(input,stepsize)
[height,width]=size(input);

r_intra=zeros(height/16,width/16)+1/256;
m=1:16;

for row=1:height/16
    for column=1:width/16
qian=dct16(input(16*(row-1)+m,16*(column-1)+m));
hou=midtread(qian,stepsize);
r_intra(row,column)=r_intra(row,column)+entr(hou);
    end
end

end
