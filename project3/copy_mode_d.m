function [d_copy]=copy_mode_d(input0,input1,stepsize)
[height,width]=size(input1);
d_copy=zeros(height/16,width/16);
%copy mode
m=1:16;
for row=1:height/16
    for column=1:width/16
    

f_dct=dct16(input0(16*(row-1)+m,16*(column-1)+m));
f_q=midtread(f_dct,stepsize);
f_idct=idct16(f_q);


d_copy(row,column)=immse(f_idct,input1(16*(row-1)+m,16*(column-1)+m));
    end
end

end
% %加第一帧
% for row=1:height/16copy_mode
%         for column=1:width/16
% hou1=dct16(input(16*(row-1)+m,16*(column-1)+m,1));
% qq1=dct16(input(16*(row-1)+m,16*(column-1)+m,1));
% qian1=midtread(qq1,stepsize);
% d_copy(row,column,1)=immse(hou1,qian1);
%         end
% end
% %

%intra mode
% for i= 1:frame  
%     for row=1:height/16
%         for column=1:width/16
% 
% qian2=dct16(input(16*(row-1)+m,16*(column-1)+m,i));
% hou2=midtread(qian2,stepsize);
% d_intra(row,column,i)=immse(hou2,qian2);
%         end
%     end
% end
% end
