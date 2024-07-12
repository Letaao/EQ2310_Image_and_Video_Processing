clear
clc
%close all

Y = yuv_import_y('mother-daughter_qcif/mother-daughter_qcif.yuv',[176 144],50);

M=176/16;
N=144/16;

rec_video = cell(50,1);
D = zeros(1,4);
R = zeros(1,4);

cnt_intra = zeros(1,4);
cnt_copy = zeros(1,4);
cnt_inter = zeros(1,4);

for l = 3:6
    % for the first frame, only intra mode can be chosen.
    cur_frame = Y{1};
    [d,r,rec_image] = intra_mode(cur_frame, 2^l);
    rec_video{1} = rec_image;

    D(l-2) = sum(d,"all");
    R(l-2) = sum(r,"all");
    cnt_intra(l-2) = 99;

    for frame = 2:50

        cur_frame = Y{frame};
        last_rec_frame = rec_video{frame-1};

        [d_intra,r_intra,rec_image_intra] = intra_mode(cur_frame, 2^l);
        [d_copy,r_copy,rec_image_copy] = copy_mode(cur_frame, last_rec_frame);
        [d_inter,r_inter,rec_image_inter] = inter_mode(cur_frame, last_rec_frame, 2^l);
        [d,r,rec_image,choice] = mode_choice(d_intra,r_intra,rec_image_intra,...
    d_copy,r_copy,rec_image_copy,...
    d_inter,r_inter,rec_image_inter,...
    2^l);
        rec_video{frame} = rec_image;
        D(l-2) = D(l-2) + sum(d,"all");
        R(l-2) = R(l-2) + sum(r,"all");
        cnt_intra(l-2) = cnt_intra(l-2) + length(find(choice==1));
        cnt_copy(l-2) = cnt_copy(l-2) + length(find(choice==2));
        cnt_inter(l-2) = cnt_inter(l-2) + length(find(choice==3));
    end
end

R = R/50/99*144*176*30/1000;
D = D/50/99;
PSNR = 10*log10(255*255./D);
plot(R,PSNR,'ro-')




function [d,r,rec_image,choice] = mode_choice(d_intra,r_intra,rec_image_intra,...
    d_copy,r_copy,rec_image_copy,...
    d_inter,r_inter,rec_image_inter,...
    steplen)

M = 176/16;
N = 144/16;

lambda = 0.2 * steplen^2;

d = zeros(9,11);
r = zeros(9,11);
rec_image = zeros(144,176);

J = zeros(9,11,3);
D = zeros(9,11,3);
R = zeros(9,11,3);
REC = zeros(144,176,3);

J(:,:,1) = d_intra + lambda * r_intra;
J(:,:,2) = d_copy + lambda * r_copy;
J(:,:,3) = d_inter + lambda * r_inter;

D(:,:,1) = d_intra;
D(:,:,2) = d_copy;
D(:,:,3) = d_inter;

R(:,:,1) = r_intra;
R(:,:,2) = r_copy;
R(:,:,3) = r_inter;

REC(:,:,1) = rec_image_intra;
REC(:,:,2) = rec_image_copy;
REC(:,:,3) = rec_image_inter;


[~, choice] = min(J,[],3);
for i = 1:M
    for j = 1:N
        d(j,i) = D(j,i,choice(j,i));
        r(j,i) = R(j,i,choice(j,i));
        rec_image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)) = REC((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16),choice(j,i));
    end
end

end


function [d,r,rec_image] = intra_mode(image,steplen)

M = 176/16;
N = 144/16;

d = zeros(N,M);
r = zeros(N,M);

dct_co = blkproc(image, [8 8], @dct2);
dct_co_qtz = midtread(dct_co,steplen);
rec_image = blkproc(dct_co_qtz, [8 8], @idct2);

for i = 1:M
    for j = 1:N
        d(j,i) = immse(rec_image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)),image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)));
        r(j,i) = 2/256 + entr(dct_co_qtz((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)));
    end
end

end

function [d,r,rec_image] = copy_mode(image,last_image)

M = 176/16;
N = 144/16;

d = zeros(N,M);
r = zeros(N,M);

% last_image denotes the reconstructed image of the last frame.
rec_image = last_image;

for i = 1:M
    for j = 1:N
        d(j,i) = immse(rec_image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)),image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)));
        r(j,i) = 2/256;
    end
end

end

function [d,r,rec_image] = inter_mode(image,last_image,steplen)

% last_image denotes the reconstructed image of the last frame.
M = 176/16;
N = 144/16;

d = zeros(N,M);
r = zeros(N,M);

[mot_vec, ~, ~] = ME_ES(image, last_image, 16, 10);
frame_pre = motionComp(last_image, mot_vec, 16);
resi = image - frame_pre;

dct_co_res = blkproc(resi, [8 8], @dct2);
dct_co_res_qtzd = midtread(dct_co_res, steplen);
rec_resi = blkproc(dct_co_res_qtzd, [8 8], @idct2);
rec_image = frame_pre + rec_resi;

for i = 1:M
    for j = 1:N
        d(j,i) = immse(rec_image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)), image((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)));
        r(j,i) = 2/256 + entr(dct_co_res_qtzd((j - 1) * 16 + (1:16),(i - 1) * 16 + (1:16)));
    end
end


end








%% FS算法：全搜索（Full Search/Exhaustive Search）
function [motionVect,blk_center,costs] = ME_ES(imgP, imgI, mbSize, dm)
% Input
%   	img : 当前帧
%   	imgI : 参考帧
%   	mbSize : MB尺寸
%   	dm : 搜索窗口大小（2dm+1）×（2dm+1）
% Ouput
%   	motionVect : 整像素精度MV
%   	EScomputations: 搜索每个宏块所需的平均点数
    [row, col] = size(imgP);
    blk_center = zeros(2, row*col/(mbSize^2)); 
    %定义每个宏块中心点位置
    motionVect = zeros(2,row*col/(mbSize^2)); %定义每个宏块运动矢量
    costs = ones(2*dm+1,2*dm+1)*20000000;
    computations = 0; %搜索的点数之和
    mb_cnt= 1;
    for i = 1:mbSize:row-mbSize+1     %当前帧起始行搜索范围，步长是块数 
        for j = 1:mbSize:col-mbSize+1 %当前帧起始列搜索范围，步长是块数
            for m= -dm: dm
                for n= -dm: dm
                    ref_blk_row = i+m; %参考帧搜索框起始行
                    ref_blk_col = j+n; %参考帧搜索框起始列
                    %超出搜索范围
                    if (ref_blk_row<1||ref_blk_row+mbSize-1>row||ref_blk_col<1||ref_blk_col+mbSize-1>col)                     
                        continue;                                                            
                    end
                    %计算SAD
                    costs(m+dm+1,n+dm+1) =...
                        costSAD(imgP(i:i+mbSize-1,j:j+mbSize-1),imgI(ref_blk_row:ref_blk_row+mbSize-1,ref_blk_col:ref_blk_col+mbSize-1));
                    computations = computations+1;
                end
            end
            blk_center(1,mb_cnt) = i+ mbSize/2-1; %记录中心点行坐标                     
            blk_center(2,mb_cnt) = j+ mbSize/2-1; %记录中心点列坐标
            [dx,dy,~]=minCost(costs); %找出有最小代价的块的下标
            motionVect(1,mb_cnt) = dx-dm-1; %垂直运动矢量
            motionVect(2,mb_cnt) = dy-dm-1; %水平运动矢量
            mb_cnt = mb_cnt+1;
            costs = ones(2*dm+1,2*dm+1)*20000000; %重新赋值
         end
    end 
end

%% 求预测帧的函数：由给定的运动矢量进行运动补偿计算预测帧
function imgComp = motionComp(imgI, motionVect, mbSize)
% Input
%       imgI : 参考帧
%       motionVect : MV（dx为垂直分量，dy为水平分量）
%   	mbSize : MB尺寸
% Ouput
%   	imgComp : 运动补偿后的图像
    [row,col]=size(imgI);
    mb_cnt=1;
    for i = 1:mbSize:row-mbSize+1                
        for j = 1:mbSize:col-mbSize+1 
             ref_blk_row=i+motionVect(1,mb_cnt); %参考帧搜索块起始行
             ref_blk_col=j+motionVect(2,mb_cnt); %参考帧搜索块起始列
             imgComp(i:i+mbSize-1,j:j+mbSize-1)=imgI(ref_blk_row:ref_blk_row+mbSize-1,ref_blk_col:ref_blk_col+mbSize-1);   
             mb_cnt=mb_cnt+1;
        end
    end
end

%% 求具有最小SAD值的函数：找出具有最小代价的块的下标
function [dx,dy,minc] = minCost(costs)
% Input
%       costs : 包含当前宏块所有运动估计误差代价的SAD矩阵
% Output
%       dx : MV的垂直分量（行位移）
%       dy : MV的水平分量（列位移）
    minc = min(min(costs));
    [dx, dy] = find(costs == minc);
    dx = dx(1);
    dy = dy(1);    
end

%% SAD计算函数：对给定的两个块计算SAD值
function cost = costSAD(currentBlk,refBlk)
% Input
%       currentBlk : 当前块
%       refBlk : 参考块
%       mbSize : MB尺寸
% Output
%       cost : 两个块之间的误差代价（SAD）
    cost=sum(sum((currentBlk-refBlk).^2)); 
end