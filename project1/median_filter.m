function [ img ] = median_filter( image, m )
%----------------------------------------------
%median
%input:
%image:original
%m:size of kernel
 
%output:
%img: the output image
%----------------------------------------------
    n = m;
    [ height, width ] = size(image);
    x1 = double(image);
    x2 = x1;
    for i = 1: height-n+1
        for j = 1:width-n+1
            mb = x1( i:(i+n-1),  j:(j+n-1) );%acquire n*n matrix
            mb = mb(:);
            mm = median(mb);%taking the median value
            x2( i+(n-1)/2,  j+(n-1)/2 ) = mm;
 
        end
    end
 
    img = uint8(x2);
end