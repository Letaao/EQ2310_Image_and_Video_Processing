function y = entr(img) % calculate the bit-rate of a given matrix
     [M, N] = size(img);
    [pixel_counts, ~] = groupcounts(img(:)); % find the counting of all levels
    pixel_probas = nonzeros(pixel_counts(:)./(M*N));
    log_pixel_probas = log2(pixel_probas);
    y = -sum(pixel_probas.*log_pixel_probas);

end
% function y=entr(img)
% [height,width] = size(img);
% value = reshape(img,[1,height*width]);
% bins = min(value):1:max(value);
% pr = hist(value(:),bins(:));
% prb = pr/sum(pr);
% y = -sum(prb.*log2(prb+eps));
% end