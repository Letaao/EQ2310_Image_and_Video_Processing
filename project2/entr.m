function y = entr(img) % calculate the bit-rate of a given matrix
    [M, N] = size(img);
    [pixel_counts, ~] = groupcounts(img(:)); % find the counting of all levels
    pixel_probas = nonzeros(pixel_counts(:)./(M*N));
    log_pixel_probas = log2(pixel_probas);
    y = -sum(pixel_probas.*log_pixel_probas);
end

