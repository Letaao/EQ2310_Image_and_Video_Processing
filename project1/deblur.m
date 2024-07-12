function f_hat = deblur(g, h, n_var)
    [g_M, g_N] = size(g);
    [h_M, h_N] = size(h);
    
    % padd pixels arround the image so the boundaries 
    g_padded = padarray(g, size(h), 'symmetric');
    
    %compute the constant K
    K=n_var/var(g(:));
    [M, N] = size(g_padded);
    %M = g_M+h_M;
    %N = g_N+h_N;
    %g_padded=g;

    % to the frequency domain
    H = fft2(h, M, N);
    G = fft2(g_padded,M,N);
    % equation of the wiener filter
    filter = (1 ./ H) .* ((abs(H).^2) ./ (abs(H).^2 + K));
    % apply wiener filter to the blurred image
    F_hat = G .* filter;
    f_hat = abs(ifft2(F_hat));
    %remove the padding
    f_hat = f_hat(ceil(double(h_M) / 2):ceil(double(h_M) / 2) + g_M - 1, ceil(double(h_N) / 2):ceil(double(h_N) / 2) + g_N - 1);
end
