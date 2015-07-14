function C = conv2nan(A, B, shape)
    % Computes MATLAB's conv2 function whilst skipping NaN's
    Anan = isnan(A);
    Bnan = isnan(B);
    
    A(Anan) = 0;
    B(Bnan) = 0;
    
    C = conv2(A,B,shape);
    
    N = conv2(real(~Anan),real(~Bnan), shape);     % normalization term
    c = conv2(ones(size(A)),ones(size(B)), shape); % correction of normalization

    C = C.*c./N;
    C(Anan) = NaN;
    
end