function f = compMehrstellenf(f, M)

    S = gallery('tridiag', M, 1, 8, 1);
    T = gallery('tridiag', M, 1, 0, 1);
    I = speye(M);
    A = kron(I, S) + kron(T, I);

    % Do not assume that f is vector-safe
    L = zeros(M*M,1);
    xpts = linspace(0, 1, M+2);
    xpts = xpts(2:end-1);
    i = 1;
    for y = xpts
        for x = xpts
            L(i) = f([x, y]);
            i = i+1;
        end
    end

    f = A * L;
    f = f / (M+1)^2 / 12;

end
