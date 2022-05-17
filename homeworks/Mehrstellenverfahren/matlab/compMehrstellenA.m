function A = compMehrstellenA(M)

    S = gallery('tridiag', M, -1, -4, -1);
    A = kron(-S, S);
    A(1:M^2+1:end) = 4 - diag(A);
    A = A/6;

end
