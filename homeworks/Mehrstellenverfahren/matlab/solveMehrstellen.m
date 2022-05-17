function u = solveMehrstellen(f, M)

    A = compMehrstellenA(M);
    f = compMehrstellenf(f, M);
    u = A \ f;

end
