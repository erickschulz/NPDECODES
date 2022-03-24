function err = compgriderr(M)

    f = @(x) sin(pi*x(1)) * sin(pi*x(2));
    u = solveMehrstellen(f, M);

    xpts = linspace(0, 1, M+2);
    xpts = xpts(2:end-1);
    i = 1;
    err = 0;
    for y = xpts
        for x = xpts
            m = abs(u(i) - f([x, y])/(2*pi^2));
            err = max(m, err);
            i = i+1;
        end
    end

end
