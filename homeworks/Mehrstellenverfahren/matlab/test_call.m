M=3;
fprintf('\n#compMehrstellenA');
A=compMehrstellenA(M)
fprintf('\n#compMehrstellenf');
func = @(x) sin(pi*x(1))*sin(pi*x(2));
f=compMehrstellenf(func,M)
fprintf('\n#solveMehrstellen');
u=solveMehrstellen(func,M)
fprintf('\n#compgriderr');
err=compgriderr(M)