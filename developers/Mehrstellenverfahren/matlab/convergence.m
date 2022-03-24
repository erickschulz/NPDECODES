M = [5, 10, 20, 40, 80, 160];
errs = [];

for m = M
    errs = [errs, compgriderr(m)];
end

loglog(M, errs, 'LineWidth', 2);
xlabel('M'); ylabel('Error'); grid on;
p = polyfit(log(M), log(errs), 1);
disp(['Order: ' num2str(-p(1))]);
