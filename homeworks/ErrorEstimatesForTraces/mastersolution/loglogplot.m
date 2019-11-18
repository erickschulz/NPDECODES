% Creating convergence plot for homework problem TraceErrorEstimates

results = csvread('results.csv');

figure('name','convergence plot');
loglog(results(:,2),results(:,1),'r+-')

p = polyfit(log(results(:,2)),log(results(:,1)),1);
xlabel('{\bf N = dimension of FE space}');
ylabel('{\bf Output error |B(u) - B(u_h)|}');
ratebox = sprintf('%s%f', 'rate: ',p(1));
legend(ratebox,'Location','southwest');
set(gca, 'YDir','reverse');

print -depsc2 'convergence_plot.eps';
