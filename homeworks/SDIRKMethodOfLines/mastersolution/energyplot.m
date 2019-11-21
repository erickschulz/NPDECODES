% Creating convergence plot for homework problem TraceErrorEstimates

energies = csvread('energies.csv');
t = linspace(0.0,1.0, length(energies));

figure('name','energy plot');
plot(t(:),energies(:),'r+-')
xlabel('{\bf Time}');
ylabel('{\bf Output energy}');
title(sprintf('m = %d, N=%d',length(energies)-1,10201));

print -depsc2 'energy_plot.eps';
