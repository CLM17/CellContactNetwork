addpath('BrainConnectivity');
A = adjacency(G);
alpha = 0.5;
[R] = randomizer_bin_und(A,alpha);

figure()
G_random = graph(R);
p1 = plot(G_random, 'XData', xNodes, 'YData', yNodes, 'markersize',2);