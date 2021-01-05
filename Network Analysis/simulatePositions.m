R = 6300;
N = 1000;
s = 100;
T = 10000;

initR = rand(1,N) * R;
initTheta = rand(1,N) * 2*pi;
r = [initR .* cos(initTheta); initR .* sin(initTheta)];

dr = rand(2,N) - 0.5;
for i = 1:N
    dr(:,i) = dr(:,i) / sqrt(dr(:,i)' * dr(:,i));
end

disp(dr);

sigma = 0.05;
edgeWell = linspace(0,2*pi,100);

figure(1)
xlim([-R/2, R/2])
ylim([-R/2, R/2])

for t = 1:T
    
    %theta = randn(N,1) * sigma;
    
    for i = 1:N
        
        %while true
        theta = randn * sigma;
        M = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        dr(:,i) = M * dr(:,i);

        if sqrt( (r(:,i) + s*dr(:,i))' * (r(:,i) + s*dr(:,i)) ) > R
            dr(:,i) = - dr(:,i);
        end
        %end

        r(:,i) = r(:,i) + s*dr(:,i);
    end
    
    plot(r(1,:), r(2,:), 'o')
    hold on
    plot(R*cos(edgeWell), R*sin(edgeWell), '-r')
    hold off
    drawnow
    %pause(0.005)
    
end