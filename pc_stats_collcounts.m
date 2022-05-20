function A = pc_stats_collcounts(a,b)
A = zeros(b,2*a+2);

for i = 1:b
    rng default;
    reset(RandStream.getGlobalStream,sum(100*clock));
    r = 10+100*rand();
    x = pc1(1000,a,10,10,r);
    
    %FOR PERCENTAGES OF EACH TRIAL%
    A(i,1:a) = x(1:a)./sum(x(1:a),2).*100;
    A(i,a+1) = r;
    A(i,end-a:end-1) = x(end-a+1:end);
    A(i,end) = sum(x(1:a),2);
end

figure(1);
for i = 1:5
    pointsize = 10;
    subplot(3,2,i);
    scatter(A(:,i+6),A(:,i),pointsize,A(:,end),'filled');
    cb = colorbar();
    title(['Percent of Collisions vs Mass of particle ', num2str(i)]);
    ylabel('Collision Percentage');
    xlabel('Mass of Particle');
end

figure(2);
for i = 1:5
    subplot(3,2,i)
    pointsize = 10;
    scatter(A(:,6),A(:,i),pointsize,A(:,end),'filled');
    cb = colorbar();
    title(['Percent of Collisions vs Right Energy '])
    ylabel('Collision Percentage')
    xlabel('Right Energy')
end