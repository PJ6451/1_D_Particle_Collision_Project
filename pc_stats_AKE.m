function A = pc_stats_AKE(a,b)
%a is the number of particles, b is the number of iterations

A = zeros(b,2*a);

for i = 1:b
    rng default;
    reset(RandStream.getGlobalStream,sum(100*clock));
    r = 1000;
    A(i,1:end) = pc_averageKE(1000,a,10,10,r);
end

figure(1);
for i = 1:a
    pointsize = 10;
    subplot(3,2,i);
    scatter(A(:,i+a),A(:,i),pointsize,'filled');
    %cb = colorbar();
    title(['Average Kinetic Energy vs Mass of particle ', num2str(i)]);
    ylabel('Average Kinetic Energy');
    xlabel('Mass of Particle');
end

% figure(2);
% for i = 1:5
%     subplot(3,2,i)
%     pointsize = 10;
%     scatter(A(:,6),A(:,i),pointsize,A(:,end),'filled');
%     cb = colorbar();
%     title(['Percent of Collisions vs Right Energy '])
%     ylabel('Collision Percentage')
%     xlabel('Right Energy')
% end
