%% 
clear all, close all, clc


figure(1)
hold all
r_mean = zeros(30,1);
var = zeros(30,1);

for n = 1:30;
    L = n .* 0.65; % contour length [nm]
    L_p = 1.5 ; %persistence length [nm]
    t = L ./ L_p;
    r = [0:0.01:0.99.*L] ./ L;
    p = ssDNA_end_to_end_dist(t, r) ;
    p_norm = p ./sum(p);
    plot(r, p_norm)
    r_mean(n) = sum(p_norm.*r);
    r2 = sum(p_norm.*r.^2);
    var(n) = r2-r_mean(n).^2;
    Lbla(n) = L;
    
end

%%
figure(2)
errorbar([1:30], r_mean.*Lbla', sqrt(var).*Lbla',  'k.'), hold on
plot(([10:50]-10), 6.*([10:50])./47 , 'r.')
plot(([10:50]), 16.*([10:50])./47 , 'g.')
plot(([10:50]-10), 27.*([10:50])./47 , 'b.')
xlabel('Number of thymine-bases'), ylabel('Mean end-to-end distance +- std [nm]')
set(gca, 'Xlim', [0 30], 'YLim', [0 10])


%%

tmp = load('/Users/jonasfunke/Documents/Typhoon_images/2014-11-13_TemplatedReaction_BranchMigration/2014-11-13_TemplatedReaction_BM_50um_03-[Cy5]_bands01/2014-11-13_TemplatedReaction_BM_50um_03-[Cy5]_bands01_data.mat')

 x = r_mean.*Lbla';
d = [ 0; x([1:10 20 30]) ];
figure(3)
plot(d, tmp.yield(1:end-1), '.-')
xlabel('Mean distance [nm]')
ylabel('Yield of reaction')
set(gca, 'YLim', [0 0.5])
%%


r_mean = sum(p_norm.*r);
r2 = sum(p_norm.*r.^2);
var = r2-r_mean.^2;

plot(r, p_norm), hold on
vline(r_mean, 'r')
vline(r_mean-sqrt(var), 'r--')
vline(r_mean+sqrt(var), 'r--')