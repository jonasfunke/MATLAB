%% 
clear all, close all, clc
run('my_prefs')

%%

k1 = 10; % 1/nM/s
k2 = 1000; % 1/nM/s
fs = @(t, y) linking_on_fs(t, y, [k1 k2]);


L_0 = [0 25 50 75 100 150 200 250 500];
FS_0 = 50;

t0 = 0;
t_end = log(100)/k1/L_0(2);




options = odeset('RelTol',1e-4,'AbsTol',[1e-2 1e-2 1e-2 1e-2 1e-2]);
yield = cell(length(L_0), 2);

for i = 1:length(L_0)
    state_0 = [L_0(i) FS_0 0 0 0]

    [T,Y] = ode45(fs, [t0 t_end], state_0, options);
    yield{i, 1} = T;
    yield{i, 2} = Y(:,5) ./ (Y(:,2) + Y(:,3) + Y(:,4) + Y(:,5));


end


%
close all

cc = varycolor(length(L_0));
myleg = cell(length(L_0),1);
hold all
for i=1:length(L_0)
    plot(yield{i,1}, yield{i,2}, 'Color', cc(i,:))
    myleg{i} = num2str(L_0(i));
end
legend(myleg)

%%


close all
subplot(2,1,1)
plot(T,Y(:,1),'-',T,Y(:,2),'-.',T,Y(:,3),'.', T,Y(:,4),'.',T,Y(:,5),'.')
legend({'Linker', 'FS', 'FS+1', 'FS+2', 'FS linked'})

subplot(2,1,2)
plot(T, 100*yield)