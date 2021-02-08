%function elf_stats_areatest

E = (-89:2:89)';
R1 = 1*randn(length(E), 30);
%R1 = repmat(E/90, [1, 30]) + 1*randn(length(E), 30);
R1_med = median(R1, 2);
R2 = 1*randn(length(E), 52);
%R2 = repmat(-E/90, [1, 52]) + 1*randn(length(E), 52);
R2_med = median(R2, 2);


%% parameters
withrepl = false; %quicker
reps = 10000;
n1 = size(R1, 2);
n2 = size(R2, 2);

test1_ind1 = randi(n1, reps);
test2_ind1 = randi(n2, reps);
test3_ind1 = randi(n1, reps);
test3_ind2 = randi(n2, reps);

if withrepl
    test1_ind2 = randi(n1, reps);
    test2_ind2 = randi(n2, reps);
else
    test1_ind2 = randi(n1-1, reps);
    test1_ind2(test1_ind2>=test1_ind1) = test1_ind2(test1_ind2>=test1_ind1) + 1;
    test2_ind2 = randi(n2-1, reps);
    test2_ind2(test2_ind2>=test2_ind1) = test2_ind2(test2_ind2>=test2_ind1) + 1;
end


D = zeros(reps, 3);

for i = 1:reps
    if mod(i, 1000)==0, disp(i), end
    
    test1_sample1 = R1(:, test1_ind1(i));
    test1_sample2 = R1(:, test1_ind2(i));
    D(i, 1) = sum(abs(test1_sample1 - test1_sample2));
    
    test2_sample1 = R2(:, test2_ind1(i));
    test2_sample2 = R2(:, test2_ind2(i));
    D(i, 2) = sum(abs(test2_sample1 - test2_sample2));
        
    test3_sample1 = R1(:, test3_ind1(i));
    test3_sample2 = R2(:, test3_ind2(i));
    D(i, 3) = sum(abs(test3_sample1 - test3_sample2));
end



DM = sum(abs(R1_med - R2_med))




%% Plot
figure(1); clf;
subplot(1, 2, 1);
plot(R1, E);
hold on;
plot(R1_med, E, 'r', 'linewidth', 3);
xlabel('Radiance environment 1');
ylabel('Elevation (\circ)');

subplot(1, 2, 2);
plot(R2, E);
hold on;
plot(R2_med, E, 'r', 'linewidth', 3);
xlabel('Radiance environment 2');
ylabel('Elevation (\circ)');

figure(2); clf;
subplot(3, 1, 1);
hist(D(:, 1), 50); xlim([0 200]);
subplot(3, 1, 2);
hist(D(:, 2), 50); xlim([0 200]);
subplot(3, 1, 3);
hist(D(:, 3), 50); xlim([0 200]);

h1 = kstest2(D(:, 3), D(:, 1))
h2 = kstest2(D(:, 3), D(:, 2))
%end




% function D = sub_area(R1, R2, E)
% 
% binwidth = median(diff(E));
% D = sum(abs(R1-R2)*binwidth);
% 
% end