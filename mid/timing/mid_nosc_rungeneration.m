% Generate two runs, 120 trials each

% load run1 from Deepu's lab
load('run1.mat')

% keep uniform distribution of conditions
cond = repelem([1,2,3], 40);
cond_rand1 = cond(randperm(numel(cond)));
cond_rand2 = cond(randperm(numel(cond)));

% add new ISIs with same properties as original
isi11 = zeros(120,1);
isi12 = zeros(120,1);
isi21 = zeros(120,1);
isi22 = zeros(120,1);

rng(1225);
for i = 1:length(isi11)
    rand_idx = randperm(75);
    isi11(i,1) = run.isi1(rand_idx(1));
    isi12(i,1) = run.isi2(rand_idx(1));
    rand_idx = randperm(75);
    isi21(i,1) = run.isi1(rand_idx(1));
    isi22(i,1) = run.isi2(rand_idx(1));
end
(sum(isi11) + sum(isi12))/60
(sum(isi21) + sum(isi22))/60


% put everything together and save
run.cond = cond_rand1;
run.isi1 = isi11;
run.isi2 = isi12;
save('run1_noloss.mat','run')
run.cond = cond_rand2;
run.isi1 = isi21;
run.isi2 = isi22;
save('run2_noloss.mat','run')


