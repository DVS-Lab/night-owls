% Generate two runs, 105 trials each

% load run1 from Deepu's lab
load('run1.mat')

% keep uniform distribution of conditions
cond = repelem([1,2,3], 35);

% add new ISIs with same properties as original
isi11 = zeros(105,1);
isi12 = zeros(105,1);
isi21 = zeros(105,1);
isi22 = zeros(105,1);

sess = 1:12;
timing = cell(length(sess),2);

rng(12345);
for ii = 1:length(sess)
    cond_rand1 = cond(randperm(numel(cond)));
    cond_rand2 = cond(randperm(numel(cond)));
    for i = 1:length(isi11)
        rand_idx = randperm(75);
        isi11(i,1) = run.isi1(rand_idx(1));
        isi12(i,1) = run.isi2(rand_idx(1));
        rand_idx = randperm(75);
        isi21(i,1) = run.isi1(rand_idx(1));
        isi22(i,1) = run.isi2(rand_idx(1));
    end

    % Check timing
    timing{ii,1} = (sum(isi11) + sum(isi12))/60;
    timing{ii,2} = (sum(isi21) + sum(isi22))/60;

    % Put everything together and save
    run.cond = cond_rand1;
    run.isi1 = isi11;
    run.isi2 = isi12;
    fileName = sprintf('ses-%d_run-1_noloss.mat', ii);
    save(fileName,'run')
    run.cond = cond_rand2;
    run.isi1 = isi21;
    run.isi2 = isi22;
    fileName = sprintf('ses-%d_run-2_noloss.mat', ii);
    save(fileName,'run')
end

maxTime = max(cellfun(@(x) max(x(:)), timing(:)));
difTime = cellfun(@(x) (x - maxTime)*60, timing, 'UniformOutput', false);
writecell(difTime,'padtimes.csv');






