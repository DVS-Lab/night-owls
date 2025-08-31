clear; close all;

% Resolve paths relative to THIS script
codedir = fileparts(mfilename('fullpath'));   % folder containing this .m file
maindir = fileparts(codedir);                 % go up one level (adjust if needed)
evdir   = fullfile(maindir,'derivatives','fsl','EVFiles');


% load sub/run list
sub = [101 103 104 105];
session = [11 12 12 12];
runs = 2;


for s = 1:length(sub)
    for ses=1:session(s)
        sesname=sprintf('ses-%02d', ses);
        for r = 1:runs
            rundir = fullfile(evdir,['sub-' num2str(sub(s))],sesname,'sharedreward',['run-' num2str(r)]);

            if ~exist(rundir, 'dir')
                subchar=num2str(s);
                seschar=num2str(ses);
                fprintf('sub %s %s run %s does not exist.\n',num2str(sub(s)),sesname,num2str(r));
                continue;
            end

            % load evs and concatenate
            misspath = fullfile(rundir,"_miss_outcome.txt");

            if exist(misspath, 'file')==2
                ev1=load(fullfile(rundir,'_outcome_computer_neutral.txt'));
                ev2=load(fullfile(rundir,'_outcome_computer_punish.txt'));
                ev3=load(fullfile(rundir,'_outcome_computer_reward.txt'));
                ev4=load(fullfile(rundir,'_outcome_stranger_neutral.txt'));
                ev5=load(fullfile(rundir,'_outcome_stranger_punish.txt'));
                ev6=load(fullfile(rundir,'_outcome_stranger_reward.txt'));
                ev7=load(fullfile(rundir,'_miss_outcome.txt'));

                all_evs = [ev1; ev2; ev3; ev4; ev5; ev6; ev7];
                all_evs = sortrows(all_evs,1,'ascend');

            else

                ev1=load(fullfile(rundir,'_outcome_computer_neutral.txt'));
                ev2=load(fullfile(rundir,'_outcome_computer_punish.txt'));
                ev3=load(fullfile(rundir,'_outcome_computer_reward.txt'));
                ev4=load(fullfile(rundir,'_outcome_stranger_neutral.txt'));
                ev5=load(fullfile(rundir,'_outcome_stranger_punish.txt'));
                ev6=load(fullfile(rundir,'_outcome_stranger_reward.txt'));


                all_evs = [ev1; ev2; ev3; ev4; ev5; ev6];
                all_evs = sortrows(all_evs,1,'ascend');
            end
            % check length of trials. everyone should have 54
            if length(all_evs) ~= 54
                fprintf('CHECK: %s %s run %s missing trials even after missed outcomes counted\n', num2str(sub(s)), sesname, num2str(r))
            end

            % extract trials and write evs
            outdir = fullfile(evdir,['sub-' num2str(sub(s))],'singletrial',sesname,'sharedreward',['run-' num2str(r)]);

            if ~exist(outdir,'dir')
                mkdir(outdir);
            end

            for t = 1:length(all_evs)
                singletrial = all_evs(t,:);
                othertrials = all_evs;
                othertrials(t,:) = []; % delete trial

                % write out single trial
                fname = sprintf('run-%d_SingleTrial%02d.txt',r,t);
                writematrix(round(singletrial, 2), fullfile(outdir, fname), 'Delimiter','tab');

                % write out other trials
                fname = sprintf('run-%d_OtherTrials%02d.txt',r,t);
                writematrix(round(othertrials, 2), fullfile(outdir, fname), 'Delimiter','tab');

            end
        end
    end
end