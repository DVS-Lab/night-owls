clear; close all;

% set up dirs
codedir = pwd; % must run from code, so this is not a good solution
cd ..
maindir = pwd;
%task='sharedreward';
evdir = fullfile(maindir,'derivatives','fsl','EVFiles');


% load sub/run list
sub = [101 103 104];%105
session=[9 12 12]; 
runs=2;
log={}

for s = 1:length(sub)
    for ses=1:session(s)
        sesname=sprintf('ses-%02d', ses);
        for r = 1:runs
            rundir = fullfile(evdir,['sub-' num2str(sub(s))],sesname,'sharedreward',['run-' num2str(r)]);

            if ~exist(rundir, 'dir')
                subchar=num2str(s);
                seschar=num2str(ses);
                log{end+1}=sprintf('sub %s %s run %s does not exist.',num2str(sub(s)),sesname,num2str(r));
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

            all_evs = [ev1; ev2; ev3; ev4; ev5; ev6;ev7];
            all_evs = sortrows(all_evs,1,'ascend');
 
            [~, miss_idx] = ismember(ev7, all_evs, 'rows');
            fprintf('sub %s %s run %s missed trial number: %s\n', num2str(sub(s)), sesname, num2str(r), mat2str(miss_idx(miss_idx > 0)))
%%%%%%take note of which trial has missing outcome event and add into run
%%%%%%script as the trial that gets skipped(continued)
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
            % check length of trials. everyone should have 64
            if length(all_evs) ~= 54
                disp(sprintf('CHECK: sub %s %s run %s missing trials even after missed outcomes counted', num2str(sub(s)), sesname, num2str(r)))
                keyboard
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
                dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f')

                % write out other trials
                fname = sprintf('run-%d_OtherTrials%02d.txt',r,t);
                dlmwrite(fullfile(outdir,fname),othertrials,'delimiter','\t','precision','%.6f')
            end
        end
    end
end
cd(codedir);



