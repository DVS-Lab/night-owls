clear; close all;

% set up dirs
codedir = pwd; % must run from code, so this is not a good solution
cd ..
maindir = pwd;
task=mid
evdir = fullfile(maindir,'derivatives','fsl','EVFiles');

% load sub/run list
sub = [101 103 104];
session=[9 12 10];
runs=2;

for s = 1:length(subrun)
    for ses=1:session(s)
        for r = 1:runs
            % load evs and concatenate
            ev1=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)],'_anticipation_neutral.txt'))
            ev2=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)],'_anticipation_reward.txt'))
            all_evs = [ev1; ev2];
            all_evs = sortrows(all_evs,1,'ascend');

            % check length of trials. everyone should have 64
            if length(all_evs) ~= 56
                disp('missing trials...')
                keyboard
            end

            % extract trials and write evs
            outdir = fullfile(evdir,['sub-' num2str(subnum)],'singletrial',['ses-0' num2str(ses(ss))],'mid',['run-' num2str(r)]);
            if ~exist(outdir,'dir')
                mkdir(outdir);
            end
            for t = 1:length(all_evs)
                singletrial = all_evs(t,:);
                othertrials = all_evs;
                othertrials(t,:) = []; % delete trial

                % write out single trial
                fname = sprintf('ses-0%drun-%d_SingleTrial%02d.txt',d,r,t);
                dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f')

                % write out other trials
                fname = sprintf('ses-0%drun-%d_OtherTrials%02d.txt',d,r,t);
                dlmwrite(fullfile(outdir,fname),othertrials,'delimiter','\t','precision','%.6f')
            end
        end
    end
end
cd(codedir);



