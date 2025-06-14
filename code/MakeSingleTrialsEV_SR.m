clear; close all;

% set up dirs
codedir = pwd; % must run from code, so this is not a good solution
cd ..
maindir = pwd;
task=sharedreward;
evdir = fullfile(maindir,'derivatives','fsl','EVFiles');

sub-101/ses-01/mid/run-1

% load sub/run list
sub = [101 103 104];
session=[9 12 10]; 
runs=2;

for s = 1:length(sub)
for ses=1:session(s)
    for r = 1:runs
        % load evs and concatenate

         ev1=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_computer_neutral.txt'))
         ev2=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_computer_punish.txt'))
         ev3=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_computer_reward.txt'))
         ev4=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_stranger_neutral.txt'))
         ev5=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_stranger_punish.txt'))
         ev6=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'sharedreward',['run-' num2str(r)],'_outcome_stranger_reward.txt'))

        all_evs = [ev1; ev2; ev3; ev4; ev5; ev6];
        all_evs = sortrows(all_evs,1,'ascend');
        
        % check length of trials. everyone should have 64
        if length(all_evs) ~= 54
            disp(sprintf('sub %s ses %s run %s smissing trials...', num2str(sub(s)), num2str(ses), num2str(r)));
            keyboard
        end
        
        % extract trials and write evs
        outdir = fullfile(evdir,['sub-' num2str(sub(s))],'singletrial',['ses-0' num2str(ses)],'mid',['run-' num2str(r)]);
        
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end

        for t = 1:length(all_evs)
            singletrial = all_evs(t,:);
            othertrials = all_evs;
            othertrials(t,:) = []; % delete trial
            
            % write out single trial
            fname = sprintf('run-0%d_SingleTrial%02d.txt',r,t);
            dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f') 
            
            % write out other trials
            fname = sprintf('run-0%d_OtherTrials%02d.txt',r,t);
            dlmwrite(fullfile(outdir,fname),othertrials,'delimiter','\t','precision','%.6f') 
        end
   end 
end
cd(codedir);



