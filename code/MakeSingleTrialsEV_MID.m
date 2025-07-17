clear; close all;

% set up dirs
codedir = pwd; % must run from code, so this is not a good solution
cd ..
maindir = pwd;
%task='mid'
evdir = fullfile(maindir,'derivatives','fsl','EVFiles');

% load sub/run list
sub = [101 103 104];
session=[9 12 12];
runs=2;
log={}

for s = 1:length(sub)
    for ses=1:session(s)
        for r = 1:runs
         rundir = fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)]);
        
         if ~exist(rundir, 'dir')
            subchar=num2str(s);
            seschar=num2str(ses);
            log{end+1}=sprintf('sub %s ses %s run %s does not exist.',num2str(s),num2str(ses),num2str(r));
            continue;
         end

            % load evs and concatenate
            ev1=load(fullfile(rundir,'_anticipation_neutral.txt'));
            ev2=load(fullfile(rundir,'_anticipation_reward.txt'));
            all_evs = [ev1; ev2];
            all_evs = sortrows(all_evs,1,'ascend');

    

            % check length of trials. everyone should have 64
            if length(all_evs) ~= 56
                disp(sprintf('sub %s ses %s run %s missing trials...', num2str(sub(s)), num2str(ses), num2str(r)));
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
                fname = sprintf('ses-0%drun-%d_SingleTrial%02d.txt',ses,r,t);
                dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f')

                % write out other trials
                fname = sprintf('ses-0%drun-%d_OtherTrials%02d.txt',ses,r,t);
                dlmwrite(fullfile(outdir,fname),othertrials,'delimiter','\t','precision','%.6f')
            end
        end
    end
end
cd(codedir);



