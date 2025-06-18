clear; close all;

% set up dirs
codedir = pwd; % must run from code, so this is not a good solution
cd ..
maindir = pwd;
<<<<<<< HEAD
%task='mid'
=======
task=mid
>>>>>>> 18f3f9a665eb565ba1486ae11d91c53bf4f60a18
evdir = fullfile(maindir,'derivatives','fsl','EVFiles');

% load sub/run list
sub = [101 103 104];
session=[9 12 10];
runs=2;
<<<<<<< HEAD
log={}

for s = 1:length(sub)
    for ses=1:session(s)
        for r = 1:runs
         rundir = fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)]);
        
         if ~exist(rundir, 'dir')
            subchar=num2str(s)
            seschar=num2str(ses)
            log{end+1}=sprintf('sub %s ses %s run %s does not exist.',num2str(s),num2str(ses),num2str(r));
            continue;
         end

            % load evs and concatenate
            ev1=load(fullfile(rundir,'_anticipation_neutral.txt'))
            ev2=load(fullfile(rundir,'_anticipation_reward.txt'))
=======

for s = 1:length(subrun)
    for ses=1:session(s)
        for r = 1:runs
            % load evs and concatenate
            ev1=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)],'_anticipation_neutral.txt'))
            ev2=load(fullfile(evdir,['sub-' num2str(sub(s))],['ses-0' num2str(ses)],'mid',['run-' num2str(r)],'_anticipation_reward.txt'))
>>>>>>> 18f3f9a665eb565ba1486ae11d91c53bf4f60a18
            all_evs = [ev1; ev2];
            all_evs = sortrows(all_evs,1,'ascend');

            % check length of trials. everyone should have 64
            if length(all_evs) ~= 56
<<<<<<< HEAD
                disp(sprintf('sub %s ses %s run %s missing trials...', num2str(sub(s)), num2str(ses), num2str(r)));
=======
                disp('missing trials...')
>>>>>>> 18f3f9a665eb565ba1486ae11d91c53bf4f60a18
                keyboard
            end

            % extract trials and write evs
<<<<<<< HEAD
            outdir = fullfile(evdir,['sub-' num2str(sub(s))],'singletrial',['ses-0' num2str(ses)],'mid',['run-' num2str(r)]);
=======
            outdir = fullfile(evdir,['sub-' num2str(subnum)],'singletrial',['ses-0' num2str(ses(ss))],'mid',['run-' num2str(r)]);
>>>>>>> 18f3f9a665eb565ba1486ae11d91c53bf4f60a18
            if ~exist(outdir,'dir')
                mkdir(outdir);
            end
            for t = 1:length(all_evs)
                singletrial = all_evs(t,:);
                othertrials = all_evs;
                othertrials(t,:) = []; % delete trial

                % write out single trial
<<<<<<< HEAD
                fname = sprintf('ses-0%drun-%d_SingleTrial%02d.txt',ses,r,t);
                dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f')

                % write out other trials
                fname = sprintf('ses-0%drun-%d_OtherTrials%02d.txt',ses,r,t);
=======
                fname = sprintf('ses-0%drun-%d_SingleTrial%02d.txt',d,r,t);
                dlmwrite(fullfile(outdir,fname),singletrial,'delimiter','\t','precision','%.6f')

                % write out other trials
                fname = sprintf('ses-0%drun-%d_OtherTrials%02d.txt',d,r,t);
>>>>>>> 18f3f9a665eb565ba1486ae11d91c53bf4f60a18
                dlmwrite(fullfile(outdir,fname),othertrials,'delimiter','\t','precision','%.6f')
            end
        end
    end
end
cd(codedir);



