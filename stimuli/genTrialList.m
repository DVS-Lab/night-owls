%function genTrialList(s)
%{

trials per run: 54
length of a trial: 5.1

trial sequence              seconds
decision phase with face	2.8	max (remaining goes into ISI)
ISI                         1.7	mean
outcome                     0.6
total                       5.1


evt1	computer_loss       11
evt2	computer_neutral	5
evt3	computer_win        11
evt4	stranger_loss       11
evt5	stranger_neutral	5
evt6	stranger_win        11


%}




maindir = pwd;
outfiles = fullfile(maindir,'event-related','params');
mkdir(outfiles);
subs = 10804

for s = subs

    subout = fullfile(outfiles,sprintf('sub-%04dsp',s));
    mkdir(subout);

    runs = 6;
    ntrials = 54;

    trial_types = [repmat([1 3 4 6],1,11) repmat([2 5],1,5)];
    ISI_distribution = repmat([0.85 1.7 2.55],1,18);
    ISI_distribution = ISI_distribution(randperm(length(ISI_distribution)));
    ITI_distribution = repmat([1 2 3 4],1,18);
    ITI_distribution = ITI_distribution(randperm(length(ITI_distribution)));

    for r = 1:runs
    
        rand_trials = randperm(ntrials);
        fname = fullfile(subout,sprintf('sub-%04dsp_run-%d_design.csv',s,r));
        fid = fopen(fname,'w');
        fprintf(fid,'Trialn,TrialType,Partner,Feedback,ITI,ISI\n');
        for t = 1:ntrials
            tt = rand_trials(t);
            switch trial_types(tt)
                case 1 %Computer Punishment
                    partner = 1;
                    feedback_mat = 1;
                case 2 %Computer neutral
                    partner = 1;
                    feedback_mat = 2;
                case 3 %Computer Reward
                    partner = 1;
                    feedback_mat = 3;
                case 4 %Stranger Punishment
                    partner = 2;
                    feedback_mat = 1;
                case 5 %Stranger Neutral
                    partner = 2;
                    feedback_mat = 2;
                case 6 %Stranger Reward
                    partner = 2;
                    feedback_mat = 3;
            end
            fprintf(fid,'%d,%d,%d,%d,%d,%d\n',t,trial_types(tt),partner,feedback_mat,ITI_distribution(tt),ISI_distribution(tt));
        end
        fclose(fid);
    end
end