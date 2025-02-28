function out = convertSharedReward2BIDSevents(subj,counterbalance)
% This function converts the raw behavioral output from psychopy into
% the BIDS *_events.tsv file format. It also collects summary information
% about the subject's data in the "out" variable.

% Example convertSharedReward2BIDSevents("10007",1)
% 
% The second batch of counterbalance order options apply to seq. pilot 2.0
% (e.g., cb 1 = cb 21

%{

run-1	run-2	run-3	run-4	run-5	run-6	counterbalance
mb1-me1	mb3-me1	mb6-me1	mb1-me4	mb3-me4	mb6-me4	1
mb3-me1	mb6-me1	mb1-me4	mb3-me4	mb6-me4	mb1-me1	3
mb6-me1	mb1-me4	mb3-me4	mb6-me4	mb1-me1	mb3-me1	5
mb1-me4	mb3-me4	mb6-me4	mb1-me1	mb3-me1	mb6-me1	2
mb3-me4	mb6-me4	mb1-me1	mb3-me1	mb6-me1	mb1-me4	4
mb6-me4	mb1-me1	mb3-me1	mb6-me1	mb1-me4	mb3-me4	6

cb1 run-1_mb-1_me-1

%}

try
    
    switch counterbalance
        case 1, acqs = {'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4'};
            acqsFull = {'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4'};
        case 3, acqs = {'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1'};
            acqsFull = {'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1'};
        case 5, acqs = {'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1'};
            acqsFull = {'mb6me1',	'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1'};
        case 2, acqs = {'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1'};
            acqsFull = {'mb1me4',	'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1'};
        case 4, acqs = {'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4'};
            acqsFull = {'mb3me4',	'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4'};
        case 6, acqs = {'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4'};
            acqsFull = {'mb6me4',	'mb1me1',	'mb3me1',	'mb6me1',	'mb1me4',	'mb3me4'};
        % Seq Pilot 2.0 
        case 21, acqs = {'mb3me4', 'mb3me4', 'mb2me4', 'mb3me1', 'mb3me3', 'mb3me3'};
                 acqsFull = {'mb3me4',	'mb3me4fa50',	'mb2me4',	'mb3me1fa50',	'mb3me3', 'mb3me3ip0'};
        case 22, acqs = {'mb3me4', 'mb2me4', 'mb3me1', 'mb3me3', 'mb3me3', 'mb3me4'};
                 acqsFull = {'mb3me4fa50',	'mb2me4',	'mb3me1fa50',	'mb3me3', 'mb3me3ip0', 'mb3me4'};
        case 23, acqs = {'mb2me4', 'mb3me1', 'mb3me3', 'mb3me3', 'mb3me4', 'mb3me4'};
                 acqsFull = {'mb2me4',	'mb3me1fa50',	'mb3me3', 'mb3me3ip0', 'mb3me4', 'mb3me4fa50'};
        case 24, acqs = {'mb3me1', 'mb3me3', 'mb3me3', 'mb3me4', 'mb3me4', 'mb2me4'};
                 acqsFull = {'mb3me1fa50',	'mb3me3', 'mb3me3ip0', 'mb3me4', 'mb3me4fa50', 'mb2me4'};
        case 25, acqs = {'mb3me3', 'mb3me3', 'mb3me4', 'mb3me4', 'mb2me4', 'mb3me1'};
                 acqsFull = {'mb3me3', 'mb3me3ip0', 'mb3me4', 'mb3me4fa50', 'mb2me4', 'mb3me1fa50'};
    	case 26, acqs = {'mb3me3', 'mb3me4', 'mb3me4', 'mb2me4', 'mb3me1', 'mb3me3'};
                 acqsFull = {'mb3me3ip0', 'mb3me4', 'mb3me4fa50', 'mb2me4', 'mb3me1fa50', 'mb3me3'};
    end
    
    % set up paths
    scriptname = matlab.desktop.editor.getActiveFilename;
    fprintf("running subject: %s counterbalance: %d \n", subj, counterbalance)
    [codedir,~,~] = fileparts(scriptname);
    cd(codedir);
    addpath(codedir);
    cd ..
    dsdir = pwd;
    
    % make default output
    out.ntrials(1) = 0;
    out.ntrials(2) = 0;
    out.nmisses(1) = 0;
    out.nmisses(2) = 0;
    out.nfiles = 0;
    
    % get relative path for source data. repos should be in same dir
    logdir = fullfile(dsdir,'stimuli','logs');
    
    for r = 1:6
        % sub-10008_task-sharedreward_run-1_mb-1_me-1_raw.csv --> sub-10008_task-sharedreward_run-1_acq-mb1me1_raw.csv
        %fname = fullfile(logdir,num2str(subj),sprintf('sub-%04d_task-sharedreward_run-%d_acq-%s_raw.csv',subj,r,acqs{r}));
        fname = fullfile(logdir,subj,sprintf('sub-%04s_task-sharedreward_run-%d_acq-%s_raw.csv',subj,r,acqs{r}));
        
        if r == 1 % only needed for first pass through
            [sublogdir,~,~] = fileparts(fname);
            sublogdir=convertStringsToChars(sublogdir);
            nfiles = dir([sublogdir '/*.csv']);
            out.nfiles = length(nfiles);
        end
        
        if exist(fname,'file')
            T = readtable(fname,'TreatAsEmpty','--');
        else
            fprintf(' ')
            fprintf(subj)
            fprintf(' ')
            fprintf(fname)
            fprintf(' ')
            fprintf('sub-%s_task-sharedreward_run-%d: No data found. Exiting...\n', subj, r)
            %exit
        end
        
        % strip out irrelevant information and missed trials
        T = T(:,{'rt','decision_onset','outcome_onset','InitFixOnset','outcome_offset','Feedback','Partner','resp'});
        goodtrials =  ~isnan(T.resp);
        T = T(goodtrials,:);
        
        if height(T) < 54
            fprintf('incomplete data for sub-%d_run-%d\n', 'subj', r)
        end
        
        start_time = T.InitFixOnset(1);
        onset_decision = T.decision_onset - start_time; % switch to outcome_onset? add regressor for decision? minimal spacing...
        onset_outcome = T.outcome_onset - start_time;
        duration = T.outcome_offset - T.outcome_onset; % outcome
        RT = T.rt;
        Partner = T.Partner;
        feedback = T.Feedback;
        response = T.resp; % Ori: Right index is 2, left index is 7
        
        out.ntrials(r) = height(T);
        out.nmisses(r) = sum(T.resp < 1);

        % output file
        fname = sprintf('sub-%04s_task-sharedreward_acq-%s_events.tsv',subj,acqsFull{r}); % need to make fMRI run number consistent with this?
        output = fullfile(dsdir,'bids',['sub-' num2str(subj)],'func');
        
        if ~exist(output,'dir')
            mkdir(output)
        end
        myfile = fullfile(output,fname);
        fid = fopen(myfile,'w');
        
        
        fprintf(fid,'onset\tduration\ttrial_type\tresponse_time\n');
        for t = 1:length(onset_decision)
            
            % Partner is Friend=3, Stranger=2, Computer=1
            % Feedback is Reward=3, Neutral=2, Punishment=1
            
            %fprintf(fid,'onset\tduration\ttrial_type\tresponse_time\n');
            if     (feedback(t) == 1) && (Partner(t) == 1)
                trial_type = 'computer_punish';
            elseif (feedback(t) == 1) && (Partner(t) == 2)
                trial_type = 'stranger_punish';
            elseif (feedback(t) == 2) && (Partner(t) == 1)
                trial_type = 'computer_neutral';
            elseif (feedback(t) == 2) && (Partner(t) == 2)
                trial_type = 'stranger_neutral';
            elseif (feedback(t) == 3) && (Partner(t) == 1)
                trial_type = 'computer_reward';
            elseif (feedback(t) == 3) && (Partner(t) == 2)
                trial_type = 'stranger_reward';
            end
            
            
            
            if response(t) == 0 %missed response
                fprintf(fid,'%f\t%f\t%s\t%s\n',onset_decision(t),2.8,'miss_decision','n/a'); % max duration with outcome as #
                fprintf(fid,'%f\t%f\t%s\t%s\n',onset_outcome(t),duration(t),'miss_outcome','n/a'); % outcome is just #
            else
                % Ori: Right index is 2, left index is 7
                if Partner(t) == 1 % computer
                    if response(t) == 2
                        fprintf(fid,'%f\t%f\t%s\t%f\n',onset_decision(t),RT(t),'guess_rightButton_computer',RT(t));
                    elseif response(t) == 7
                        fprintf(fid,'%f\t%f\t%s\t%f\n',onset_decision(t),RT(t),'guess_leftButton_computer',RT(t));
                    end
                elseif Partner(t) == 2 % stranger (face)
                    if response(t) == 2
                        fprintf(fid,'%f\t%f\t%s\t%f\n',onset_decision(t),RT(t),'guess_rightButton_face',RT(t));
                    elseif response(t) == 7
                        fprintf(fid,'%f\t%f\t%s\t%f\n',onset_decision(t),RT(t),'guess_leftButton_face',RT(t));
                    end
                end
                fprintf(fid,'%f\t%f\t%s\t%s\n',onset_outcome(t),duration(t),['outcome_' trial_type],'n/a');
            end
            
            
        end
        fclose(fid);
    end
    cd(codedir);
    
catch ME
    disp(ME.message)
    disp(['check line: ' num2str(ME.stack(1).line) ]);
    %keyboard
end

