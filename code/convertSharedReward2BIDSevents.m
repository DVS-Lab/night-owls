function out = convertSharedReward2BIDSevents(subj,ses)
% This function converts the raw behavioral output from psychopy into
% the BIDS *_events.tsv file format. It also collects summary information
% about the subject's data in the "out" variable.

% Example convertSharedReward2BIDSevents(101,1)

try

    % set up paths
    scriptname = matlab.desktop.editor.getActiveFilename;
    [codedir,~,~] = fileparts(scriptname);
    [dsdir,~,~] = fileparts(codedir);

    % make default output
    out.ntrials(1) = 0;
    out.ntrials(2) = 0;
    out.nmisses(1) = 0;
    out.nmisses(2) = 0;
    out.nfiles = 0;

    % get relative path for source data. repos should be in same dir
    logdir = fullfile(dsdir,'stimuli','sharedreward','logs');

    for r = 1:2
        fname = fullfile(logdir,['sub-' num2str(subj)],sprintf('sub-%03d_task-sharedreward_ses-%d_run-%d_raw.csv',subj,ses,r));

        if r == 1 % only needed for first pass through
            [sublogdir,~,~] = fileparts(fname);
            sublogdir=convertStringsToChars(sublogdir);
            nfiles = dir([sublogdir '/*.csv']);
            out.nfiles = length(nfiles);
        end

        if exist(fname,'file')
            T = readtable(fname,'TreatAsEmpty','--');
        else
            fprintf('sub-%d_ses-%d_task-sharedreward_run-%d: No data found. Exiting...\n', subj, r)
        end

        % strip out irrelevant information and missed trials
        T = T(:,{'rt','decision_onset','outcome_onset','InitFixOnset','outcome_offset','Feedback','Partner','resp'});
        goodtrials =  ~isnan(T.resp);
        T = T(goodtrials,:);

        if height(T) < 54
            fprintf('incomplete data for sub-%d_ses-%d_run-%d\n', subj, ses, r)
        end

        onset_decision = T.decision_onset;
        onset_outcome = T.outcome_onset;
        duration = T.outcome_offset - T.outcome_onset; % outcome
        RT = T.rt;
        Partner = T.Partner;
        feedback = T.Feedback;
        response = T.resp;

        out.ntrials(r) = height(T);
        out.nmisses(r) = sum(T.resp < 1);

        % output file
        fname = sprintf('sub-%03d_ses-%02d_task-sharedreward_run-%d_events.tsv',subj,ses,r);
        output = fullfile(dsdir,'bids',['sub-' num2str(subj)],sprintf('ses-%02d',ses),'func');
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
    keyboard
end

