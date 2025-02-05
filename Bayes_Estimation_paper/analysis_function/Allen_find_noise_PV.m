

% obtain the manifest table of experimental sessions
sessions = bot.listSessions('VisualCoding', 'Ophys');

% select experiments from the most orientation/direction-selective areas of visual cortex 
sessions = sessions(ismember(sessions.targeted_structure_acronym, categorical(["VISp"])), :);

% select experiments with GCaMP6f expression enriched in Layer 2/3 and Layer 4 of cortex
sessions = sessions(sessions.cre_line == "Pvalb-IRES-Cre", :);

driftGratings = {};
varResponse = {};
meanResponse = {};
sessSuccess = {};    

for i=1:size(sessions,1)
    disp(['Working on ' int2str(i) ' of ' int2str(size(sessions,1))]);
    driftGratings{i} = [];
    varResponse{i} = [];
    meanResponse{i} = [];
    sessSuccess{i} = [];    
    try,
        try,
            sess = bot.getSessions(sessions(i,:));
        catch,
            sessSuccess{i} = lasterr;
        end;
        driftGratings{i} = sess.getStimulusTable('drifting_gratings');
    end

    if ~isempty(driftGratings{i}),
        traces = sess.fluorescence_traces_dff;

        % use a utility function to estimate single-trial responses for each stimulus
        responses = bot.util.StimulusAlignedResp(driftGratings{i}, traces);
        
        % locate "blank" stimulus presentations
        blankPresentations = driftGratings{i}.blank_sweep == 1;
        
        % extract neural responses to blank presentations
        blankResponses = mean(responses(blankPresentations, :),'omitnan');
        
        % subtract blank responses
        responsesCorrected = responses - blankResponses;

        % get a list of unique stimulus parameters and assign a stimulus ID to each presentation
        driftingGratingsTbl = timetable2table(driftGratings{i},"ConvertRowTimes",false);
        [uniqueStimuli, ~, presentationStimID] = unique(driftingGratingsTbl(:, 1:2), 'rows');
        
        % find NAN stimuli and remove them
        nanStimuli = isnan(uniqueStimuli.temporal_frequency);
        uniqueStimuli = uniqueStimuli(~nanStimuli, :);
        numStimuli = size(uniqueStimuli, 1);
        numROIs = size(responsesCorrected, 2);
        
        % preallocate arrays
        singleTrialResponses = cell(numStimuli, 1);
        meanStimulusResponses = nan(numStimuli, numROIs);
        varStimulusResponses = nan(numStimuli, numROIs);
        
        % collect single-trial responses for each stimulus
        for stimulusID = 1:size(uniqueStimuli, 1)
            % find the presentations of this stimulus
            theseTrials = presentationStimID == stimulusID;
            
            % collect all the fluorescence responses for this stimulus into a cell array
            singleTrialResponses{stimulusID} = responsesCorrected(theseTrials, :);
            
            % estimate response statistics for this stimulus
            meanStimulusResponses(stimulusID, :) = mean(singleTrialResponses{stimulusID},'omitnan');
            varStimulusResponses(stimulusID, :) = var(singleTrialResponses{stimulusID},'omitnan');
        end
        varResponse{i} = varStimulusResponses;
        meanResponse{i} = meanStimulusResponses;
    end;

end;

PV_varResponse = varResponse;
PV_meanResponse = meanResponse;