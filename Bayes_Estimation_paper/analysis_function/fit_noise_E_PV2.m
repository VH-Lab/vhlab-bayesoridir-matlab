% obtain all sessions in the ephys dataset
sessions = bot.listSessions('VisualCoding', 'Ephys');

% obtain all sessions with the desired session_type and then select the first one
desiredSessions = sessions(sessions.session_type == 'brain_observatory_1.1', :);
sess = bot.getSessions(desiredSessions(1, :));


units = sess.units;


% obtain all probes in the session
probes = sess.probes;

% obtain the brain region that each probe examines
disp(probes.ephys_structure_acronyms)

% the stimulus_epochs property displays the sequence of visual stimuli
sess.stimulus_epochs

probes.ephys_structure_acronyms(1)
probe = bot.getProbes(probes(1, :));

% select all units of the probe
probeUnits = probe.units;

% select all units in the AM using logical indexing
probeUnits = probeUnits(probeUnits.ephys_structure_acronym == 'VISam', :);

% obtain a `unit` object for the first unit in the AM
unit = bot.getUnits(probeUnits(1, :));


% get a full table of individual stimulus presentations
stimuli = sess.stimulus_presentations;

