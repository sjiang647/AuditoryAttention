%% Objective: Test effects of attention on auditory ensemble perception

clear all;
close all;
clc;

%% Screen setup

Screen('Preference', 'SkipSyncTests', 1);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
[window, rect] = Screen('OpenWindow', 0); 
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
HideCursor();

windowX = rect(3);
windowY = rect(4);
center = [windowX/2, windowY/2];
%% Constants and global variables

% Experiment
numTrials = 5;
tonePause = 0.300;
trialPause = 0.500;

% Auditory tone generation
numTones = 7;
meanRange = 48:72;
toneRange = [2 4 6];
outlierRange = [6 8 10 12];

% Auditory frequency generation
fs = 44100;
toneDuration = 0.300;
toneLength = 0:1/fs:toneDuration;
freqRamp = 1/(2*.01);
rampVector = 1:441;

% Data saving
data = zeros(1, numTrials);
% 1. first name, 2. last name, 3. gender, 4. age [1x1]
% 5. mean tone, 6. noise type, 7. outlier tone, 8. outlier position [1xnumTrials]
subjectData = cell(1, 6);

%% Load/trim auditory noise stimuli

%% Generate auditory tone stimuli

offset = (1 + sin(2 * pi * freqRamp * rampVector ./ fs + (pi/2))) / 2;
onset = (1 + sin(2 * pi * freqRamp * rampVector ./ fs + (-pi/2))) / 2;
frequencies = cell(1, 127);

for k = 1:127 
    toneFrequency = 440 * 2 ^ ((k - 69)/12);
    midiTones = sin(2 * pi * toneFrequency * toneLength);
    midiTones(1:441) = onset .* midiTones(1:441); 
    midiTones((end - 440):end) = offset .* midiTones((end - 440):end);
    frequencies{k} = repmat(midiTones, 2, 1);
end

%% Counterbalance conditions

highLow = mod(randperm(numTrials), 2); % 1 if high, 0 if low
outlierDiff = outlierRange(mod(randperm(numTrials), 4) + 1);
outlierPos = mod(randperm(numTrials), 7) + 1;

for i = 1:numTrial
    if highlow(i) == 0
    outlierDiff(i) = -outlierDiff(i);
    end
end

counterbalancing = [outlierDiff; outlierPos];
subjectData{5} = counterbalancing;
% matrix of 3xNumTrials
%numtrials 1 column --> what to ask 
% 2nd column --> high or lower than th test tone 
% 3rd column --> what they focus on 

%% Subject data input

subjectData{1} = Ask(window, 'First Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{2} = Ask(window, 'Last Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{3} = Ask(window, 'Gender(M/F): ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{4} = str2double(Ask(window, 'Age: ', [],[], 'GetChar', RectLeft, RectTop, 25));

%% Task instructions

%% Stimuli display (experiment)

%% End experiment

ShowCursor();
Screen('CloseAll');

%% Data analysis

%% Save data

if ~isdir(['participant_data/', subjectData{1}])
    mkdir(['participant_data/', subjectData{1}]);
end

cd(['participant_data/', subjectData{1}]);
save('data', 'subjectData');
cd('..');
cd('..');
