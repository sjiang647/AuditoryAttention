
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

window_w = rect(3);
window_h = rect(4);
center_x = window_w / 2;
center_y = window_h / 2;
cd AudioStimuli
names = dir('*.wav');
audios = cell(length(names));
for i = 1:length(names)
    audios{i} = audioread(names(i).name);
end
cd ..;
% toneLength = 0:1/44100:.300;
% freqRamp = 1/(2*(.01));
% rampVector = [1:441];
% fs = 44100;
% offset = (1+sin(2*pi*freqRamp*rampVector./fs + (pi/2)))/2;
% onset = (1+sin(2*pi*freqRamp*rampVector./fs + (-pi/2)))/2;
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

askWhat = mod(randperm(numTrials), 2); % 1 if mean, 0 if word
highLow = mod(randperm(numTrials), 2); % 1 if high, 0 if low
focusWhat = mod(randperm(numTrials), 2); % 1 if mean, 0 if word

for i = 1:numTrials
    counterbalancing(1,i) = askWhat(i);
    counterbalancing(2,i) = highLow(i);
    counterbalancing(3,i) = focusWhat(i);
end

subjectData{5} = counterbalancing;

%% Subject data input

subjectData{1} = Ask(window, 'First Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{2} = Ask(window, 'Last Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{3} = Ask(window, 'Gender(M/F): ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{4} = str2double(Ask(window, 'Age: ', [],[], 'GetChar', RectLeft, RectTop, 25));

%% Task instructions

%% Stimuli display (experiment)

handle = PsychPortAudio('Open', [], [], 0, 44100, 2); 

for trial = 1:numTrial
    Screen('Flip', window);
    
    meanDiff = 2;
    
    % 1. what to ask 2. high/low % 
    trialSettings = counterbalancing(:, trial);
    
    % Randomly shuffle tones to be played
    meanTone = randsample(meanRange, 1);
    tones = randsample([-toneRange toneRange], numTones); 
    toneVectors = freq(allTones + tones);
    
    % Display instructions
    if trialSettings(1)
        Screen('DrawText', window, 'Focus to the 6 tones.', center_x - 150, center_y);
    else
        Screen('DrawText', window, 'Focus on the words.', center_x - 150, center_y);
    end
    
    % Loop through and play all tones
    for toneNum = 1:numTones
        PsychPortAudio('FillBuffer', handle, toneVectors{toneNum});
        PsychPortAudio('Start', handle, 1, 0, 1);
        WaitSecs(tonePause);
        PsychPortAudio('Stop', handle);
    end
    
    if trialSettings(3)
        % Audio task instructions
        Screen('DrawText', window, 'You will now hear a test tone.', center_x - 250, center_y - 25);
        Screen('DrawText', window, 'Press any key to continue.', center_x - 250, center_y);
        Screen('Flip', window);
        
        % Play audio tone
        offtone = meanTone + meanDiff * round((counterbalancing(2) - 0.5) * 2);
        PsychPortAudio('FillBuffer', handle, toneVectors{toneNum});
        PsychPortAudio('Start', handle, 1, 0, 1);
        WaitSecs(tonePause);
        PsychPortAudio('Stop', handle);
        
        % Keyboard instructions
        Screen('DrawText', window, 'Press h if the test tone was higher than the mean.', center_x - 250, center_y - 25);
        Screen('DrawText', window, 'Press l if the test tone was lower than the mean.', center_x - 250, center_y);
        Screen('Flip', window);
        
        % Check keyboard presses
        KbName('UnifyKeyNames');
        while true
            [keyDown, secs, keyCode, deltaSecs] = KbCheck(-1); 
            key = KbName(find(keyCode));

            if strcmp(key, 'h')
                response = 'h';
                break;
            end
            if strcmp(key, 'l')
                response = 'l';
                break;
            end
        end
        % Check accuracy of response
        if (response == 'h' && trialSettings(2)) || (response == 'l' && ~trialSettings(2))
            data(trial) = 1;
        end
    else
        % Ask for words
        
    end

    WaitSecs(trialPause);
end

%% End experiment

PsychPortAudio('Close', handle);
ShowCursor();
Screen('CloseAll');

%% Data analysis

%% Save data

if ~isdir(['participant_data/', subjectData{1}])
    mkdir(['participant_data/', subjectData{1}]);
end

function playAudio(m)
handle = PsychPortAudio('Open', [], [], 0, 44100, 2); 

toneLength = 0:1/44100:.300;
freqRamp = 1/(2*(.01));
rampVector = [1:441];
fs = 44100;
offset = (1+sin(2*pi*freqRamp*rampVector./fs + (pi/2)))/2;
onset = (1+sin(2*pi*freqRamp*rampVector./fs + (-pi/2)))/2;
    if ~isscalar(m)
        if size(m, 2) < size(m, 1)
            m = m';
        end
        newm = repmat(m,2,1);
        PsychPortAudio('FillBuffer', handle, newm);
        PsychPortAudio('Start', handle, 1, 0, 1);
        WaitSecs(.3)
        PsychPortAudio('Stop', handle);
    else
        toneFrequency = 440*2^((m-69)/12);
        midiTone = sin(2*pi* toneFrequency * toneLength);%creating the tones in terms of frequency
        midiTone(1:441) = onset .* midiTone(1:441);
        midiTone(end - 440: end) = offset .* midiTone(end - 440: end);
        newm = repmat(midiTone, 2, 1); %duplicates the sound in order to hear through headphones
        PsychPortAudio('FillBuffer', handle, newm);
        PsychPortAudio('Start', handle, 1, 0, 1);
        WaitSecs(.3)
        PsychPortAudio('Stop', handle);
    end

end