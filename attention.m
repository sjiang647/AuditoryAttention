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

%% Load auditory noise stimuli

cd('AudioStimuli');
names = dir('*.wav');
audios = cell(1, length(names));
for i = 1:length(names)
    audios{i} = audioread(names(i).name);
    audios{i} = audios{i} * .5;
end
cd('..');

%% Constants and global variables

% Experiment
numTrials = 11;
tonePause = 0.300;
trialPause = 0.500;

% Auditory tone generation
numTones = 6;
meanRange = 48:72;
toneRange = [1 3 5];
testRange = [2 4 6];

% Data saving
data = zeros(1, numTrials);
% 1. first name, 2. last name, 3. gender, 4. age [1x1]
% 5. mean tone, 6. noise type, 7. outlier tone, 8. outlier position [1xnumTrials]
subjectData = cell(1, 6);

%% Counterbalance conditions

sequence = randperm(numTrials);
askWhat = mod(sequence, 2); % 1 if mean, 0 if word
highLow = floor(mod(sequence, 4) ./ 2); % 1 if high, 0 if low
focusWhat = floor(mod(sequence, 8) ./ 4); % 1 if mean, 0 if word

testDist = testRange(floor(mod(sequence, 24) ./ 8) + 1); % 1, 2, or 3
testDist = testDist .* round(2 * (highLow - 0.5));

counterbalancing = [askWhat; testDist; focusWhat];
subjectData{5} = counterbalancing;

%% Subject data input

subjectData{1} = Ask(window, 'First Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{2} = Ask(window, 'Last Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{3} = Ask(window, 'Gender(M/F): ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{4} = str2double(Ask(window, 'Age: ', [],[], 'GetChar', RectLeft, RectTop, 25));

%% Task instructions

Screen('DrawText', window, 'You will listen to various audio tones. Pay attention to the various auditory stimuli.', center(1)-420, center(2)-60);
Screen('DrawText', window, 'At the end of each trial, you will be asked to make an input based on a question asked.', center(1)-450, center(2)-30);
Screen('DrawText', window, 'Press "Return" to continue', center(1)-150, center(2));
Screen('Flip', window);
KbWait([], 2);

%% Stimuli display (experiment)

handle = PsychPortAudio('Open', [], [], 0, 44100, 2);

for trial = 1:numTrials
    Screen('Flip', window);
    % 1. what to ask 2. test tone distance % 3. what to focus on
    trialSettings = counterbalancing(:, trial);
    
    meanTone = randsample(meanRange, 1); % Randomly shuffle tones to be played
    tones = randsample([-toneRange toneRange], numTones);
    setSounds = randsample(length(audios), numTones, true); % Creating set of sounds

    if trial < 4
        %% Only tones

        Screen('DrawText', window, 'Focus on the tones.', center(1) - 150, center(2));
        Screen('Flip', window);
        KbWait();
        for toneNum = 1:numTones
            playAudio(tones(toneNum) + meanTone);
            WaitSecs(0.45);
        end
        
        Screen('DrawText', window, 'You will now hear a test tone.', center(1) - 250, center(2) - 25);
        Screen('DrawText', window, 'Press "Return" to continue.', center(1)- 250, center(2));
        Screen('Flip', window);
        KbWait();

        playAudio(meanTone + trialSettings(2));
        
        % Keyboard instructions
        Screen('DrawText', window, 'Press h if the test tone was higher than the mean.', center(1) - 250, center(2) - 25);
        Screen('DrawText', window, 'Press l if the test tone was lower than the mean.', center(1) - 250, center(2));
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
    elseif trial < 7
        %% Only words
        
        Screen('DrawText', window, 'Focus on the words.', center(1) - 150, center(2));
        Screen('Flip', window);
        KbWait();
        for toneNum = 1:numTones
            playAudio(audios{setSounds(toneNum)});
            WaitSecs(.3);
        end
        
        % Ask for number of times words played
        while true
            res = Ask(window, ['How  many times was ' names(randsample(3,1)).name ' played (1-6): '], [],[], 'GetChar', RectLeft, RectTop, 25);
            if res=='1' || res=='2' || res=='3' || res=='4' || res=='5' || res=='6' || res=='7' || res=='8' || res=='9'
                break;
            end
        end
    else 
        %% Main experiment
        
        % Display instructions
        if trialSettings(1)
            Screen('DrawText', window, 'Focus on the tones.', center(1) - 150, center(2));
            Screen('Flip', window);
            KbWait();
        else
            Screen('DrawText', window, 'Focus on the words.', center(1) - 150, center(2));
            Screen('Flip', window);
            KbWait();
        end
        
        % Loop through and play all tones
        for toneNum = 1:numTones
            playAudio(audios{setSounds(toneNum)});
            WaitSecs(0.15);
            playAudio(tones(toneNum) + meanTone);
            WaitSecs(0.45);
        end
        
        if trialSettings(3)
            % Audio task instructions
            Screen('DrawText', window, 'You will now hear a test tone.', center(1) - 250, center(2) - 25);
            Screen('DrawText', window, 'Press "Return" to continue.', center(1)- 250, center(2));
            Screen('Flip', window);
            KbWait();
            playAudio(meanTone + trialSettings(2));
            
            % Keyboard instructions
            Screen('DrawText', window, 'Press h if the test tone was higher than the mean.', center(1) - 250, center(2) - 25);
            Screen('DrawText', window, 'Press l if the test tone was lower than the mean.', center(1) - 250, center(2));
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
            % Ask for number of times words played
            while true
                res = Ask(window, ['How  many times was ' names(randsample(3,1)).name ' played (1-6): '], [],[], 'GetChar', RectLeft, RectTop, 25);
                if res=='1' || res=='2' || res=='3' || res=='4' || res=='5' || res=='6' 
                    break;
                end
            end
        end
    end
    WaitSecs(trialPause);
end

%% End experiment

PsychPortAudio('Close', handle);
ShowCursor();
Screen('CloseAll');

%% Save data

if ~isdir(['participant_data/', subjectData{1}])
    mkdir(['participant_data/', subjectData{1}]);
end

%% Functions

function playAudio(m)
    handle = PsychPortAudio('Open', [], [], 0, 44100, 2);
    
    fs = 44100;
    toneLength = 0:1/fs:.450;
    freqRamp = 1/(2*.01);
    rampVector = 1:441;
    
    offset = (1 + sin(2*pi * freqRamp * rampVector ./ fs + (pi/2))) / 2;
    onset = (1 + sin(2*pi * freqRamp * rampVector ./ fs + (-pi/2))) / 2;
    
    if ~isscalar(m)
        if size(m, 2) < size(m, 1)
            m = m';
        end
        newm = repmat(m,2,1);
    else
        toneFrequency = 440 * 2^((m - 69)/12);
        midiTone = sin(2*pi * toneFrequency * toneLength);
        midiTone(1:441) = onset .* midiTone(1:441);
        midiTone(end - 440: end) = offset .* midiTone(end - 440: end);
        newm = repmat(midiTone, 2, 1); 
    end
    
    PsychPortAudio('FillBuffer', handle, newm);
    PsychPortAudio('Start', handle, 1, 0, 1);
end
