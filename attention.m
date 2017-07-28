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
audios = cell(length(names));
for i = 1:length(names)
    audios{i} = audioread(names(i).name);
    audios{i} = audios{i} * .7;
end
cd('..');

%% Constants and global variables

% Experiment
numTrials = 200;
numTests = 4;
tonePause = 0.350;
wordPause = 0.300;
trialPause = 0.500;

% Auditory tone generation
numTones = 4;
meanRange = 48:72;
toneRange = [1 3 5];
testRange = [2 4 6];

% Data saving
data = zeros(1, numTrials);
% 1. first name, 2. last name, 3. gender, 4. age, 5. counterbalancing,
% 6. tone accuracy, 7. names asked, 8. sounds played, 9. word responses
subjectData = cell(1, 9);

%% Generate conditions

% Generate counterbalanced conditions

sequence = randperm(numTrials - 2 * numTests);
askWhat = [randi([0, 1], 1, 2 * numTests) mod(sequence, 2)]; % 1 if mean, 0 if word
highLow = [randi([0, 1], 1, 2 * numTests) floor(mod(sequence, 4) ./ 2)]; % 1 if high, 0 if low
focusWhat = [randi([0, 1], 1, 2 * numTests) floor(mod(sequence, 8) ./ 4)]; % 1 if mean, 0 if word

testDist = [testRange(randi([1, 3], 1, 2 * numTests)) testRange(floor(mod(sequence, 24) ./ 8) + 1)]; % 1, 2, or 3
testDist = testDist .* round(2 * (highLow - 0.5));

counterbalancing = [askWhat; testDist; focusWhat];
subjectData{5} = counterbalancing;

% Generate random conditions

nameIndices = randi(length(audios), 1, numTrials);
setSounds = zeros(numTrials, numTones);
for i = 1:numTrials
    setSounds(i,:) = randsample(length(audios), numTones, true);
end

subjectData{7} = nameIndices;
subjectData{8} = setSounds;

%% Subject data input

% Input general info
subjectData{1} = Ask(window, 'First Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{2} = Ask(window, 'Last Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{3} = Ask(window, 'Gender(M/F): ', [],[], 'GetChar', RectLeft, RectTop, 25);
subjectData{4} = str2double(Ask(window, 'Age: ', [],[], 'GetChar', RectLeft, RectTop, 25));

% Initialize response cells
subjectData{6} = repmat(-1, 1, numTrials);
subjectData{9} = repmat(-1, 1, numTrials);

%% Task instructions

Screen('DrawText', window, 'You will listen to various audio tones. Pay attention to the various auditory stimuli.', center(1)-windowX/3.5, center(2));
Screen('DrawText', window, 'At the end of each trial, you will be asked to make an input based on a question asked.', center(1)-windowX/3, center(2)+windowX/35);
Screen('DrawText', window, 'Press enter to continue', center(1)-windowX/13.5, center(2) + windowX/13);
Screen('Flip', window);
KbWait([], 2);

%% Stimuli display (experiment)

handle = PsychPortAudio('Open', [], [], 0, 44100, 2);

for trial = 1:numTrials
    Screen('Flip', window);
    
    if trial <= numTests
        %% Audio only
        meanTone = randsample(meanRange, 1);
        tones = randsample([-toneRange toneRange], numTones);

        showSingleInstructions(window, numTones, rect, 'audio tones');
        
        for toneNum = 1:numTones
            playAudio(tones(toneNum) + meanTone);
            WaitSecs(tonePause);
        end
        
        audioTaskInstructions(window, rect, meanTone + testDist(trial));
        data(trial) = analyzeHighLow(testDist(trial));
    elseif trial <= numTests * 2
        %% Words only
        showSingleInstructions(window, numTones, rect, 'words');
        
        for toneNum = 1:numTones
            playAudio(audios{setSounds(trial, toneNum)});
            WaitSecs(tonePause);
        end
        
        % Ask for number of times words played
        res = wordTaskInstructions(window, names(nameIndices(trial)).name, numTones);
        subjectData{9}(trial) = str2double(res);
    else
        %% Main Experiment
        showSingleInstructions(window, numTones, rect, 'sets of words and tones');
        
        % Randomly shuffle tones to be played
        meanTone = randsample(meanRange, 1);
        tones = randsample([-toneRange toneRange], numTones);
        
        % Display instructions
        if focusWhat(trial)
            showFocusInstructions(window, rect, 'Focus on the tones.');
        else
            showFocusInstructions(window, rect, 'Focus on the words.');
        end
        
        % Loop through and play all tones
        for toneNum = 1:numTones
            playAudio(audios{setSounds(trial, toneNum)});
            WaitSecs(wordPause);
            playAudio(tones(toneNum) + meanTone);
            WaitSecs(tonePause);
        end
        
        if askWhat(trial) % Play test tone
            audioTaskInstructions(window, rect, meanTone + testDist(trial));
            subjectData{6}(trial) = analyzeHighLow(testDist(trial));
        else % Ask for number of times words played
            res = wordTaskInstructions(window, names(nameIndices(trial)).name, numTones);
            subjectData{9}(trial) = str2double(res);
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

cd(['participant_data/', subjectData{1}]);
save('data', 'subjectData');
cd('..');
cd('..');

%% Info Instructions

function showSingleInstructions(window, numTones, rect, type)
    msg = [num2str(numTones) ' ' type ' will be played.'];
    showFocusInstructions(window, rect, msg);
end

function showFocusInstructions(window, rect, msg)
    windowX = rect(3);
    windowY = rect(4);
    center = [windowX/2, windowY/2];
    
    Screen('DrawText', window, msg, center(1) - windowX/11, center(2)); 
    Screen('DrawText', window, 'Press ENTER to continue.', center(1) - windowX/11, center(2)+windowY/13);
    Screen('Flip', window);
    KbWait();
    Screen('DrawText', window, msg, center(1) - windowX/11, center(2));
    Screen('Flip', window);
end

%% Task Instructions

function audioTaskInstructions(window, rect, toneToPlay)
    windowX = rect(3);
    windowY = rect(4);
    center = [windowX/2, windowY/2];
    
    % Audio task instructions
    Screen('DrawText', window, 'You will now hear a test tone.', center(1) - windowX/12, center(2));
    Screen('DrawText', window, 'Press ENTER to continue.', center(1) - windowX/11, center(2) + windowY/13);
    Screen('Flip', window);
    KbWait();
    playAudio(toneToPlay);

    % Keyboard instructions
    Screen('DrawText', window, 'Press H if the test tone was higher than the mean.', center(1) - windowX/5, center(2) - windowY/13);
    Screen('DrawText', window, 'Press L if the test tone was lower than the mean.', center(1) - windowX/5.1, center(2));
    Screen('Flip', window);
end

function response = wordTaskInstructions(window, nameToAsk, numTones)
    while true
        res = Ask(window, ['How  many times was ' nameToAsk ' played (0-' num2str(numTones) '): '], [],[], 'GetChar', RectLeft, RectTop, 25);
        if ismember(str2double(res), 0:numTones)
            break;
        else 
            res = Ask(window, ['How  many times was ' nameToAsk ' played (0-' num2str(numTones) '): '], [],[], 'GetChar', RectLeft, RectTop, 25);
        end 
    end
    response = res;
end

%% Response Analysis

function correct = analyzeHighLow(position)
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
    if (response == 'h' && position > 0) || (response == 'l' && position < 0)
        correct = 1;
    else
        correct = 0;
    end
end

%% Audio playback

function playAudio(m)
    handle = PsychPortAudio('Open', [], [], 0, 44100, 2);
    
    fs = 44100;
    toneLength = 0:1/fs:.350;
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

