%% general outlines of autastic in 2017 2.0

clear all;
close all;
clc;

% 1. first name, 2. last name, 3. gender, 4. age, 5. counterbalancing,
% 6. tone accuracy, 7. names asked, 8. sounds played, 9. word responses

%% Psychometric curve analysis

celerey = cell(1, 5);
paramatrixa = zeros(5, 2);
paramatrixb = zeros(5, 2);

for i = 1:4
    % Create a matrix for each data
    celerey{i} =  load(['new_results/data_' num2str(i) '.mat']);
    celerey{i}.subjectData
    
    % Organize data so that outlier distance accounts for both +/-
    asked = celerey{i}.subjectData{5}(1,:);
    distanced = celerey{i}.subjectData{5}(2,:);
    focused = celerey{i}.subjectData{5}(3,:);
    responded = celerey{i}.subjectData{6};

    testPos = zeros(2, 48); % row 1 = focus, row 2 = no focus
    response = zeros(2, 48); % row 1 = focus, row 2 = no focus
    focusIndex = 1;
    ignoreIndex = 1;
    
    for thing = 9:200
       if asked(thing) % if tone was asked
           if focused(thing) % if subject focused on tone
               testPos(1, focusIndex) = distanced(thing);
               if distanced(thing) > 0 
                   response(1, focusIndex) = responded(thing);
               else % flip data if negative dist
                   response(1, focusIndex) = 1 - responded(thing);
               end
               focusIndex = focusIndex + 1;
           else % if subject focused on word
               testPos(2, ignoreIndex) = distanced(thing);
               if distanced(thing) > 0 
                   response(2, ignoreIndex) = responded(thing);
               else % flip data if negative dist
                   response(2, ignoreIndex) = 1 - responded(thing);
               end
               ignoreIndex = ignoreIndex + 1;
           end
       end
    end
    
    [a_cond1, b_cond1] = j_fit(testPos(1,:)', response(1,:)','logistic1',2);
    [a_cond2, b_cond2] = j_fit(testPos(2,:)', response(2,:)','logistic1',2);
    paramatrixa(i,:) = [a_cond1 a_cond2];
    paramatrixb(i,:) = [b_cond1 b_cond2];
end

%% Basic accuracy analysis

datarray = cell(1, 5);
% focus/ask: 1 = word/word, 2 = tone/word, 3 = word/tone, 4 = tone/tone
accuracies = zeros(4, 5);

for i = 1:4
    data = load(['new_results/data_' num2str(i) '.mat']);
    datarray{i} = data.subjectData;
    % Ignore first 8 trials
    for j = 9:200
        % subjectData{5} is the counterbalancing matrix
        asked = datarray{i}{5}(1, j);
        focused = datarray{i}{5}(3, j);
        index = (focused + 1) + 2 * asked;
        if asked % tone was played
            accuracies(index, i) = accuracies(index, i) + datarray{i}{6}(j);
        else % word was played
            num = sum(datarray{i}{8}(j,:) == datarray{i}{7}(j));
            accuracies(index, i) = accuracies(index, i) + (num == datarray{i}{9}(j));
        end
    end
end

% Calculate mean
accuracies = accuracies ./ 48;
