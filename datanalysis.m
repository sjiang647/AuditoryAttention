% general outlines of autastic in 2017 2.0

%% Psychometric curve analysis
 %subjectData{1} = Ask(window, 'First Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
% subjectData{2} = Ask(window, 'Last Name: ', [],[], 'GetChar', RectLeft, RectTop, 25);
% subjectData{3} = Ask(window, 'Gender(M/F): ', [],[], 'GetChar', RectLeft, RectTop, 25);
% subjectData{4} = str2double(Ask(window, 'Age: ', [],[], 'GetChar', RectLeft, RectTop, 25));
%   counterbalancing = [askWhat; testDist; focusWhat];
%   askWhat = [randi([0, 1], 1, 2 * numTests) mod(sequence, 2)]; % 1 if mean, 0 if word
%   highLow = [randi([0, 1], 1, 2 * numTests) floor(mod(sequence, 4) ./ 2)]; % 1 if high, 0 if low
%   focusWhat = [randi([0, 1], 1, 2 * numTests) floor(mod(sequence, 8) ./ 4)]; % 1 if mean, 0 if word
% subjectData{5} = counterbalancing;
% subjectData{7} = nameIndices;
% subjectData{8} = setSounds;
% subjectData{9} = repmat(-1, 1, numTrials);



for i = 1:length(names)
    %% 1. Cleaning data
    
    % Iteratively call in individual subject data
    
    % Create a matrix for each data
    celerey{i} =  load(['subject_results/data_' num2str(i) '.mat']);
    
    % Organize data so that outlier distance accounts for both +/-
    % subjectData{4}(1,:) is outlier offset
    counterbalancing = celerey{i}.subjectData{5};
    
    % Calculate accuracy
    % subjectData{5} is right/wrong
    accuracy = celerey{i}.subjectData{5};
    
    % Create 'all_data' matrix that combines all data
    all_data = [counterbalancing; accuracy];
    
    % Apply certain row/column to j_fit to compare
    
    % ??????
    
    %% 2. Flipping data
    
    % Since we want to measure % that outlier is higher than mean, we
    % need to flip data for negative values.
    
    % (Hit = 1) in negative outliers are saying ?lower than mean? so we
    % want to flip that. Vice versa for 0s.
    
    % So, for negative outlier distances, we want to flip accuracies of 0s
    % to 1s and 1s to 0s
    
    for j = 1:length(all_data)
        if all_data(1, j) < 0
            all_data(3, j) = 1 - all_data(3, j);
        end
    end
    
    outlier_diffs = [-16 -14 -10 -6 6 10 14 16];
    accuracy_percentage = zeros(1, 8);
    for j = 1:length(outlier_diffs)
        indices = find(all_data(1,:) == outlier_diffs(j));
        results = all_data(3, indices);
        accuracy_percentage(j) = mean(results);
        if outlier_diffs(j) == 14
            sendAccuracies{i} = results;
        end
    end
    
%     for thing = 1:length(sendAccuracies)
%         sendAccuracies{thing} = sendAccuracies{thing};
%         
%     end
    %% 3. Calling jfit
    
    % Make sure j_fit.m is in same folder as your analysis.m
    % Call in j_fit within your analysis.m code
    
    % ** You do not have to directly make changes on j_fit.m file
    
    
    [a_cond1, b_cond1] = j_fit(all_data(1,:)', all_data(3,:)','logistic1',2);
    
    sendPvalues(i) = b_cond1;
    
    sendAll = {sendPvalues; sendAccuracies};
    if ~isdir(['Group5Send/', names{i}])
        mkdir(['Group5Send/', names{i}]);
    end
    
    cd(['Group5Send/', names{i}]);
    save('data', 'sendAll');
    cd ..
    cd ..
end

