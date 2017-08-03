%% general outlines of autastic in 2017 2.0

clear all;
close all;
clc;

% 1. first name, 2. last name, 3. gender, 4. age, 5. counterbalancing,
% 6. tone accuracy, 7. names asked, 8. sounds played, 9. word responses

%% Psychometric curve analysis

celerey = cell(5);

for i = 1:2
    %% 1. Cleaning data
    
    % Create a matrix for each data
    celerey{i} =  load(['new_results/data_' num2str(i) '.mat']);
    
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
    
%     sendPvalues(i) = b_cond1;
%     
%     sendAll = {sendPvalues; sendAccuracies};
%     if ~isdir(['Group5Send/', names{i}])
%         mkdir(['Group5Send/', names{i}]);
%     end
%     
%     cd(['Group5Send/', names{i}]);
%     save('data', 'sendAll');
%     cd ..
%     cd ..
end

