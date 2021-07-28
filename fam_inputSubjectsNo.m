% @date 200101 @author D. Iggena
% @update 210527 

% function asks to enter two strings containing numbes. 
% The first subject-number and the last
% subject-number of your dataset. 
% Required too loop through data-folders. 

% examples
% input: 33001 33999 (string)
% output: firstSubject, lastSubject (numbers)

function [firstSubject,lastSubject]=fam_inputSubjectsNo()


% 1st subject
invalidSubjectStartNo = true;
while invalidSubjectStartNo
    firstSubject = str2double(input('Enter first subjectNo: ','s'));
    if ~isreal(firstSubject)
        clear('firstSubject');
    else
        invalidSubjectStartNo = false;
    end
end

% last subject
invalidSubjectEndNo = true;
while invalidSubjectEndNo
    lastSubject = str2double(input('Enter last subjectNo: ','s'));
    if ~isreal(firstSubject)
        clear('firstSubject');
    else
        invalidSubjectEndNo = false;
    end
end

end

