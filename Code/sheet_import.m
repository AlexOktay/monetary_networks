function [data] = sheet_import(filename)
%SHEET_IMPORT Imports a multi-sheet excel file

% Get sheet names
[~,sheet_name]=xlsfinfo(filename);

% Import each sheet and place it in a vector
for k=1:numel(sheet_name)
    data{k}=xlsread(filename,sheet_name{k});
end

end