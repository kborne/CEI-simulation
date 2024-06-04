function [data] = load_csv_data(directory_path)


% Enter the directory path that you want to loop through

% Get a list of all files in the directory
files = dir(fullfile(directory_path, '*.csv'));

% Initialize an empty array to store the data

for i = 1:15;
    data.CoM_mom{i} = zeros(length(files),3);
end

% Loop through each file in the directory
for i = 1:length(files)
    % Check if the current item is a file (not a directory)
    if ~files(i).isdir
        % Load the file using csvread
        file_data = csvread(fullfile(directory_path, files(i).name));
%         file_data = reshape(file_data,[3,15])';
        
    % Append the file data to the data array

    for j = 1:15    
        data.CoM_mom{j}(i,:)= file_data(j,:);

    end

    end
end

% Display the concatenated data array
disp(data);



end