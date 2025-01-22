 % Read the image
image = imread('zap000.jpg');  % Replace 'your_image.jpg' with the actual image file and the image is 0.5x of binning4 tomo

% Display the image
imshow(image);

% Allow the user to draw a freehand line on the image
h = imfreehand;
wait(h); % Wait for the user to finish drawing

% Get the coordinates of the freehand line
freehandLine = h.getPosition;

% Adjust y-coordinates to reverse the direction
freehandLine(:, 2) = size(image, 1) - freehandLine(:, 2) + 1;

% Read coordinates from the text file (assuming the file is 'coordinates.txt')
coordinates = readmatrix('p04_ori.xlsx');

% Find the nearest point on the freehand line for each point in the text file
nearestPoints = zeros(size(coordinates, 1), 2);
distances = zeros(size(coordinates, 1), 1);

for i = 1:size(coordinates, 1)
    % Extract x and y coordinates from the 17th and 23rd columns and divided by 2 
    x = coordinates(i, 1)/2;
    y = coordinates(i, 2)/2;
    
    % Find the nearest point on the freehand line
    [Idx(i),distances(i)] = knnsearch(freehandLine, [x, y]); 
end    

% Calculate the distance along the line from the start point to the nearest point
for n = 1:length(Idx)
    distanceAlongLine(n) = sum(sqrt(sum(diff(freehandLine(1:Idx(n), :)).^2, 2)));
end

for k =1:size(coordinates, 1)  
Coorts_receptor(k,:)= [2*distanceAlongLine(k), coordinates(k,3)]; % distancealongline mutiplied by 2 to go back to binning 4 tomo
end

% Export Coorts_receptor to an Excel file
writematrix(Coorts_receptor, 'Coorts_receptor.xlsx');
disp('Coorts_receptor exported to Coorts_receptor.xlsx');
