%Define number of batch processing
num = 18;


%%%%%

%Read Files
columnNum = 1;
for file=1:num
    XFile = file + "_trans" + ".xlsx";
    YFile = "PSD_" + file + ".xlsx";
    ZFile = file + "_bg" + ".txt";
    columnNum = processing(XFile, YFile, ZFile,file, columnNum);
end
"Done Processing"
%%%%%

%process files
function columnNum = processing(XFile,YFile,ZFile, file, columnNum)
    %Importing the vesical value data (X)
    Xx = readmatrix(XFile,'Range','A:A');
    %Xx(1,:) = [];
    Xy = readmatrix(XFile,'Range','B:B');
    %Xy(1,:) = [];
    X = [Xx, Xy];
    [numVesicles, y] = size(X);
    
    %Importing the gray value data (Y) (Possible Coding improvements)
    inputY = readmatrix (YFile);
    [nrows, ncols]=size(inputY); %stores the number of rows and columns in input
    values=[]; %output matrix
    for r=1:nrows
        for c=1:ncols
            values=[values;r c inputY(r,c)];% keep appending [r,c,input] to new columns
        end
    end
    

    %Import backgroud value
    bg_value = importdata(ZFile);
    TF = values(:,3)<bg_value;
    M3 = values(TF, :);
    M = M3(:,3);

    %Finding top 25% value threshold
    [m, n] = size(M);
    TTP = ceil(m*n*.25);%find the top 25 percent value
    [val,ind] = sort(M);
    Threshold = val(TTP);%the value of the top 25 percent
    
    %Clustring using DBSCAN
    M_Thres = values(values(:,3)<Threshold,1:2);
    idx = dbscan(M_Thres,6,32); 
    %gscatter(M_Thres(:,1),M_Thres(:,2),idx)

    %Delete points which belong to cluster -1 
    N = [M_Thres idx];
    idx2 = idx(idx~= -1);
    numGroups = length(unique(idx2));
    clr = hsv(numGroups);
    N1=N(N(:,3) ~= -1,:);
%     h = figure();
%     set(h,'Visible', 'off');
%     gscatter(N(N(:,3) ~= -1,1),N(N(:,3) ~= -1,2),idx2,clr);
% 
%     filename = strcat("DataSet",num2str(file));
%     axis image;
%     xlabel("X length");
%     ylabel("Y length");
%     title(filename);
%     filename = strcat(filename,".tif");
%     exportgraphics(h,filename,'Resolution',400)
     
    %Finding the centroids of the clusters and area of each cluster
    t = unique(idx2); % the unique values in the idx2
    maxDistances = zeros(length(t), 1); % initialize an array to store max distances
    for i = 1:length(t)
%     tf = N1(:,3)==i; %i is the No. of clusters.
%     p=N1(tf,:);
%     x=p(:,1);
%     y=p(:,2);
%     dt = delaunayTriangulation(x,y);
%     k = convexHull(dt);
%     area(i,1)=abs(trapz(x(k),y(k)));
  
 
    Y(length(t),2) = 0;%matrix holding mean value
        M5_Temp = N1(N1(:,3)==i,:);
        Y(i,1) = mean(M5_Temp(:,1));%mean of rows
        Y(i,2) = mean(M5_Temp(:,2));%mean of columns

       % Compute distances between centroid and cluster pixels
    distance = sqrt((M5_Temp(:,1) - Y(i,1)).^2 + (M5_Temp(:,2) - Y(i,2)).^2); 
    % Find the maximum distance
    maxDistances(i) = max(distance);
    
    end
    
    D_maxD=[Y maxDistances];

    % Number of rows in X and Y
m = size(X, 1);
n = size(D_maxD, 1);

% Initialize a matrix to store distances
Distances = zeros(m, n);

% Calculate distances
for i = 1:m
    for j = 1:n
        % Calculate distance using Euclidean distance formula
        d = norm(X(i, 1:2) - D_maxD(j, 1:2));
        
        % Divide by maxDistances
        Distances(i, j) = d / D_maxD(j,3);
    end
end

% Find the minimum value for each row
min_distances = min(Distances, [], 2);


%Importing the random vesical value range
    %XtempL = readmatrix(XFile,'Range','F:F');
    %XLength = XtempL(1,1);
    %XtempH = readmatrix(XFile,'Range','E:E');
    %XHeight = XtempH(1,1);
    randXx = (nrows).*rand(numVesicles,1);
    randYx = (ncols).*rand(numVesicles,1);
    XRand = [randXx, randYx];

    % Initialize a matrix to store distances
rand_distances = zeros(m, n);

% Calculate distances
for u = 1:m
    for v = 1:n
        % Calculate distance using Euclidean distance formula
        d = norm(XRand(u, 1:2) - D_maxD(v, 1:2));
        
        % Divide by square root of the third element in Y(v,:)
        rand_distances(u, v) = d / D_maxD(v,3);
    end
end

% Find the minimum value for each row
DRand = min(rand_distances, [], 2);


    [Knn,D] = knnsearch(Y,X);  %Nearest neighbour distance (Possible Algorithm improvement)


   
        
    %Output all significant values
    Threshold;
    X;
    Y;
    D;
    DRand;
    min_distances;
    write(X,XRand,Y,D,DRand,min_distances,file,columnNum);
    columnNum = columnNum + numVesicles;
end

% %randomize
% function [XRand, DRand] = randomize(XFile,numVesicles,Y)
%     %Importing the random vesical value range
%     XtempL = readmatrix(XFile,'Range','F:F');
%     XLength = XtempL(2,1);
%     XtempH = readmatrix(XFile,'Range','E:E');
%     XHeight = XtempH(2,1);
%     randXx = (XHeight).*rand(numVesicles,1);
%     randYx = (XLength).*rand(numVesicles,1);
%     XRand = [randXx, randYx];
% 
%     [Knn,DRand] = knnsearch(Y,XRand); 
% end
%%%%%

%write files
function write(X,XRand,Y,D,DRand,min_distances,file,columnNum)
    filename = 'MaxD_norm.xlsx';
    title = ["X Values", "", "Random X Values", "", "Y Values", "", "D Values","Randomized min Values", "min_distances"];
    sheetName = ['Data Set',num2str(file)];
    writematrix(title, filename,'Sheet',sheetName,'Range','A1'); 
    writematrix(X, filename,'Sheet',sheetName,'Range',['A',num2str(2)]);
    writematrix(XRand, filename,'Sheet',sheetName,'Range',['C',num2str(2)]);
    writematrix(Y, filename,'Sheet',sheetName,'Range',['E',num2str(2)]);
    writematrix(D, filename,'Sheet',sheetName,'Range',['G',num2str(2)]);
    writematrix(DRand, filename,'Sheet',sheetName,'Range',['H',num2str(2)]);
    writematrix(min_distances, filename,'Sheet',sheetName,'Range',['I',num2str(2)]);

    %['A',num2str(columnNum)]
    writematrix(min_distances, filename, 'Sheet', "min and Random min Values", 'Range',['A',num2str(columnNum)]);
    writematrix(DRand, filename, 'Sheet', "min and Random min Values", 'Range',['B',num2str(columnNum)]);
end
%%%%%
