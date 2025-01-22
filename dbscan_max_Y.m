%Define number of batch processing
num = 5;

%%%%%

%Read Files
for file=1:num
    XFile = file + "_PSD" + ".xlsx";
    YFile = file + "_bg" + ".txt";
    processing(XFile, YFile,file);
end
"Done Processing"
%%%%%

%process files
function processing(XFile,YFile,file)
%Importing the bg value
bg_value = importdata(YFile);

%Importing the gray value data (X) (Possible Coding improvements)
    inputX = readmatrix (XFile);
    [nrows, ncols]=size(inputX); %stores the number of rows and columns in input
    values=[]; %output matrix
    for r=1:nrows
        for c=1:ncols
            values=[values;r c inputX(r,c)];% keep appending [r,c,input] to new columns
        end
    end
    

%Import backgroud value
TF = values(:,3)<bg_value;
M3 = values(TF, :);
M3c = M3(:,3);

%Finding top 20% value threshold
[m, n] = size(M3c);
TN = m*n;
TTP = TN*.2;
TTP = ceil(TTP);
[val,ind] = sort(M3c);
Thres = val(TTP);

%Clustring using DBSCAN
    M_Thres = values(values(:,3)<Thres,1:2);
    idx = dbscan(M_Thres,7,32); 

%Finding the centroids of the clusters
    M_Cluster = [idx,M_Thres];
    numOfIdx = size(unique(idx))-1;%the number of valid indexes
    Y(numOfIdx(1),2) = 0;%matrix holding mean values
    for i= 1:numOfIdx(1)
        M5_Temp = M_Cluster(M_Cluster(:,1)==i,:);
        Y(i,1) = mean(M5_Temp(:,2));%mean of rows
        Y(i,2) = mean(M5_Temp(:,3));%mean of columns
    end

%Delete points which belong to cluster -1 and plot
    N = [M_Thres idx];
    idx2 = idx(idx~= -1);
    numGroups = length(unique(idx2));
    clr = hsv(numGroups);
    N1=N(N(:,3) ~= -1,:);
    h = figure();
    set(h,'Visible', 'off');
    gscatter(N(N(:,3) ~= -1,1),N(N(:,3) ~= -1,2),idx2,clr);

filename = strcat("DataSet",num2str(file));
    axis image;
    xlabel("X length");
    ylabel("Y length");
    title(filename);
    filename = strcat(filename,".tif");
    legend('off');
    exportgraphics(h,filename,'Resolution',400)

%Count the number of pixels in each cluster 
t = unique(idx2); % the unique values in the idx2     
 for i = 1:length(t)
   counts(i,1) = sum(idx2==t(i)); % number of times each unique value is repeated
 end
 % t(i) is repated count(i) times

%area of each cluster
for i = 1:length(t)
tf = N1(:,3)==i; %i is the No. of clusters.
p=N1(tf,:);
x=p(:,1);
y=p(:,2);
dt = delaunayTriangulation(x,y);
k = convexHull(dt);
%plot(x,y, '.', 'markersize',10);
%hold on;
%plot(x(k), y(k), 'r');
% Perimeter = sqrt(diff(x(k))*diff(x(k))'+ diff(y(k))*diff(y(k))'); % Perimeter
area(i,1)=abs(trapz(x(k),y(k)));  
end

% Nearest distance between clusters
for e = 1:length(t)
    for f = 1:length(t)
D_temp(f,1)=sqrt((Y(e,1)-Y(f,1)).^2+(Y(e,2)-Y(f,2)).^2);
    end
Z=sort(D_temp); D(e,1)=Z(2);
end

%Finding maximum Y coordinate distance of each cluster
for j = 1:length(t)
   W(j,1) = max(N(N(:,3)==j,2))-min(N(N(:,3)==j,2));
end

%Output all significant values
    Thres;
    area;
    D;
    W;
    write(area,D,W,file);

end

%write files
function write(area,D,W,file)
filename = 'DBSCAN_width_distance.xlsx';
title = ["Area_Values","", "D Values", "", "Width_max_Y_distance"];
sheetName = num2str(file);
writematrix(title, filename,'Sheet',sheetName,'Range','A1');
writematrix(area, filename,'Sheet',sheetName,'Range',['A',num2str(2)]);
writematrix(D, filename,'Sheet',sheetName,'Range',['C',num2str(2)]);
writematrix(W, filename,'Sheet',sheetName,'Range',['E',num2str(2)]);

end
%%%%%
