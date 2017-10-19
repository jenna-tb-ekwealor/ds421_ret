%% PATCH ELONGATION VS. AREA
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73rep.mat');
state2 = kron(state,ones(1,2)); % duplicate all columns since cells are twice as large in the y direction as
L = logical(state2); % convert state2 array to logical array for binary image tools
[labeled,numRidges] = bwlabel(L,4); % label each ridge (ones) %The parameter 4, passed to the bwlabel function, means that pixels must touch along an edge to be considered connected. For more information about the connectivity of objects, see Pixel Connectivity.

ridgedata = regionprops(labeled,'basic'); % calculate their areas = structure type
ridge_areas = ridgedata.Area;
cc = bwconncomp(L, 4); 

%% one column to test
% first_col = L(:,1); %first column, all rows
% [labeled_col,numRidges_col] = bwlabel(first_col,4);% to label each ridge (1s)
% ridgedata_col = regionprops(labeled_col,'area', 'centroid', 'perimeter'); %calculate their areas = structure type
% sum_length_m_col = 5*(sum([ridgedata_col.Area])); % multiply by 5 m for l
% disp(sum_length_m_col)
% 

%% Loop over each column and row WHEN only grabbing one ridge to get length and width data PER RIDGE
%columns
[dx, dy] = size(L);
E_Patch = [];
Area_ridge = [];
for i = 1:numRidges
total_sum_width_m_col = 0;
total_ridges_cols = 0;
mean_width_m_col = 0;
L_ridge = false(size(L)); % make the image the same dimensions as L but with only the selected ridge showing
L_ridge(cc.PixelIdxList{i}) = true; % select ridge 
%figure, imshow(L_ridge); % uncomment out to see which ridge is being analyzed!
    for n = 1:dy % dy is 372, could also do 372 but keeping it in case different dimensions on another plot
        col = L_ridge(:,n); % column n
        [labeled_col,numRidges_col] = bwlabel(col,4);  % label each "ridge" (slice) in the column
        ridgedata_col = regionprops(labeled_col,'area', 'centroid', 'perimeter'); % calculate the ridge slice areas = structure type
        sum_width_m_col = 5*(sum([ridgedata_col.Area])); % convert to meters
        total_sum_width_m_col = total_sum_width_m_col + sum_width_m_col; % running sum of ridge slice widths
        total_ridges_cols = total_ridges_cols + numRidges_col;
        mean_width_m_cols = total_sum_width_m_col/total_ridges_cols;
    end
%% repeat for rows
total_sum_length_m_row = 0;
total_ridges_rows = 0;
mean_length_m_row = 0;
    for n = 1:dx % dx is 254, could also do 254 but keeping it in case different dimensions on another plot
        row = L(n,:); % row n
        [labeled_row,numRidges_row] = bwlabel(row,4);  % label each "ridge" (slice) in the row
        ridgedata_row = regionprops(labeled_row,'area', 'centroid', 'perimeter'); %calculate the ridge slice areas in pixels = structure type
        sum_length_m_row = 5*(sum([ridgedata_row.Area])); %convert to meters
        total_sum_length_m_row = total_sum_length_m_row + sum_length_m_row;
        total_ridges_rows = total_ridges_rows + numRidges_row;
        mean_length_m_rows = total_sum_length_m_row/total_ridges_rows;
    end
E = mean_length_m_rows/mean_width_m_cols; % mean length to width ratio for ith patch
E_Patch = [E_Patch, E]; % Add patch E to matrix of patch elongation ratios

%% Find area of each ridge
area_ridge = ridgedata(i).Area; % find the area of the ith ridge, in square m
area_ridge25 = (area_ridge)*25;
Area_ridge = [Area_ridge, area_ridge25]; % add area of the ith ridge to the area matrix for plotting

end

%% Plot ridge elongation vs log of ridge area
figure, plot(Area_ridge, log(E_Patch),'ko', 'MarkerFaceColor', 'k','MarkerSize',4)
%showplottool('on','propertyeditor')
xlabel('Patch size (m^2)','FontSize',14,'FontWeight','Bold');%x axis title
ylabel('Log Patch Elongation Ratio','FontSize',14,'FontWeight','Bold');%y axis title
%axis([(1e2),(1e6),(-0.5),(4)])%axes matching Casey et al.
title('Patch Elongation, Corrected Pixel Size - rep','FontSize',18,'FontWeight','Bold');
