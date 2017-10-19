load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73rep.mat');
state2 = kron(state,ones(1,2)); % duplicate all columns since cells are twice as large in the y direction
L = logical(state2); % convert state2 array to logical array for binary image tools

[labeled,numRidges] = bwlabel(L,4); % label each ridge (1s)  % The parameter 4, passed to the bwlabel function, means that pixels must touch along an edge to be considered connected.
ridgedata = regionprops(labeled,'area', 'centroid', 'perimeter')%to calculate their areas = structure type
perim = bwperim(labeled); % outlines perimeter of objects
Area = [ridgedata.Area].';
Perimeter = [ridgedata.Perimeter].';

figure, imshowpair(imrotate(L,90),imrotate(perim,90),'montage')
showplottool('on','propertyeditor')
title('RASCAL Binary Image and Ridges Outlined, rep','FontSize',18,'FontWeight','Bold');

