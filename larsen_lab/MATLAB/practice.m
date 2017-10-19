load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat');
state2 = kron(state,ones(1,2)); %to duplicate all columns since cells are twice as large in the y direction as
xx = kron(x,ones(1,2));%y are 1/2x so need to duplicate all columns
L = logical(state2);%to convert state2 array to logical array for binary image tools
D = state2;
DD = D - repmat((mean(D)),size(D,1),1);%subtract mean of matrix from every cell
radialPSD(DD,5); %this function does the power spectra density and fits power law; may have to change axes manually
