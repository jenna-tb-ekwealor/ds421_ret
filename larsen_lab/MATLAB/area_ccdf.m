%% AREA CCDF
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73rep.mat');
state2 = kron(state,ones(1,2)); %to duplicate all columns since cells are twice as large in the y direction as
L = logical(state2);%to convert state2 array to logical array for binary image tools
[labeled,numRidges] = bwlabel(L,4);%to label each ridge (ones) %The parameter 4, passed to the bwlabel function, means that pixels must touch along an edge to be considered connected. For more information about the connectivity of objects, see Pixel Connectivity.

ridgedata = regionprops(labeled,'area', 'centroid', 'perimeter')%to calculate their areas = structure type
perim = bwperim(labeled);%outlines perimeter of objects

Area = [ridgedata.Area].';
Perimeter = [ridgedata.Perimeter].';

Perimeter5 = Perimeter*5;
Area25 = Area*25;

%% Calculate CCDF values
[ycdf,xcdf] = cdfcalc(Area);%calculate ccdf values for ridge area
xccdf = xcdf;
yccdf = 1-ycdf(1:end-1);

%% Plot CCDF of area, uncorrected pixel size
figure, loglog(xccdf, yccdf, 'ko', 'MarkerFaceColor', 'k')%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (25 m^2)','FontSize',14,'FontWeight','Bold');%x axis title
ylabel('Probability Density (P(x))','FontSize',14,'FontWeight','Bold');%y axis title
%axis([(1e2),(1e6),(1e-3),(1e0)])%axes matching Casey et al.
title('Area CCDF, Log-Log Scale','FontSize',18,'FontWeight','Bold');

%% Fit lognormal and generalized pareto models to data
hold on

pd_l = fitdist(xccdf, 'Lognormal');%(or Area; same values but xccdf in order and no duplicates)
pd_g = fitdist(xccdf,'GeneralizedPareto');

y_l = cdf(pd_l, xccdf);%calculate cdf of lognormal fit
y_g = cdf(pd_g, xccdf);%calculate cdf of gp fit

%% Plot model fittings on same figure
plot(xccdf,1-y_g)%1-cdf for ccdf
plot(xccdf,1-y_l)
legend('Area CCDF','Lognormal', 'Generalized Pareto');
hold off

%% Plot CCDF of area, corrected pixel size (1 pixel is 5 m x 5 m)
[ycdf2,xcdf2] = cdfcalc(Area25);%calculate ccdf values for ridge area with corrected pixel siezs
xccdf2 = xcdf2;
yccdf2 = 1-ycdf2(1:end-1);

figure, loglog(xccdf2, yccdf2, 'ko', 'MarkerFaceColor', 'k','MarkerSize',4)%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (m^2)','FontSize',14,'FontWeight','Bold');%x axis title
ylabel('Probability Density (P(x))','FontSize',14,'FontWeight','Bold');%y axis title
%axis([(1e2),(1e6),(1e-2),(1e0)])%axes matching Casey et al.
axis([(1e2),(1e6),(1e-2),(1e0)])
title('Area CCDF, Log-Log Scale, Corrected Pixel Size, rep','FontSize',18,'FontWeight','Bold');

%fit lognormal and generalized pareto to data
hold on

pd_l2 = fitdist(xccdf2, 'Lognormal');%(or Area25; same values but xccdf in order and no duplicates)
pd_g2 = fitdist(xccdf2,'GeneralizedPareto');

y_l2 = cdf(pd_l2, xccdf2);%calculate cdf of lognormal fit
y_g2 = cdf(pd_g2, xccdf2);%calculate cdf of gp fit

plot(xccdf2,1-y_g2,'LineWidth',2)%1-cdf for ccdf
plot(xccdf2,1-y_l2,'LineWidth',2)
legend('Patch','Lognormal', 'Generalized Pareto');
hold off

%:)