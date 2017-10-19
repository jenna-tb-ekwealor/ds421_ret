
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73rep.mat');
state2 = kron(state,ones(1,2)); %to duplicate all columns since cells are twice as large in the y direction as
xx = kron(x,ones(1,2));%y are 1/2x so need to duplicate all columns
L = logical(state2);%to convert state2 array to logical array for binary image tools
[labeled,numRidges] = bwlabel(L,4);%to label each ridge (1s)
%The parameter 4, passed to the bwlabel function, means that pixels must touch along an edge to be considered connected. For more information about the connectivity of objects, see Pixel Connectivity.

ridgedata = regionprops(labeled,'area', 'centroid', 'perimeter')%to calculate their areas = structure type
perim = bwperim(labeled);%outlines perimeter of objects

Area = [ridgedata.Area].';
Perimeter = [ridgedata.Perimeter].';
figure(4), loglog(Area, Perimeter, 'ko', 'MarkerFaceColor', 'k')
showplottool('on','propertyeditor') 
axis([(1e0),(1e6),(1e0),(1e6)])
xlabel('Patch size (25 m^2)');
ylabel('Perimeter (5 m)');
title('Perimeter to Area Scaling; Uncorrected Pixel Size, rep','FontSize',18,'FontWeight','Bold');

Perimeter5 = Perimeter*5;
Area25 = Area*25;

figure(5), loglog(Area25, Perimeter5, 'ko', 'MarkerFaceColor', 'k','MarkerSize',4)
showplottool('on','propertyeditor') 
title('Perimeter to Area Scaling; Corrected Pixel Size, rep','FontSize',18,'FontWeight','Bold');
xlabel('Patch size (m^2)','FontSize',14,'FontWeight','Bold');
ylabel('Perimeter (m)','FontSize',14,'FontWeight','Bold');
axis([(min(Area25)),(1e6),(min(Perimeter5)),(1e5)])
hold on

lArea25 = log(Area25);
lPerimeter5 = log(Perimeter5);

plin = polyfit(lArea25, lPerimeter5, 1); %find linear fit coefs
x_fitlin = linspace(min(lArea25),max(lArea25),2); %how to choose?
y_fitlin = plin(1)*x_fitlin + plin(2);
x_fitlinlog = exp(x_fitlin);
y_fitlinlog = exp(y_fitlin);
hold on
plot(x_fitlinlog,y_fitlinlog,'-b', 'LineWidth',2) %plot linear fit

pquad = polyfit(lArea25, lPerimeter5, 2); %find quad fit coefs
x_fitquad = linspace(min(lArea25),max(lArea25),100);
%x_fitquad = linspace(0,12.4,100);
y_fitquad = pquad(1)*(x_fitquad.^2) + pquad(2)*x_fitquad + pquad(3);
x_fitquadlog = exp(x_fitquad);
y_fitquadlog = exp(y_fitquad);
plot(x_fitquadlog, y_fitquadlog,'-r','LineWidth',2) %plot quad fit
legend('Patch','Linear', 'Quadratic');
hold off



%:)