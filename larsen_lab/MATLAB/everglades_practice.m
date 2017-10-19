% %to open file
% load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat');
% state2 = kron(state,ones(1,2)); %to duplicate all columns since cells are twice as large in the y direction as
% %they are in the x direction
% %state2=state(ceil((1:2*size(state,1))/2), :); %another way
% xx = kron(x,ones(1,2));%y are 1/2x so need to duplicate all columns
% %figure, pcolor(y,x,state'), shading('flat'), colormap('gray'), axis('equal')%make figure with original array
% %figure, pcolor(y,xx,state2'), shading('flat'), colormap('gray'), axis('equal')%make figure with new array
% L = logical(state2);%to convert state2 array to logical array for binary image tools
% %figure, imshow(imrotate(L,90))%to display the logical array as a figure
% [labeled,numRidges] = bwlabel(L,4);%to label each ridge (1s)
% %The parameter 4, passed to the bwlabel function, means that pixels must touch along an edge to be considered connected. For more information about the connectivity of objects, see Pixel Connectivity.
% %labeled = double matrix, bwconncomp can make one too of class uint8 which
% %is more memory efficient
% numRidges
% %44
% 
% 
% 
% 
% 
% %% try to plot elongation (mean length of a ridge over mean width of a
% % ridge) with area of ridge
% 
% % % 
% % % %%%%patch elongation(log-transformed length to width ratio) vs log scale of patch area
% % 
% % [labeled_all,numRidges] = bwlabel(L,4);  % label each ridge in the site
% % ridgedata_all = regionprops(labeled_all,'area', 'centroid', 'perimeter'); %calculate the ridge slice areas in pixels = structure type
% % 
% % perim = bwperim(labeled);%outlines perimeter of objects
% % 
% % Area = [ridgedata_all.Area].';
% % Perimeter = [ridgedata_all.Perimeter].';
% % 
% % Perimeter5 = Perimeter*5;
% % Area25 = Area*25;
% 
% %figure, scatter(E,Area25, 'ko', 'MarkerFaceColor', 'k')
% 
% 
% 
% %%










load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat');
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
title('Perimeter to Area Scaling; Uncorrected Pixel Size','FontSize',18,'FontWeight','Bold');

Perimeter5 = Perimeter*5;
Area25 = Area*25;

figure(5), loglog(Area25, Perimeter5, 'ko', 'MarkerFaceColor', 'k','MarkerSize',4)
showplottool('on','propertyeditor') 
title('Perimeter to Area Scaling; Corrected Pixel Size, r2','FontSize',18,'FontWeight','Bold');
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




% % % 
% % % %%%%patch elongation(log-transformed length to width ratio) vs log scale of
% % % %%%%patch area
% % % 
% % % p_elong = log(length)/log(width);
% % % Area25
% % % figure, scatter(p_elong,Area25, 'ko', 'MarkerFaceColor', 'k')




figure(7), imshowpair(imrotate(L,90),imrotate(perim,90),'montage')
showplottool('on','propertyeditor')
title('RASCAL Binary Image and Ridges Outlined, r2','FontSize',18,'FontWeight','Bold');

[ycdf,xcdf] = cdfcalc(Area);%calculate ccdf values for ridge area
xccdf = xcdf;
yccdf = 1-ycdf(1:end-1);

figure(8), loglog(xccdf, yccdf, 'ko', 'MarkerFaceColor', 'k')%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (25 m^2)','FontSize',14,'FontWeight','Bold');%x axis title
ylabel('Probability Density (P(x))','FontSize',14,'FontWeight','Bold');%y axis title
%axis([(1e2),(1e6),(1e-3),(1e0)])%axes matching Casey et al.
title('Area CCDF, Log-Log Scale, r2','FontSize',18,'FontWeight','Bold');

%fit lognormal and generalized pareto to data
hold on

pd_l = fitdist(xccdf, 'Lognormal');%(or Area; same values but xccdf in order and no duplicates)
pd_g = fitdist(xccdf,'GeneralizedPareto');

y_l = cdf(pd_l, xccdf);%calculate cdf of lognormal fit
y_g = cdf(pd_g, xccdf);%calculate cdf of gp fit

plot(xccdf,1-y_g)%1-cdf for ccdf
plot(xccdf,1-y_l)
legend('Area CCDF','Lognormal', 'Generalized Pareto');
hold off

%r2eat with corrected pixel size
[ycdf2,xcdf2] = cdfcalc(Area25);%calculate ccdf values for ridge area with corrected pixel siezs
xccdf2 = xcdf2;
yccdf2 = 1-ycdf2(1:end-1);

figure(9), loglog(xccdf2, yccdf2, 'ko', 'MarkerFaceColor', 'k','MarkerSize',4)%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (m^2)','FontSize',14,'FontWeight','Bold');%x axis title
ylabel('Probability Density (P(x))','FontSize',14,'FontWeight','Bold');%y axis title
%axis([(1e2),(1e6),(1e-2),(1e0)])%axes matching Casey et al.
axis([(1e2),(1e6),(1e-2),(1e0)])
title('Area CCDF, Log-Log Scale, Corrected Pixel Size, r2','FontSize',18,'FontWeight','Bold');

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








figure(7), imshowpair(imrotate(L,90),imrotate(perim,90),'montage')
showplottool('on','propertyeditor')
title('RASCAL Binary Image and Ridges Outlined','FontSize',18,'FontWeight','Bold');





F = fft2(DD); %shift for viewing %discrete 2D fourier transform of binary image where D is the size-corrected double (state2)
figure, imshow(F) %unshifted fourier spectrum
af1=log(1+abs(F)); %log function
fm=max(af1(:)); %find max
M = im2uint8(af1/fm); %show fourier correcting for max??
figure, imshow(M)
title('Image FFT2 Magnitude, log & normalized by max')
showplottool('on','propertyeditor')

%use magnitude to find r spectra
%casey et al +- 10 degrees from perpendicular to flow
%scatter plot of all columns, direction perpendicular to flow (flow is
%horizontal)


%MAGNITUDE(F) = SQRT( REAL(F)^2+IMAGINARY(F)^2 ) ====> r????



figure(1),imshow(imrotate(L,90))
showplottool('on','propertyeditor') 
title('Size-Corrected Logical Array');
figure(2),imshow(abs(fft2(imrotate(L,90))))
showplottool('on','propertyeditor') 
title('FFT2 Magnitude of Size-Corrected Logical Array');
plot(log(abs(fftshift(fft2(L)))))
imagesc(log(abs(fftshift(fft2(L)))))




showplottool('on','propertyeditor') 
title('Plot of Magnitude of Size-Corrected Logical Array');
figure(3),imshow(imrotate(labeled,90))
showplottool('on','propertyeditor') 
title('Ridge-Labeled Size-Corrected Double Array');
impixelregion %to examine pixels and labels

ridgedata = regionprops(labeled,'area', 'centroid', 'perimeter')%to calculate their areas = structure type
max_area = max([ridgedata.Area])%find highest ridge area
%8148
max_area_m2 = 25*(max([ridgedata.Area])) %find highest ridge area

sum_area_m2 = 25*(sum([ridgedata.Area]))

sum_perim = sum([ridgedata.Perimeter])
% 7.2002e+03

sum_perim_m = 5*(sum([ridgedata.Perimeter]))
% 3.6001e+04


mean_perim = mean([ridgedata.Perimeter])%find mean perimeter of all ridges
%163.6414 %is in pixel lengths, where each pixel is 5x5 m^2. multiply by 5
%for m length

mean_perim_m = 5*(mean([ridgedata.Perimeter]))%find mean perimeter of all ridges
%  818.2068


mean_area = mean([ridgedata.Area])%find mean of all ridge areas
%996.8636
mean_area_m2 = 25*(mean([ridgedata.Area]))%find mean of all ridge areas
%   2.4922e+04

[dx, dy] = size(L)

site_area_m2 = 25*(dx*dy)
%   2362200

ridge_density = sum_area_m2/site_area_m2
%0.4642

edge_density = sum_perim_m/site_area_m2
%    0.0152














perim = bwperim(labeled);%outlines perimeter of objects

Area = [ridgedata.Area].';
Perimeter = [ridgedata.Perimeter].';
figure(4), loglog(Area, Perimeter, 'ko', 'MarkerFaceColor', 'k')
showplottool('on','propertyeditor') 
axis([(1e0),(1e6),(1e0),(1e6)])
xlabel('Patch size (25 m^2)');
ylabel('Perimeter (5 m)');
title('Perimeter to Area Scaling; Uncorrected Pixel Size');


Perimeter5 = Perimeter*5;
Area25 = Area*25;
figure(5), loglog(Area25, Perimeter5, 'ko', 'MarkerFaceColor', 'k')
showplottool('on','propertyeditor') 
title('Perimeter to Area Scaling; Corrected Pixel Size');
xlabel('Patch size (m^2)');
ylabel('Perimeter (m)');
axis([(1e0),(1e6),(1e0),(1e6)])


figure(6), imshow(imrotate(perim,90))%show outlines and rotate
showplottool('on','propertyeditor')
title('Ridges Outlined on Size-Corrected Logical Array');

figure(7), imshowpair(imrotate(L,90),imrotate(perim,90),'montage')
showplottool('on','propertyeditor')
title('RASCAL Binary Image and Ridges Outlined on Size-Corrected Logical Array');


Area = [ridgedata.Area].';
[ycdf,xcdf] = cdfcalc(Area);%calculate ccdf values for ridge area
xccdf = xcdf;
yccdf = 1-ycdf(1:end-1);


figure(8), loglog(xccdf, yccdf, 'ko', 'MarkerFaceColor', 'k')%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (25 m^2)');%x axis title
ylabel('Probability Density (P(x))');%y axis title
axis([(1e2),(1e6),(1e-3),(1e0)])%axes matching Casey et al.
title('Area CCDF, Log-Log Scale');

%fit lognormal and generalized pareto to data
hold on

pd_l = fitdist(xccdf, 'Lognormal');%(or Area; same values but xccdf in order and no duplicates)
pd_g = fitdist(xccdf,'GeneralizedPareto');

y_l = cdf(pd_l, xccdf);%calculate cdf of lognormal fit
y_g = cdf(pd_g, xccdf);%calculate cdf of gp fit

plot(xccdf,1-y_g)%1-cdf for ccdf
plot(xccdf,1-y_l)
legend('Area CCDF','Lognormal', 'Generalized Pareto');
hold off




%r2eat with corrected pixel size
[ycdf2,xcdf2] = cdfcalc(Area25);%calculate ccdf values for ridge area with corrected pixel siezs
xccdf2 = xcdf2;
yccdf2 = 1-ycdf2(1:end-1);

figure(9), loglog(xccdf2, yccdf2, 'ko', 'MarkerFaceColor', 'k')%plot ccdf of area
showplottool('on','propertyeditor')
xlabel('Patch size (m^2)');%x axis title
ylabel('Probability Density (P(x))');%y axis title
axis([(1e2),(1e6),(1e-3),(1e0)])%axes matching Casey et al.
title('Area CCDF, Log-Log Scale, Corrected Pixel Size');

%fit lognormal and generalized pareto to data
hold on

pd_l2 = fitdist(xccdf2, 'Lognormal');%(or Area25; same values but xccdf in order and no duplicates)
pd_g2 = fitdist(xccdf2,'GeneralizedPareto');

y_l2 = cdf(pd_l2, xccdf2);%calculate cdf of lognormal fit
y_g2 = cdf(pd_g2, xccdf2);%calculate cdf of gp fit

plot(xccdf2,1-y_g2)%1-cdf for ccdf
plot(xccdf2,1-y_l2)
legend('Area CCDF','Lognormal', 'Generalized Pareto');
hold off

%:)
































%%r spectra


%subtract mean of the matrix from the matrix before applying the fft

%L is size corrected logical array of binary image
%subtract mean of matrix from matrix
%D = double(L);
D = state2;
DD = D - r2mat((mean(D)),size(D,1),1);%subtract mean of matrix from every cell
%figure, imshow(DD)%see mean-corrected image
F = fft2(DD); %shift for viewing %discrete 2D fourier transform of binary image where D is the size-corrected double (state2)
%inverse = ifft2(F); %testing
%figure, imshow(inverse) %testing
figure, imshow(F) %unshifted fourier spectrum
af1=log(1+abs(F)); %log function
fm=max(af1(:)); %find max
M = im2uint8(af1/fm); %show fourier correcting for max??
figure, imshow(M)
title('Image FFT2 Magnitude, log & normalized by max')
showplottool('on','propertyeditor')

%only need fftshift for viewing fourier transform; may want to keep in
%unshifted to plot
%display magnitude and phase of 2d fft, where F is the fft of the
%normalized double image

%figure, imshow(log(abs(F)),[24 100000]), colormap gray
%The parameters [24 100000] are the limits of the colormap. 
%Values below 24 are displayed as black, and values above 100000 are displayed as white, 
%while values in between are displayed as grayscale values.
 %magnitude, log transformed

%use magnitude to find r spectra
%casey et al +- 10 degrees from perpendicular to flow
%here going to do just perpendicular first
%scatter plot of all columns, direction perpendicular to flow (flow is
%horizontal)

s


%use power spectral density?
n = 254;
m = 372;
a = 5;
PSD = (a^2/n*m)*(abs(log(abs(F))));
figure, imshow(PSD), colormap gray
showplottool('on','propertyeditor')
%where a is pixel width and n & m are number of pixels in each direction.

%tic;Y=fft(PSD);toc %tic toc to calculate time elapsed!
%fft of each column
%fft(X,[],1) operates along the columns of X and returns the Fourier transform of each column.
%fft(X,[],2) operates along the rows of X and returns the Fourier transform of each row.
%fft of array will operate column by column by default
maybe = fft(PSD,[],1);%fft of each column of power spectrum


%M = mean(A,dim) returns the mean along dimension dim. 
%For example, if A is a matrix, then mean(A,2) is a column vector containing the mean of each row.

%average all columns or rows??

a1 = mean(maybe,1);
a2 = mean(maybe,2);
a1abs = (abs(a1));
a2abs = (abs(a2));








figure, imshow(angle(F),[-pi pi]), colormap gray
showplottool('on','propertyeditor')
%I chose [-pi pi] as the phase limits so that the center of the color spectrum of the graph would be at zero phase. 
%The phase limits could also have been from 0 to 2*pi. Either one is fine.
figure, imshow(angle(F)) %phase
showplottool('on','propertyeditor')
figure, imagesc(angle(F)), colormap gray
title('Image FFT2 Phase')
showplottool('on','propertyeditor')



imagesc(abs(F),[24 100000]), colormap gray
title('Image FFT2 Magnitude')
showplottool('on','propertyeditor')


figure, imshow(DD)
showplottool('on','propertyeditor')
figure, imshow(F)
showplottool('on','propertyeditor')
figure, imshow(abs(F))
showplottool('on','propertyeditor')

%visualize discrete fourier transform
F2 = log2(F);
figure, imshow(F2)
showplottool('on','propertyeditor')
figure, imshow(abs(F2))
showplottool('on','propertyeditor')

%create 2D periodogram
%pL = periodogram(D); %includes fourier transform
%figure, plot(pL)










state2norm = state2 - r2mat((mean(state2)),size(state2,1),1);
img_fft = fft2(state2norm);
img_spec = abs(img_fft);
img_spec_ln = log(img_spec);
img_shift = fftshift(img_spec_ln);

%Y = log(X) returns the natural logarithm ln(x) of each element in array X.
%The log function's domain includes negative and complex numbers, which can lead to unexpected results if used unintentionally. For negative and complex numbers z = u + i*w, the complex logarithm log(z) returns
%log(abs(z)) + 1i*angle(z)

%[q,p] = find(img_shift == max(img_shift(:)));

%find the peak (may be the DC value)
%if want to find the peak that's not the DC value then do:
idx = find(img_shift == max(img_shift(:)));
img_shift(idx) = NaN;
%followed by:
%to get both answers: [q,p] = find(img_shift == max(img_shift(:)));
[q,p] = find(img_shift == max(img_shift(:)), 1);

%since shifted::::
%However, this isn't with respect to the centre as you have shifted the image, so once you're done with the above, you need to do:
[rows, cols] = size(img_fft); %rows = 254 cols = 372
qn = q - floor(rows/2);
pn = p - floor(cols/2);


r = sqrt(pn^2 + qn^2) %magnitude
theta = atan2(pn, qn) %phase

%all standard FT frequency images are magnitude not phase

%Here, q would be the horizontal location of your image while p would be the vertical position. 


%M = mean(A,dim) returns the mean along dimension dim. 
%For example, if A is a matrix, then mean(A,2) is a column vector containing the mean of each row.
Mrows = mean(img_shift,2); % vector containing mean of each row %all -Inf %q??
Mcols = mean(img_shift,1); % vector containing mean of each column %p??

qq = Mrows;
pp = Mcols;

r_values = sqrt(pp.^2 + qq.^2);
theta_values = atan2(pp, qq);



ps1=abs(img_fft).^2;% power spectrum using fft
histogram(img_shift)

figure;
pcolor(log(1+abs(img_fft)));
axis ij;

surfl(img_shift)
mesh(img_shift)

%MAGNITUDE(F) = SQRT( REAL(F)^2+IMAGINARY(F)^2 ) ====> r????
%PHASE(F) = ATAN( IMAGINARY(F)/REAL(F) ) ====> theta???

%The MAGNITUDE tells "how much" of a certain frequency component is present and 
%the PHASE tells "where" the frequency component is in the image.




%or for unshifted
%find the peak (may be the DC value)
%if want to find the peak that's not the DC value then do:
idx2 = find(img_spec_ln == max(img_spec_ln(:)));
img_spec_ln(idx2) = NaN;
%followed by:
%to get both answers: [q,p] = find(img_spec_ln == max(img_spec_ln(:)));
[q2,p2] = find(img_spec_ln == max(img_spec_ln(:)), 1);

r2 = sqrt(p2^2 + q2^2) %5.3852
theta2 = atan2(p2, q2) %0.3805


%?? not getting the same results for unshifted












%dummy data
test = load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/test_array.txt');

I = test - r2mat((mean(test)),size(test,1),1);

F=fft2(I);
figure, imshow(F)
%S=fftshift(F);
L=log2(F);
figure, imshow(L)
A=abs(L);
figure, imshow(A)
imagesc(L)
imagesc(A)



























%r2EAT FOR SECOND IMAGE
%to open file
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat')

state2 = kron(state,ones(1,2));
whos 
xx = kron(x,ones(1,2));
figure, pcolor(y,x,state'), shading('flat'), colormap('gray'), axis('equal')
figure, pcolor(y,xx,state2'), shading('flat'), colormap('gray'), axis('equal')
state2logi = logical(state2);
%figure, imshow(state2logi)

figure, imshow(imrotate(state2logi, 90))
[labeled,numRidges] = bwlabel(state2logi,4);

numRidges

figure, imshow(imrotate(labeled,90))
impixelregion

ridgedata = regionprops(labeled,'basic')

maxarea = max([ridgedata.Area])

meanarea = mean([ridgedata.Area])
perim = bwperim(labeled);

Area = [ridgedata.Area].';


%convert axes to log scale
set(gca, 'xscale','log')
set(gca, 'yscale','log')

load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat')

test = load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/test_array.txt')

I = test;
figure, imshow(I)
F = fft2(I);

FC = fftshift(F); % Center FFT

FCA = abs(FC); % Get the magnitude
FCAL = log(FCA+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
FCALG = mat2gray(FCAL); % Use mat2gray to scale the image between 0 and 1
imshow(FCALG,[]);

img_shift = fftshift(FCALG);
plot(img_shift)

test = load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/test_array.txt')






