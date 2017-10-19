function raPsd2d(img,res)
% function raPsd2d(img,res)
%
% Computes and plots radially averaged power spectral density (power
% spectrum) of image IMG with spatial resolution RES.
%
% (C) E. Ruzanski, RCG, 2009

%% Load image and process size information
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat');
state2 = kron(state,ones(1,2)); %to duplicate all columns since cells are twice as large in the y direction as
xx = kron(x,ones(1,2));%y are 1/2x so need to duplicate all columns
L = logical(state2);%to convert state2 array to logical array for binary image tools
D = state2;
DD = D - repmat((mean(D)),size(D,1),1);%subtract mean of matrix from every cell
img = DD;
res =5; %may need to correct to 1 to match casey et al!

[N M] = size(img);

%% Compute power spectrum
imgf = fftshift(fft2(img));
imgfp = (abs(imgf)/(N*M)).^2;                                               % Normalize

%% Adjust PSD size
dimDiff = abs(N-M);
dimMax = max(N,M);
% Make square
if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
end

halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)

%% Compute radially average power spectrum
[X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho);
i = cell(floor(dimMax/2) + 1, 1);
for r = 0:floor(dimMax/2)
    i{r + 1} = find(rho == r);
end
Pf = zeros(1, floor(dimMax/2)+1);
for r = 0:floor(dimMax/2)
    Pf(1, r + 1) = nanmean( imgfp( i{r+1} ) );
end

%% Setup plot
fontSize = 14;
maxX = 10^(ceil(log10(halfDim)));
f1 = linspace(1,maxX,length(Pf));                                           % Set abscissa

% Find axes boundaries
%xMin = 0;                                                                   % No negative image dimension
%xMax = ceil(log10(halfDim));
%xRange = (xMin:xMax);
%yMin = floor(log10(min(Pf)));
%yMax = ceil(log10(max(Pf)));
%yRange = (yMin:yMax);

% Create plot axis labels
%xCell = cell(1:length(xRange));
%for i = 1:length(xRange)
   % xRangeS = num2str(10^(xRange(i))*res);
  %  xCell(i) = cellstr(xRangeS);
%end
%yCell = cell(1:length(yRange));
%for i = 1:length(yRange)
 %   yRangeS = num2str(yRange(i));
   % yCell(i) = strcat('10e',cellstr(yRangeS));
%end

%% Generate plot
figure
loglog(f1',Pf,'ko', 'MarkerFaceColor', 'k') %need spectral density vs wavenumber
xlabel('Wavenumber (1/km)');%x axis title
ylabel('Power');%y axis title
axis([(1e0),(1e3),(1e-6),(1e-2)])%axes matching Casey et al.
title('Radially averaged power spectrum');
hold on

%% fit power law


% %How good a fit is the power law distribution?
% %     rt = sort(Res_times{n});
% %     rt = rt(1:find(rt>=cutoffs(n), 1, 'first')-1);
% %     [alpha(n), xmin, L] = plfit(rt, 'xmin', min(rt));
% %     [p(n), gof(n)] = plpva(rt, xmin, 'xmin', min(rt));
% %     plplot(rt, xmin, alpha(n))

legend('Spectral Power','Generalized Pareto');
hold off
%% other stuff for changing font size etc. remove eventually
%set(gcf,'color','white')
%set(gca,'FontSize',fontSize,'FontWeight','bold','YMinorTick','off',...
    %'XTickLabel',xCell,'XGrid','on','YAxisLocation','right','XDir','reverse');
%xlabel('Wavelength (km)','FontSize',fontSize,'FontWeight','Bold');
%ylabel('Power','FontSize',fontSize,'FontWeight','Bold');
%title('Radially averaged power spectrum','FontSize',fontSize,'FontWeight','Bold')
end