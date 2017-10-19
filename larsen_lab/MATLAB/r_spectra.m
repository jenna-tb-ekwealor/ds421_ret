%% Omnidirectional r spectrum
load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/laststepRSL73r2.mat'); % load file
state2 = kron(state,ones(1,2)); % duplicate all columns since cells are twice as large in the y direction as
D = state2;
DD = D - repmat((mean(D)),size(D,1),1); % subtract mean of matrix from every cell

img = DD;
res =25;

%% Process image size information
[N, M] = size(DD);

%% Compute power spectrum
imgf = fftshift(fft2(DD));
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
[X, Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta, rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
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
xMin = 0;                                                                   % No negative image dimension
xMax = ceil(log10(halfDim));
xRange = (xMin:xMax);
yMin = floor(log10(min(Pf)));
yMax = ceil(log10(max(Pf)));
yRange = (yMin:yMax);

% Create plot axis labels
xCell = cell(1:length(xRange));
for i = 1:length(xRange)
    xRangeS = num2str(10^(xRange(i))*res);
    xCell(i) = cellstr(xRangeS);
end
%% Generate plot
figure
loglog(f1,Pf,'ko', 'MarkerFaceColor', 'k','MarkerSize',4)
set(gcf,'color','white')
axis([(1e0),(1e4),(1e-8),(1e-2)])
xlabel('Wavenumber (1/km)','FontSize',14,'FontWeight','Bold');
ylabel('Power','FontSize',14,'FontWeight','Bold');
title('Radially averaged power spectrum, r2','FontSize',18,'FontWeight','Bold')
hold on

%% fit radial r spectra to power law and check how good a fit it is

% f1sort = sort(f1);
% %f1sort = f1sort(1:find(f1sort>=cutoffs(n), 1, 'first')-1);
% [alpha, xmin, L] = plfit(f1sort, 'xmin', min(f1sort));
% [p, gof] = plpva(f1sort, xmin, 'xmin', min(f1sort));
% plplot(f1sort, xmin, alpha);
% set(gcf,'color','white');
