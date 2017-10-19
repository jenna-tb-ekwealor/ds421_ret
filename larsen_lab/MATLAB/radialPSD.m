function radialPSD(img,res)
% function raPsd2d(img,res)
%
% Computes and plots radially averaged power spectral density (power
% spectrum) of image IMG with spatial resolution RES.
%
% (C) E. Ruzanski, RCG, 2009
% % % % Using modified version! % % % % 

%% Process image size information
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

%% Generate plot
figure
loglog(f1',Pf,'ko', 'MarkerFaceColor', 'k') %need spectral density vs wavenumber
xlabel('Wavenumber (1/km)');%x axis title
ylabel('Power');%y axis title
axis([(1e0),(1e3),(1e-6),(1e-2)])%axes to matching Casey et al.?
title('Radially averaged power spectrum');
hold on

%% fit power law
[alpha, xmin, L] = plfit(f1', 'xmin', min(f1'));
[p, gof] = plpva(f1', xmin, 'xmin', min(f1'));
plplot(f1', xmin, alpha)
%    The output 'alpha' is the maximum likelihood estimate of the scaling
%    exponent, 'xmin' is the estimate of the lower bound of the power-law
%    behavior, and L is the log-likelihood of the data x>=xmin under the
%    fitted power law.
%    
% %How good a fit is the power law distribution?
% %     
% %     [alpha(n), xmin, L] = plfit(rt, 'xmin', min(rt));
% %     [p(n), gof(n)] = plpva(rt, xmin, 'xmin', min(rt));
% %     plplot(rt, xmin, alpha(n))
hold off
end