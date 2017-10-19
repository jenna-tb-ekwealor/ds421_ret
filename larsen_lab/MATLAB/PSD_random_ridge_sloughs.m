%PSD_random_ridge_sloughs.m
%This routine randomly generates coherent ridges and sloughs (1's and 0's)
%that are spaced at a mean distance of 150 m, though with varying degrees
%of noise. Pixels are in meters, and the 1-D domain is 6 km. The code then
%generates a power spectrum, which is the average power spectrum of 1000
%realizations of these 1D slices through the landscape. 

noise = [25 50 75 100 125 150]; %In meters
for n = 1:6
Psdx = [];
for jj = 1:1000
    Test = [];
    N = 6000;
    while length(Test)<N
        len0 = max(0,round(150+noise(n)*randn(1)));
        len1 = max(0,round(150+noise(n)*randn(1)));
        Test = [Test; zeros(len0,1); ones(len1,1)];
    end
    Test = Test(1:N);
    xdft = fft(Test);
    xdft = xdft(1:N/2+1);
    psdx = abs(xdft).^2; %Normalized would be x1/N
    psdx(2:end-1) = 2*psdx(2:end-1); %Final power spectrum up to Nyquist freq
    Psdx = [Psdx, psdx]; %Add this to the matrix of power spectra
end
    freq = (0:1/N:1/2)*1000; %cycles/km
    subplot(2,3,n)
    loglog(freq, mean(Psdx,2), 'k.')
    xlabel('Wave number, cycles/km')
    ylabel('PSD')
    title(sprintf('%s%u', 'Noise level ', noise(n)))
    set(gca, 'XLim', [1e-1 1e3])
end
