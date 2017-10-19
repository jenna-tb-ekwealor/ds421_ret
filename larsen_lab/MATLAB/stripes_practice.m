I = load('/Users/jennabaughman/Documents/berkeley/ds421/larsen_lab/test_array.txt');
%I = logical(test);%convert to logical array
figure, imshow(I)%to display the logical array as a figure
mean2(I);%mean of matrix
%F = fftshift(fft2(I));
F = fftshift(fft2(I-mean2(I)));%fast fourier transform of mean-corrected matrix
%figure, imshow(F)
af1=log(1+abs(F));
fm=max(af1(:));
figure, imshow(im2uint8(af1/fm));
%figure, imshow(mat2gray(log(1+abs(af1))));




%The parameter 'n' in fft(x,n) is the length of the FFT. If you leave the parameter out, Matlab performs an fft of
%length(x). To pad the length out to a power of 2, you have to provide n
%on your own.  So x will be
%padded with (512-497) zeros. FFT's gain their efficiency by using lengths
%that are powers of two. Speed of algorithm may not be an issue for me. 

N = length(I)
plot(abs(fft2(I)))


F2 = fft2(I-mean2(I));
P = log(abs(F2)); %power spectrum "The resulting 2-D periodogram (i.e. spectral density) is a grid
%representing the magnitude of cosine and sine waves of possible wavenumbers (i.e. spatial
%frequencies), and orientations to the spectrum"
figure, imshow(P)

%PP = log(abs(FF));%not mean-corrected
%max(max(PP))%not mean-corrected
%10.7674 %not mean-corrected
%min(min(PP))%not mean-corrected
%-37.0546 %not mean-corrected

max(max(P))%use max and min to rescale values between 0 and 1
%10.2619
min(min(P))
%-Inf after subtracting mean
figure, imshow(P);%frequency spectrum once scaled to 0-1
plot(P)
hist(P)

%%%%%%%%%%%%%%%alternate
img_fft = fft2(I);
img_spec = abs(img_fft);
img_shift = fftshift(img_spec);
idx = find(img_shift == max(img_shift(:)));
img_shift(idx) = NaN;
[u,v] = find(img_shift == max(img_shift(:)),1);%remove 1 to show both results, change to 2 for other result
[rows, cols] = size(img_fft);

u = u - floor(rows/2);
v = v - floor(cols/2);
r = sqrt(u^2 + v^2)
%36.0139
theta = atan2(v, u)
%-1.5430
























%OR?

FA = abs(F); % Get the magnitude
figure, imshow(FA)
plot(FA)
hist(FA)
FAL = log(FA+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
plot(FAL)
figure, imshow(FAL,[]);
hist(FAL)

F2 = log(FA);%visualize discrete fourier transform??
figure, imshow(F2)
hist(F2)
plot(F2)





