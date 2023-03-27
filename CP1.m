%Why this homework was challenging for me 
% - The biggest issue was not knowing the difference between fft and ffftn
% (the latter is for arrays of higher dimensions) 
% - Also not knowing if the autograder expected a shifted transform or
% unshifted transform 
% - Understanding what a tensor is and what the Kraken data file actually
% represented also took me a long time 
% - Finding the max values was also challenging (process is to find the
% maximum values in the whole matrix and then to find the corresponding
% indexes in the x,y,z/ ks,ks,ks domains


% Clean workspace
 clear all; close all; clc
 load('Kraken.mat') %load data matrix
 L = 10; % spatial domain
 n = 64; % Fourier modes
 x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x; %create 3D axis arrays with 64 points
 k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; %create frequency array and rescale them to be 2pi periodic 
 ks = fftshift(k); %shift values to order them correctly
 x1 = linspace(-L,L,64);
 %Create 3D grids for both spatial domain and frequency domain 
 [X,Y,Z] = meshgrid(x,y,z);
 [Kx,Ky,Kz] = meshgrid(ks,ks,ks);
 threshold = 0.3; %lower bound for frequencies that are recognized
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% sum all realizations in frequency space
transformsSum = zeros(n,n,n);
for j = 1:49
    Un(:, :, :) = reshape(Kraken(:, j), n, n, n);
    transformsSum= transformsSum  + fftn(Un);
end
%transformsSum = fftshift(transformsSum);

 % after your for loop ends, save the sum as variable A1
 A1 = transformsSum;


 % Average the sum over the 49 realizations (i.e., A1/49) and save as A2
 A2 = A1/49;

 % store largest element in Un (tensor)
 M = max(abs(A2), [], 'all');
 %plot = 2D plane, isosurface = 3D plane, plots points above threshold 
 figure(1)
 close all, isosurface(Kx, Ky, Kz, abs(A2)/M,threshold)
 title('Averaged Frequenices of Seismograph Realizations','FontSize',15)
xlabel('k');
ylabel('k');
zlabel('k');
 axis([-10 10 -10 10 -10 10]), grid on, drawnow

 % find the peak frequencies in x, y, and z directions; i.e., find the
 % max in each direction of the normalized sum A2.
 % save these variables as A3, A4, and A5
%sumX = sum(abs(A2/M),1);
%[maxValueX,maxIndexX] = max(sumX(:));
%sumY = sum(A2/M,2);
%[maxValueY,maxIndexY] = max(sumY(:));
%sumZ = sum(A2/M,3);
%[maxValueZ,maxIndexZ] = max(sumZ(:));

% use max function to find max value and linear index in all of A2
A2_analysis = fftshift(A1)/49;
[maxVal,maxIndex] = max(A2_analysis(:));
% use ind2sub to obtain 3D array index of maxvalue 
[i,j,k]=ind2sub(size(A2), maxIndex);
% save x coordinate
A3 = ks(i);
% save y coordinate
A4 = ks(j);
%save z coordinate
A5 = ks(k);

%create an appropriate Gaussian filter and save it as A6
k0 = 0;
%filter = (exp(-tau*((Kx+2.5*A5).^2 + (Ky-A4*0.5).^2 +(Kz-1.15*A3).^2)))/M;
tau = 0.025;
filter = exp(-tau*((Kx+A3).^2 + (Ky+A4).^2 + (Kz+A5).^2));
isosurface(Kx,Ky,Kz,filter);


A6 = filter;
 

% Using the peak frequencies for the filtered signal, estimate the x, y, and z coordinates of the Kraken over time and save as A7, A8, A9
A7 = []; 
A8 = [];
A9 = [];
for c =1:49
    %obtain realization
    Un(:,:,:) = reshape(Kraken(:,c),n,n,n);
    %apply gaussain filter to fourier transform of realization
    Ufreqfilt = A6.*fftshift(fftn(Un));
    % obtain inverse fourier transform of filtered realization transform 
    Unfilt = ifftn(Ufreqfilt);

    % obtain max value element (this should be estimate of location)
    [maxValfilt,maxIndexfilt] = max(abs(Unfilt(:)));
    [i,j,k] = ind2sub(size(Unfilt),maxIndexfilt);
    % obtain and save x location estimate 
    A7(c) = x(i); 
    % obtain and save y location estimate
    A8(c) = y(j);
    % obtain and save z location estimate 
    A9(c) = z(k);
end

% Plot the location in x-y-z space over time for your report (not for the autograder)
% Plot the projection onto the x-y plane for your reprot (not for the autograder)
z_zeros = zeros(size(A7));
%z_zeros(:,:)= -10;
figure(2)
plot3(A7,A8,A9);grid on; 
hold on
title('Estimated Location of Vibration Source','FontSize',15)
xlabel('X');
ylabel('Y');
zlabel('Z');

figure(3)
plot(A7,A8); grid on; 
title('Estimated Location of Vibration Source on XY Plane','FontSize',15);
xlabel('X');
ylabel('Y');

% Save a table of the x-y-z coordinates for your report (not for the autograder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include all helper functions below since the autograder can only accept
% one file.