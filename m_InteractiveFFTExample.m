%% Q3FoundationOfElectronMicroscopy
% By Robert J Scales HT 2020

clear
clc
close all

%% Settings

singlePixelPerAtom = false;
SavingLocation = 'C:\Users\rober\OneDrive - Nexus365\DPhil\PhD Lectures\Problem Sheets\Matlab Q3\Noise Without ColorScale';
NoisyTF = true;
MaskSideLength = 10;
ZoomedInSideLength = 30;
TakeSpotsOnFFT_TF = false;
FigureWaitTime = 0.5;

SavingLocation = horzcat(SavingLocation,'\');
MaskSideLength = RoundEven(MaskSideLength);
ZoomedInSideLength = RoundEven(ZoomedInSideLength);

%% STEM Image Creation

STEM = zeros(512,512);

gridSize = 10;
Y = 56;

N = nan(Y,1);
for i = 0:Y
    N(i+1) = (i*gridSize)-(i-1)+3;
    disp(N(i+1))
end

for j = 1:length(N)
    CurrRow = N(j);
    for k = 1:length(N)
        CurrCol = N(k);
        STEM(CurrRow,CurrCol) = 0.5;

        if singlePixelPerAtom == false
            STEM(CurrRow-1,CurrCol-1) = 0.125;
            STEM(CurrRow+1,CurrCol-1) = 0.125;
            STEM(CurrRow-1,CurrCol+1) = 0.125;
            STEM(CurrRow+1,CurrCol+1) = 0.125;
            STEM(CurrRow-1,CurrCol-1) = 0.125;
            STEM(CurrRow-1,CurrCol) = 0.25;
            STEM(CurrRow+1,CurrCol) = 0.25;
            STEM(CurrRow,CurrCol-1) = 0.25;
            STEM(CurrRow,CurrCol+1) = 0.25;
        end
    end
end
disp('Done creating array');

CurrRow = 256;
CurrCol = 256;
if singlePixelPerAtom == false
    STEM(CurrRow,CurrCol) = 1;
    STEM(CurrRow-1,CurrCol-1) = 0.25;
    STEM(CurrRow+1,CurrCol-1) = 0.25;
    STEM(CurrRow-1,CurrCol+1) = 0.25;
    STEM(CurrRow+1,CurrCol+1) = 0.25;
    STEM(CurrRow-1,CurrCol-1) = 0.25;
    STEM(CurrRow-1,CurrCol) = 0.5;
    STEM(CurrRow+1,CurrCol) = 0.5;
    STEM(CurrRow,CurrCol-1) = 0.5;
    STEM(CurrRow,CurrCol+1) = 0.5;
end

ZoomRange = 256-(ZoomedInSideLength/2):1:256+(ZoomedInSideLength/2);

filename_nonoise = horzcat(SavingLocation,'STEM_NoNoise.png');
filename_noise = horzcat(SavingLocation,'STEM_Noise.png');

data = uint8(STEM*255);
disp('Image creating');
% AutoFigure(STEM,SavingLocation,'STEM_NoNoise',true,300);
% AutoFigure(STEM,SavingLocation,'STEM_NoNoise_color',false,300);
figure('Name','Original STEM','WindowState','Maximized');
imshow(STEM);
axis image
imwrite(data, filename_nonoise,'bmp');
pause(FigureWaitTime);

%     AutoFigure(STEM(ZoomRange,ZoomRange),SavingLocation,'STEM_NoNoiseZoomed',true,300);
%     AutoFigure(STEM(ZoomRange,ZoomRange),SavingLocation,'STEM_NoNoiseZoomed_color',false,300);
figure('Name','STEM Zoomed','WindowState','Maximized');
imshow(STEM(ZoomRange,ZoomRange));
imwrite(STEM(ZoomRange,ZoomRange),horzcat(SavingLocation,'\ZoomedIn_NoNoise.bmp'),'bmp');
pause(FigureWaitTime);

data = uint8(STEM);


if NoisyTF == true
    disp('Adding Noise');
    I = imread(filename_nonoise);
    STEM_noise = imnoise(I,'gaussian',0,0.005);
    NoisyData = uint8(STEM_noise);
%     NoisyData = STEM_noise;
%     AutoFigure(NoisyData,SavingLocation,'STEM_Noise',true,300);
%     AutoFigure(NoisyData,SavingLocation,'STEM_Noise_color',false,300);
    figure('Name','Noisy STEM','WindowState','Maximized');
    imshow(STEM_noise);
    imwrite(uint8(STEM_noise),filename_noise,'bmp');
    pause(FigureWaitTime);

%     AutoFigure(NoisyData(ZoomRange,ZoomRange),SavingLocation,'STEM_NoiseZoomed',true,300);
%     AutoFigure(NoisyData(ZoomRange,ZoomRange),SavingLocation,'STEM_NoiseZoomed_color',false,300);
    figure('Name','Noisy STEM Zoomed','WindowState','Maximized');
    imshow(NoisyData(ZoomRange,ZoomRange));
    imwrite(NoisyData(ZoomRange,ZoomRange),horzcat(SavingLocation,'\ZoomedIn_Noise.bmp'),'bmp');
    pause(FigureWaitTime);
    
    data = uint8(STEM_noise);
end

commandwindow;
disp('Done Making STEM Images');
close all;
%% FFT Section
clc
Y1 = (fft2(data));
Y2 = fftshift(fft2(data));
Y = Y2;

% FFT = figure('Name','FFT','WindowState','Maximized');
AutoFigure(abs(Y),SavingLocation,'FFT',false,300);
% savename = horzcat(SavingLocation,'\FFT.png');
% imagesc(abs(Y));
% axis image
% exportgraphics(FFT,savename,'Resolution',300);
pause(FigureWaitTime);

AutoFigure(abs(log2(Y)),SavingLocation,'FFTLog2',false,300);
% Log2FFT = figure('Name','Log2 FFT','WindowState','Maximized');
% savename = horzcat(SavingLocation,'\FFTLog.png');
% imagesc(abs(log2(Y)));
% axis image
% exportgraphics(Log2FFT,savename,'Resolution',300);
pause(FigureWaitTime);

if TakeSpotsOnFFT_TF == true
%     fh = findall(groot, 'Type', 'figure', 'Name', 'abc');
    fh = findobj( 'Type', 'Figure', 'Name', 'FFT' );
    figure(fh);
end

MaskRelPxlRng = (-(MaskSideLength/2)+1):1:(MaskSideLength/2);

SafeZonePoints = nan(1,2);
EndSpotZone = false;
i = 1;
while EndSpotZone == false
    try SafeZonePoints(1,:) = ginput(1);
        SafeZonePoints = floor(SafeZonePoints);
        SafeZonePointX = SafeZonePoints(1,1);
        SafeZonePointY = SafeZonePoints(1,2);
        hold on
        plot(SafeZonePointX,SafeZonePointY,'rx','MarkerSize',12);

        [A,B] = meshgrid(SafeZonePointX+(MaskRelPxlRng),SafeZonePointY+(MaskRelPxlRng));
        c=cat(2,A',B');
        CurrSafeZone=reshape(c,[],2);
        InputRow = ((i-1)*100)+1;
        if i ==1
            SafeZones = CurrSafeZone;
        else
            SafeZones = vertcat(SafeZones,CurrSafeZone);
        end
        commandwindow;
        str = input('If finished type Y, if type in anything: ','s');
        if strcmp(str,'Y') == true || strcmp(str,'y')
            EndSpotZone = true;
        else
            i = i+1;
        end
    catch
        disp('Do not press enter!');
    end
end

commandwindow;
disp('Done FFT Section');
close all;
%% Masking

Mesh = zeros(size(Y));

for i=1:size(SafeZones,1)
    CurrentPixel = SafeZones(i,:);
    Mesh(CurrentPixel(2),CurrentPixel(1)) = 1;
end

MaskedY = Y.*Mesh;

AutoFigure(abs(log2(MaskedY)),SavingLocation,'FFTLog2Masked',false,300);
pause(FigureWaitTime);
AutoFigure(abs(MaskedY),SavingLocation,'FFTMasked',false,300);
pause(FigureWaitTime);

% FFTLog2Masked = figure('Name','Masked FFT','WindowState','Maximized');
% savename = horzcat(SavingLocation,'\FFTLog2Masked.png');
% imagesc(abs(log2(MaskedY)));
% axis image
% exportgraphics(FFTLog2Masked,savename,'Resolution',300);

commandwindow;
disp('Done Masking FFT');
close all;
%% Inverse FFT

inverseFFT = ifft2(MaskedY);
AutoFigure(abs(inverseFFT),SavingLocation,'invFFT',true,300);
pause(FigureWaitTime);
AutoFigure(abs(inverseFFT),SavingLocation,'invFFT_color',false,300);
pause(FigureWaitTime);
% figure('Name','Inverse FFT Post-Masks','WindowState','Maximized');
% savename = horzcat(SavingLocation,'\invFFT.bmp');
% % imshow(uint8(abs(inverseFFT)));
% imagesc(abs(inverseFFT));
% colobar
% imwrite(uint8(abs(inverseFFT)), horzcat(SavingLocation,'\invFFT.bmp'),'bmp');
% axis image

ZoomedIn = abs(inverseFFT(ZoomRange,ZoomRange));
AutoFigure(ZoomedIn,SavingLocation,'ZoomedIn_invFFT',true,300);
pause(FigureWaitTime);
AutoFigure(ZoomedIn,SavingLocation,'ZoomedIn_invFFT_color',false,300);
pause(FigureWaitTime);
% figure('Name','Inverse FFT Post-Masks Zoomed','WindowState','Maximized');
% imshow(ZoomedIn);
% imwrite(ZoomedIn, horzcat(SavingLocation,'\invFFTZoomed.bmp'),'bmp');
disp('Done with Matlab code');

%% Function

function y=RoundEven(x)
    y = 2*floor(x/2);
end

function AutoFigure(Data,SavingLocation,SaveID,imshowTF,Resolution)

    figure('Name',SaveID,'WindowState','Maximized');
    savename = sprintf('%s%s%s',SavingLocation,SaveID,'.png');
    if imshowTF == true
        imshow(uint8(Data));
        colorbar
        axis image
        imwrite(uint8(Data), savename,'png');
        fprintf('Saved figure "%s" with imwrite method\n',SaveID);
    else
        imagesc(Data);
        colorbar
        axis image
        exportgraphics(gcf,savename,'Resolution',Resolution);
        fprintf('Saved figure "%s" with exportgraphics method\n',SaveID);
    end

end