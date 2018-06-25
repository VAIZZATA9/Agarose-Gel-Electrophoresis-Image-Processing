clear all
close all

%%%%% NOTE: DM denotes Daniel Mar edit %%%%%
%%%%% JLC denotes Jackson Chin comments; edits are listed with date %%%%%

%DM "filePath" looks for a folder that is in the same location as this MatLab
%file. "fileName" looks for the image file within that folder that you want
%to analyze. Note that you may have to change the extension depending on
%the format of the file (e.g. JPG vs PNG)

filePath = './Gel Images/';
fileName = '06-21-18 CryoCore21 Liver Samples PIXUL mixed rpt.JPG';

im = imread([filePath,fileName]);
im = rgb2gray(im);
im = double(im)/256;

stds = 1.0; % standard deviations to report about the mean. Originally 1.0
o = 0.5;    % lambda used in box-cox transformation, must be positive, good range probably from .1 to 1
            % originally o = 0.25

% JLC 5/3/2018: Added another base pair at 0 to improve range of window
calBasePairs = [1650,1000,850,650,500,400,300,200,100]; % base pairs in ladder

% DM Creating Arrays for waterfall plot and 200-600bp percentage calculation
ArrayGrid = []; 
ArrayGrid2 = []; 
ArrayGridNormalized = [];
ArrayGridNorm2 = [];
BPwindow = [];
BPtotal = [];
BPmean = [];
BPplus = [];
BPminus = [];

% JLC 4/4/2018: Replaced user input box with an automatic algorithm that
% fits a box to the data. Uses the function listed at the bottom of the
% code to find this box.
g = BoxFinderv2(im);

% JLC - Crops input image file to the rectangle specified by getBox.
im = im(g(1):g(1)+g(3)-1,g(2):g(2)+g(4)-1);
im = impyramid(im,'expand'); % might want to remove this if image already high-res

% JLC - Plots cropped image in grayscale. 
figure(1);
set(1,'pos',[1036,77,793,605]);
set(1,'paperpositionMode','auto');
subplot(2,3,[1,2]);
hold off;
imagesc(im);
colormap(bone)
hold on;
cmap = get(gca,'colororder');
xlabel('pixels');
ylabel('pixels');
title('Original PAGE Image');

% JLC - Computes the mean intensity of every column of pixels within a well,
% then returns the index of said column as I. 'dummy' is not used.
[dummy,I] = max(mean(im(1:end,1:150),1));  %DM edit
% DM note: need to narrow down by row and column so that it selects the
% correct ladder. eg. mean(im(row1:row2,col1:col2))
% JLC - Uses correlation coefficients to identify wells in figure.
% Correlation coefficients identify linear dependence between output value
% and relative vertical position
c = corrcoef(im);
c = c(I,:);
I = find(c>=0.875);  %DEFAULT .975 - this controls the width of the ladder, change if ladder bounds look wrong
% JLC -  Editted to 0.875 to better fit well width

% JLC - Provides reference bars for well width, and calculates the overall
% well width
plot([I(1),I(1)],ylim,'color',cmap(2,:));
plot([I(end),I(end)],ylim,'color',cmap(2,:));

wellWidth = I(end)-I(1);

% JLC - Plots the mean intensity of each pixel row in the far left well
cal = mean(im(:,I),2);
figure(2);
subplot(1,3,[1,2]);
hold off;
plot(cal);
hold on;

% JLC 3/27/18: Editted smoothing to 3 to match calibration curve closer.
% Dots on final plot now pull from cal curve, not cal2 curve
% cal2 = smooth(cal,3);
% JLC 4/25/18: Swapped 'MinPeakHeight' to "MinPeakProminence' to account
% for lower peaks in some plots
[pks,locs] = findpeaks(cal,'MinPeakProminence',0.025); %DM cal changed to cal2
%DM added MinPeakHeight to further remove false peaks (based on visual
%cutoff. May need to adjust this number manually, depending on the gel

% JLC 3/27/18: Uncertain of following 2 lines of code. Needed to make locs
% and calBasePairs share common dimensions. Is there a better way to do
% this? Are early peaks false positives?
pks = pks(end-length(calBasePairs)+1:end);
locs = locs(end-length(calBasePairs)+1:end);
plot(locs,cal(locs),'.');

% JLC - Plots abbreviated peak locations list against calibration base pairs
subplot(1,3,3);
hold off;
plot(locs,calBasePairs,'.');
hold on;

% JLC - Stores base pair indices as their corresponding peak locales and
% interpolates
bpIndeces = locs(1):locs(end);
bpValues = interp1(locs,calBasePairs,bpIndeces,'cubic');

% JLC - Plots interpolated base pair curve over scatter plot above
plot(bpIndeces,bpValues);

% JLC - Establishes limits to plot containing base pair indeces and their
% respective values
figure(1);plot(xlim,[locs(1),locs(1)],'color',cmap(2,:));
figure(1);plot(xlim,[locs(end),locs(end)],'color',cmap(2,:));

% JLC - Plots Base Pair values against their respective indeces
subplot(2,3,[4,5]);
hold off;
pcolor(1:size(im,2),bpValues,im(bpIndeces(1):bpIndeces(end),:));
shading flat;
axis xy;
xlabel('pixels');
ylabel('base pairs');
title('PAGE Image Rescaled Using Base Pair Calibration');
hold on;

%%

% JLC - Initializes neccesary variables needed to run loop
g = 0;
n = 1;
imageCaption = [];

% JLC - While loop that processes each selected data point. Runs only if g
% is not empty, which happens if the user presses 'enter' instead of
% selecting a new point.
while ~isempty(g)
    % JLC - Selects figure with picture of wells, and requests user to
    % input a point in the middle of a well. Sets I to be equivalent to indeces
    % of all points within the selected well.
    figure(1);
    subplot(2,3,[1,2]);
    disp('Select a point in the center of the well you want to analyze. Return if done...');

    g = round(ginput(1));
    
    if ~isempty(g)
        I = g(1)-wellWidth:g(1)+wellWidth;
        I = I(I>=1);
        I = I(I<=size(im,2));
        
        % JLC 3/28/18: Altered below code to better match MatLab syntax.
        % imTemp below is a cropped form of im, which focuses only on the
        % indeces of the base pairs and the well width
        imTemp = im(bpIndeces,I);
        
        % JLC - Buffers the following m and I vectors into matrices. Each
        % matrix has wellWidth vectors with length wellWidth - 1
        m = buffer(mean(imTemp),wellWidth,wellWidth-1,'nodelay');
        I = buffer(I,wellWidth,wellWidth-1,'nodelay');
        
        % JLC - Calculates the product of each column in the m matrix and
        % returns the values as a row. Finds the location of the maximum
        % value, and returns the index as II.
        [dummy,II] = max(prod(m));
        
        % JLC 3/28/18: Adjusted syntax below. Updates imTemp to include
        % only the column of values containing the maximum intensity.
        imTemp = im(bpIndeces,I(:,II));
        
        % JLC - Uses Cox-Box (classical curve-fitting algorithm) to match  a
        % normal curve to the generated data. Plots adjusted x-axis against
        % the mean intensity for each pixel in the column with the highest
        % intensity in each well
        x = (bpValues(:).^o-1)/o;
        y = mean(imTemp,2);
        figure(n+2);
        set(gcf,'pos',[8   407   665   278]);
        % subplot(1,3,[1,2]);
        set(gcf,'paperpositionMode','auto');
        hold off;
        
        ArrayGrid(n,:) = (x*o+1).^(1/o);  %DM WATERFALL
        ArrayGrid2(n,:) = y;  %DM WATERFALL

        
        %for loop1 = 1:length(y)                     %DM WATERFALL
        %ArrayGrid(loop1,:) = ((x*o+1).^(1/o));    %DM WATERFALL
        %end                                     %DM WATERFALL
        [dummy,I] = max(y);
        ft = fittype( 'a+b*exp(-(x-c)^2/(2*d^2))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.1 .5 x(I) 1/o];    % might want to play with these if fit looks really bad
                                                % default: opts.StartPoint
                                                % = [0.5 .5 x(I) 1/o]
        
        [fitresult, gof] = fit( x, y, ft, opts );
        
        % JLC 4/19/18: Editted base pair intensity values to normalize and
        % start at 0. ylim sets axis to fit adjusted curve.
        y = y - min(fitresult(x));
        plot((x*o+1).^(1/o),y);
        ylim([0 1.25*max(y)])
        hold on;
        
        ArrayGrid2(n,:) = fitresult(x); %DM Use this line to have the waterfall plot show the best-fit instead of the raw data
        AGbaseline = min(ArrayGrid2(n,:));
        %ArrayGridNormalized(n,:) = fliplr(ArrayGrid2(n,:)) - AGbaseline; 
        %DM 11-09-17: add the division to see if figure displays better
        %The following code refines how the waterfall plot is displayed.
        %Basically, we fix the starting point of each curve at 1 by
        %dividing every value by that starting value.
        ArrayGridNormalized(n,:) = (fliplr(ArrayGrid2(n,:)) - AGbaseline);
        ArrayGridNorm2(n,:) = ArrayGrid2(n,:)/ArrayGrid2(n,end);
        ArrayGridNorm2(n,:) = fliplr(ArrayGridNorm2(n,:));
        ArrayGridNorm2(n,:) = ArrayGridNorm2(n,:) - 1;
                
        %Use these lines to determine ratio of 200bp-600bp area to 100-1000bp area
        %NOTE: must rescale with each graph
        BP1 = round((length(ArrayGridNormalized)/9),0);
        BP2 = round((5*length(ArrayGridNormalized)/9),0);
        BP3 = length(ArrayGridNormalized);
        BPwindow(n) = trapz(BP1:BP2, ArrayGridNormalized(n,BP1:BP2));
        BPtotal(n) = trapz(1:BP3,ArrayGridNormalized(n,1:BP3));
        BPratio(n) = BPwindow(n)/BPtotal(n);
        %The problem with these calculations is that we don't have a good
        %baseline value for the area. My solution is to create a new
        %ArrayGrid called ArrayGridNormalized (flipped left-to-right to
        %correct for orientation), in which we subtract out the minimum
        %value of each row (effectively setting baseline at 0 for each
        %gel).
        
        % JLC 4/19/18: Plots best-fit Gaussian curve. Editted to normalize
        % and start curve at 0.
        plot((x*o+1).^(1/o),fitresult(x)-min(fitresult(x)),'linewidth',2);
        
        m = (fitresult.c*o+1)^(1/o);
        r = ([fitresult.c-stds*fitresult.d,fitresult.c+stds*fitresult.d]*o+1).^(1/o);
        yl = ylim;
        %Plotting vertical lines in graph
        plot([m,m],yl,'-','color',cmap(3,:));
        plot([r(1),r(1)],yl,'-','color',cmap(4,:));
        plot([r(2),r(2)],yl,'-','color',cmap(4,:));
        %
        ylim(yl);
        
        %DM 02-09-18: add output for mean +/- stdev
        BPmean(n) = m;
        BPplus(n) = r(2);
        BPminus(n) = r(1);
        
        %DM: use the following (commented-out) line of code if you want to display the standard
        %deviation metrics:
        %s = sprintf('well %d, mean: %0.1f, +-%0.1f std: %0.1f, %0.1f, BPratio: %%%0.1f',n,m,stds,r(1),r(2),BPratio(n)*100); 
        s = sprintf('well %d, mean: %0.1f, BPratio: %%%0.1f',n,m,BPratio(n)*100);
        xlabel('base pairs');
        ylabel('normalized intensity');
        title(s);
        legend('image intensity','Gaussian fit with Box-Cox','mean',sprintf('+-%0.1f std range',stds),'Location','EastOutside');
        n = n+1;
        
        imageCaption = sprintf('%s%s\n',imageCaption,s);
        
        figure(1);
        subplot(2,3,[4,5]);
        plot([g(1)-wellWidth/2,g(1)+wellWidth/2],[m,m],'-','color',cmap(2,:),'linewidth',2);
        plot([g(1)-wellWidth/4,g(1)+wellWidth/4],[r(1),r(1)],'-','color',cmap(2,:),'linewidth',2);
        plot([g(1)-wellWidth/4,g(1)+wellWidth/4],[r(2),r(2)],'-','color',cmap(2,:),'linewidth',2);
        
    end
end

% JLC - Displays data derived from each well on the original image of the
% wells
subplot(2,3,[4,5]);
text(max(xlim)+0.05*diff(xlim),max(ylim),imageCaption);

%DM Export BP ratios to text file
BPfile = fopen('BPratios.txt','a');
fprintf(BPfile,'%s',fileName);
fprintf(BPfile,'\r\n');
for BPprintloop = 1:n-1
fprintf(BPfile,'%%%0.1f \t',BPratio(BPprintloop)*100);
end
fprintf(BPfile,'\r\n');
fclose('all');

%DM Export BP means and range to text file
BPfile2 = fopen('BPmeans-range.txt','a');
fprintf(BPfile2,'%s',fileName);
fprintf(BPfile2,'\r\n');
for BPprintloop = 1:n-1
fprintf(BPfile2,'%0.1f \t',BPminus(BPprintloop));
end
fprintf(BPfile2,'\r\n');
for BPprintloop = 1:n-1
fprintf(BPfile2,'%0.1f \t',BPmean(BPprintloop));
end
fprintf(BPfile2,'\r\n');
for BPprintloop = 1:n-1
fprintf(BPfile2,'%0.1f \t',BPplus(BPprintloop));
end
fprintf(BPfile2,'\r\n');
fclose('all');


%%

[dummy,fn] = fileparts(fileName);

print ([fn,' All Wells'],'-f1','-djpeg90','-r250');

for loop1 = 1:(n-1)
    figNum = loop1+2;
    print (sprintf('%s Well %02d',fn,loop1),sprintf('-f%d',figNum),'-djpeg90','-r250');
end

%%%%%DM WATERFALL%%%%%
figure(100);

    %To Display non-normalized version: 
    %ArrayGrid2 = fliplr(ArrayGrid2); % To correct orientation
    %ArrayPlot = waterfall(ArrayGrid2);
    
%To Display normalized version:
ArrayPlot = waterfall(ArrayGridNorm2);
set(ArrayPlot,'LineWidth',2); %DM change width of waterfall lines here

%DM Set x-axis labels so that they accurately reflect base length instead
%of array length, which varies by graph
%Also have the option to fix z-axis or have it scale by max value
xlabelquarter = length(ArrayGridNorm2)/4;
xlabelhalf = length(ArrayGridNorm2)/2;
xlabelthreequarters = 3*length(ArrayGridNorm2)/4;
xlabelfull = length(ArrayGridNorm2);
zmax = max(ArrayGridNorm2(:)) + 0.01;
set(gca,'xlim',[0 xlabelfull]);
set(gca,'XTick',[0 xlabelquarter xlabelhalf xlabelthreequarters xlabelfull]);
set(gca,'XTickLabel',[100 325 550 775 1000]);
%set(gca,'zlim',[0.97 1.2]);
set(gca,'zlim',[0 zmax]);

print([fn,' Well Array'],'-f100','-djpeg90','-r250');       %DM WATERFALL