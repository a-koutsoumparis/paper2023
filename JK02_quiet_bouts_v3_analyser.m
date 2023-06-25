% % Quiescence Bout Analyser_V3 % %
% by JK
% % 
% -> Analyses Quiescence Bouts 
% Needed format of files: 'intensity_[STRAIN].txt'
% Supplied by 'IMAGE SUBTRACTOR_V2'

close all
clc;  % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
clear

% Define a starting folder.
topLevelFolder = uigetdir('','Where are your files?');
if topLevelFolder == 0
	return;
end
cd(topLevelFolder)

%common name of files to work with
common_name='intensity_*.txt';
list_info=dir(common_name);


% %Dialog Box: Worm strain
% Sub-Part: Get list
for list_num=1:length(list_info);
    list_info2=list_info(list_num).name;
    list_info3=cellstr(list_info2);
    list_info4=regexprep(list_info3,'intensity_','');
    list_info5=regexprep(list_info4,'.txt','');
    list_info_x(list_num,:)=list_info5;
        
end    
% Sub-Part: Dialog Box
[selection,ok] = listdlg('ListString',list_info_x,...
    'SelectionMode','single',...
    'ListSize',[300 150],...
    'Name','Strains:',...
    'PromptString','Choose Strain to analyse:',...
    'OKString','Go for it!',...
    'CancelString','Nahh, I changed my mind.');


%Dialog Box: Analysing specifications
prompt = {'Smoothing method (loess/rloess):','Smooth span (frames):','Threshold (0-1):','Minimum quiescence bout length (second)'};
dlg_title = 'Analyse';
num_lines = 1;
defaultans = {'loess','40','0.3','180'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

%input_excel_file='S1b_48h.xlsx'
strain=list_info_x{selection};
input_file=strcat('intensity_',strain,'.txt');
smooth_fit=answer{1,1};
span=str2num(answer{2,1});
thresh=str2num(answer{3,1});
min_bout_length=str2num(answer{4,1});

%get your data
    A1=readtable(input_file);
    A=table2array(A1);
    B=A(2:end,:);
    %round timestamp to seconds 
%     B(:,1)=round(B(:,1)/1000);
%     B(:,1)=round(B(:,1));
    B_index(:,1) = B(:,1);
%     B_index(:,1) = A(2:end,:);
    
% for ind=1:size(C,2); %temporary
for ind=1:size(B,2)-1;
    ind
    
    %find max/min of worm  
    B_smooth(:,ind) = smooth(B(:,ind+1), span, smooth_fit);
    
    %plot standardized raw data
    max1 = max(B_smooth(:,ind));
    min1 = min(B_smooth(:,ind));
    D=(B(:,ind+1)-min1)/(max1-min1);
    
    figure
    plot(B_index, D(:,1), 'x');
    hold on
    
    %plot smoothed data
    smooth_data =(B_smooth(:,ind)-min1)/(max1-min1);
    plot(B_index, smooth_data);    
    
    %build a table with values are above ('thresh') and below threshold ('0')
    thresh_data(:,1)=B_index;
    for ind2=1:size(B,1);    
        if  smooth_data(ind2) < thresh
            thresh_data(ind2,2) = 0;
        else
            thresh_data(ind2,2) = thresh;
        end
    end
    
    %plot threshold
    plot(B_index, thresh_data(:,2), 'k-')
    hold off
    
    %search quite bouts:
    i=1;
    j=0;
    k=1;
        while i < size(thresh_data,1)
            if thresh_data(i,2)==0
                %oh a quiet bout...better count the length
                i = i+1;
                j = j+1; 
            elseif and(thresh_data(i,2)==thresh, j>0)
                %that's where the shit is happening
                i = i+1;
                bout_length = thresh_data(i,1)-thresh_data(i-j,1);
                
                if bout_length > min_bout_length
                    bouts_quiet(k,ind) = bout_length;
                    k = k+1;
                else
                end
                j = 0;
                
            elseif and(thresh_data(i,2)==thresh, j==0)
                %continued thresh, no quiet bout...walk along
                i = i+1;
                j = 0;
                

                
            end
        end
    
if and(i==size(thresh_data,1), size(bouts_quiet,2)<ind)
%all quiescence bouts were to short, fill bouts_quiet-list
%with zeros
    bouts_quiet(:,ind) = 0;    
end

    %search activity bouts:
    i=1;
    j=0;
    k=1;
        while i < size(thresh_data,1)
            if thresh_data(i,2) == thresh
                %oh a wake bout...better count the length
                i = i+1;
                j = j+1;  
                
            elseif and(thresh_data(i,2)==0, j>0)
                %that's where the shit is happening
                i = i+1; 
                bout_length = thresh_data(i,1)-thresh_data(i-j,1);
                
                if bout_length > min_bout_length
                    bouts_activity(k,ind) = bout_length;
                    k = k+1;
                else
                end
                j = 0;
                
            elseif and(thresh_data(i,2)==0, j==0)
                %continued quiet, no activity...walk along
                i = i+1;
                j = 0;
            end
        end
        
if and(i==size(thresh_data,1), size(bouts_activity,2)<ind)
%all quiescence bouts were to short, fill bouts_quiet-list
%with zeros
    bouts_activity(1,ind) = B_index(end,1);    
end
        
    % % Calculating the difference/distance between "mean" and the normalized
    % % "50% value" of smoothed fit
    mean_data_smooth = mean(smooth_data);
    mean_dif = 0.5 - mean_data_smooth ;
    mean_distance = sqrt(mean_dif .* mean_dif);
    % % Calculating the standard deviation of raw data set
    standard_deviation=std(B_smooth(:,ind));      
    
    %tables for 'false positive' check (see below)
      tab_check(1,ind) = mean_distance;
      tab_check(2,ind) = standard_deviation;
      tab_check(3,ind) = min1;
    
    %check for 'false positive' sleep bouts in non-sleeping worms and
    %remove (you can adjust min1 to "<63" to include them often) standard
    %is "min1<59"
        if or(or(mean_distance>0.15 , standard_deviation>6), min1<62)
            bouts_quiet(:,ind)=bouts_quiet(:,ind);
            bouts_activity(:,ind)=bouts_activity(:,ind);
        else
            bouts_quiet(:,ind)=[0];
            bouts_activity(1,ind)=max(B_index);
            bouts_quiet(2:end,ind)=[0];
        end


end

%replace '0' with 'NaN' in  the bouts tables    
bouts_quiet(~bouts_quiet)=NaN;
bouts_activity(~bouts_activity)=NaN;

%write files
quiet_array=array2table(bouts_quiet);
writetable(quiet_array,[strain,'_all_quiet_bouts.txt'],'Delimiter','\t');

move_array=array2table(bouts_activity);
writetable(move_array,[strain,'_all_move_bouts.txt'],'Delimiter','\t');

