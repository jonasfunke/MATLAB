%% Process loaded RT-PCR-MX-P3005 data
%averagedData is averaged RT-PCR data (averaged measurement, well+filter)
%backgroundWells is array of background wells (relative numbering 1,2,3...)
%backgroundNumber is number of background wells
%averageBackground is array of background data of each filter
%averagedDataBack is averagedData with background substracted in non-background-wells
%averagedDataBackAvg is averagedDataBack with multiple identical wells averaged

%% subtract background, divide by background, multiply by mean background (to keep magnitude of intensity constant)

button='No';
while strcmp(button,'No')                                   %select background wells

    backgroundWellsString=inputdlg('comma-separated background well nrs');
    backgroundWells=str2num(backgroundWellsString{1})
    backgroundNumber=length(backgroundWells);

    figure
    for i=1:backgroundNumber                                %display selected background wells
        plot(averagedData(:,(backgroundWells(i)-1)*filterNumber+1:(backgroundWells(i)*filterNumber)))
        hold on
    end
    hold off
    
    button = questdlg('wells ok?','is it?' ,'No','Yes', 'Yes');
    clf
end

averageBackground=zeros(size(averagedData,1),filterNumber); %averaged backgrounds for filters
                                                            
for i=1:backgroundNumber                                    %calculate averaged backgrounds
    for j=1:filterNumber
        averageBackground(:,j)=averageBackground(:,j)+averagedData(:,(backgroundWells(i)-1)*filterNumber+j);
    end
end
averageBackground=averageBackground./double(backgroundNumber);

averagedDataBack=averagedData;
for i=1:wellNumber                                          %subtract average background from non-background data, divide by average background, multiply with mean average background
    for j=1:filterNumber
        if ismember(i,backgroundWells)
            averagedDataBack(:,(i-1)*filterNumber+j)=averagedData(:,(i-1)*filterNumber+j);
        else
            averagedDataBack(:,(i-1)*filterNumber+j)=((averagedData(:,(i-1)*filterNumber+j)-averageBackground(:,j))./averageBackground(:,j)).*mean(averageBackground(:,j));
        end
    end
end

%% Plot filter data

for i=1:filterNumber
    figure
    plot(averagedDataBack(:,i:filterNumber:wellNumber*filterNumber))
    legend(num2str([wellNames{1,i:filterNumber:wellNumber*filterNumber}]'))
    title(filters(i,:));
end

%% join multiple identical samples

button = questdlg('multiple wells?','multiple wells?' ,'No','Yes', 'Yes');

if strcmp(button,'Yes')
    numberDuplicatsString=inputdlg('number of identical samples');%select number of identical samples
    numberDuplicats=str2num(numberDuplicatsString{1})

    averagedDataBackAvg=zeros(size(averagedData,1),size(averagedData,2)/numberDuplicats);

    for i=1:wellNumber/numberDuplicats                          %average over multiples of same sample
        for j=1:filterNumber
            for k=1:numberDuplicats
                averagedDataBackAvg(:,(i-1)*filterNumber+j)=averagedDataBackAvg(:,(i-1)*filterNumber+j)+averagedDataBack(:,((i-1)*numberDuplicats+k-1)*filterNumber+j);
            end
            averagedDataBackAvg(:,(i-1)*filterNumber+j)=averagedDataBackAvg(:,(i-1)*filterNumber+j)./numberDuplicats;
        end
    end    
end


%% Plot averagedDataBackAvg data

if strcmp(button,'Yes')
    for i=1:filterNumber
        figure
        plot(averagedDataBackAvg(:,i:filterNumber:wellNumber*filterNumber/numberDuplicats))
        legend(num2str([1:wellNumber/numberDuplicats]'))
        title(filters(i,:));
    end
end
% %% normalize data 
% 
% for i=1:wellNumber/numberDuplicats                                                 %average over duplicats of same sample
%     for j=1:filterNumber
%         averagedData_back_avg_norm(:,(i-1)*filterNumber+j)=(averagedData_back_avg(:,(i-1)*filterNumber+j)-min(averagedData_back_avg(1:end,(i-1)*filterNumber+j)))./(max(averagedData_back_avg(1:end,(i-1)*filterNumber+j))-min(averagedData_back_avg(1:end,(i-1)*filterNumber+j)));
%     end
% end
% 
% 
% %% Plot filter data
% 
% for i=1:filterNumber
%     figure
%     plot(averagedData_back_avg_norm(:,i:filterNumber:wellNumber*filterNumber/numberDuplicats))
%     title(filters(i,:));
% end

%% save background-corrected data

file=fopen([pathname filename(1:length(filename)-38) '_minBackground.txt'], 'w');%open file to write
                                                    
for i=1:wellNumber*filterNumber                             %write 'wellNr_filterType' in each column header
        fprintf(file,[num2str(wellNames{1,i}) '_' wellNames{2,i} '\t']);
end
fprintf(file,'\n');
fclose(file);
                                                            %save averaged data to file
dlmwrite([pathname filename(1:length(filename)-38) '_minBackground.txt'], averagedDataBack, 'delimiter', '\t','-append')

%% save background-corrected averaged data    

if strcmp(button,'Yes')
    file=fopen([pathname filename(1:length(filename)-38) '_minBackgroundAvg.txt'], 'w');%open file to write                                                    

    for i=1:wellNumber*filterNumber/numberDuplicats             %write 'wellNr_filterType' in each column header
            fprintf(file,[num2str(ceil(double(i)/double(filterNumber))) '_' wellNames{2,i} '\t']);
    end
    fprintf(file,'\n');
    fclose(file);
                                                                %save averaged data to file
    dlmwrite([pathname filename(1:length(filename)-38) '_minBackgroundAvg.txt'], averagedDataBackAvg, 'delimiter', '\t','-append')
end

