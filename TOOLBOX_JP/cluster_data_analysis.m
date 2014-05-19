%% select data
%selects data to open
%filename{N} are filenames of files
%pathname is pathnames of files

[filename pathname]=uigetfile('*','feed me cluster data human','MultiSelect','on');


%% load data and analyze
%opens file, copies data into parameters and currentData, analyzes,
%then opens next file
%parameters are cells with parameter strings
%currentData is array of traceData x structures, with empty lines between
%traces omited

numStructures=15;                       %number of structures simulated
timeSteps=7200+1;                       %number of steps printed
traces=10;                              %number of traces per structure and parameter set

numParameterLines=28;                   %current output length before data starts
temperatureParameterLine=20;            %location of temperature value in beginning output
parameters=cell(1+numStructures+numParameterLines,1);%saves parameters text, each cell one line

numFiles=size(filename,2);              %number of files loaded
currentData=zeros(timeSteps*traces,numStructures+1);%temporary storage current file

results=cell(numFiles+1,1+1+numStructures);

for i=1:numFiles
    
    currentFile=fopen([pathname filename{i}],'r');      %load current file
    filename{i}
    for j=1:1+numStructures+numParameterLines
        parameters{j}=fgetl(currentFile);               %load parameter text
    end

    for j=0:(traces-1)                                  %load simulation data
        fgetl(currentFile);                             
        for k=1:timeSteps
            currentData(k+j*timeSteps,:)=fscanf(currentFile,'%f',numStructures+1);
        end
    end   

    fclose(currentFile);
    
    %analysis part
    
    %fill results with starting data
    results{i+1,1}=filename{i};
    tempLocation=strfind(filename{i}, '_t');
    results{i+1,2}=filename{i}(tempLocation+2:tempLocation+3);
    
    for j=1:numStructures
        results{i+1,j+2}=traces;
    end
    %check if structures are unfolded at end of trace
    for j=1:traces
        for k=1:numStructures
            if sum(currentData(j*timeSteps-30:j*timeSteps,k+1))~=0
                results{i+1,k+2}=results{i+1,k+2}-1;
            end
        end
    end
    
end

%write column names from parameters text block
results{1,1}='filename';
results{1,2}='0temperature';
results{1,numStructures+3}=0;

for i=1:numStructures
    tempLocation=strfind(parameters{i+1}, '/');
    results{1,i+2}=parameters{i+1}(tempLocation(size(tempLocation,2))+1:size(parameters{i+1},2));
end
    
%sort according to temperature and number of folded structs
for i=1:numFiles
    results{i+1,numStructures+3}=sum([results{i+1,3:numStructures+2}]);
end
results=sortrows(results,[2,numStructures+3]);
results{1,numStructures+3}='total_folded';


%% save results to  results.txt

currentFile=fopen([pathname filesep 'results.txt'],'w');

for i=1:numFiles+1
    for j=1:numStructures+3
        if ischar(results{i,j})
            fprintf(currentFile,'%s\t',results{i,j});
        else
            fprintf(currentFile,'%f\t',results{i,j});
        end
    end
    fprintf(currentFile,'\n');
end

fclose(currentFile);


%% visualisation

resultsMat=zeros(numFiles,numStructures);

for i=1:numFiles
    for j=1:numStructures
        resultsMat(i,j)=results{i+1,j+2};
    end
end

resultsMat=resultsMat(:,[1,2,4,3,7,8,11,5,12,9,10,13,6,15,14]);

currentT=46;
previousPos=1;
currentPos=1;
while currentT==str2num(results{currentPos+1,2})
    currentPos=currentPos+1;
end

fig46=figure;
barPlot46=bar3(resultsMat(previousPos:currentPos-1,:));
for k = 1:length(barPlot46)
    zdata = get(barPlot46(k),'ZData');
    set(barPlot46(k),'CData',zdata,...
             'FaceColor','interp')
end
view(0, 90);
title('46 C');
set(gca, 'XTickLabel',{'v1','v2','v4','v3','v7','v8','v14','v5','v15','v11','v13','ssp','v6','RR','RRv3'})
xlabel('structure')
ylabel('parameter Set')

previousPos=currentPos;
currentT=50;

while currentT==str2num(results{currentPos+1,2})
    currentPos=currentPos+1;
end

fig50=figure;
barPlot50=bar3(resultsMat(previousPos:currentPos-1,:));
for k = 1:length(barPlot50)
    zdata = get(barPlot50(k),'ZData');
    set(barPlot50(k),'CData',zdata,...
             'FaceColor','interp')
end
view(0, 90);
title('50 C');
set(gca, 'XTickLabel',{'v1','v2','v4','v3','v7','v8','v14','v5','v15','v11','v13','ssp','v6','RR','RRv3'})
xlabel('structure')
ylabel('parameter Set')

previousPos=currentPos;

fig54=figure;
barPlot54=bar3(resultsMat(previousPos:numFiles,:));
for k = 1:length(barPlot54)
    zdata = get(barPlot54(k),'ZData');
    set(barPlot54(k),'CData',zdata,...
             'FaceColor','interp')
end
view(0, 90);
title('54 C');
set(gca, 'XTickLabel',{'v1','v2','v4','v3','v7','v8','v14','v5','v15','v11','v13','ssp','v6','RR','RRv3'})
xlabel('structure')
ylabel('parameter Set')

























