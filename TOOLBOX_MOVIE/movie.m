classdef movie < handle
    %TRACE_SELECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frames; %all images to be read
        N_read = 300; % number of images to read in one portion
        counter = 1; % internal counter for reading the movie
        
        sequence; % sequence to be read, e.g. 101010
        first; % first image to be read
        last; % last image to be read
        
        pname; %pathname of file location
        fname; %filename of file
        
        sizeX; % number of pixel in X-dimension
        sizeY; % number of pixel in Y-dimension
        mov_length; % number of frames in thw whole movue
        
        info; %fits info
        h_min; %minimal heigth for peak fidning
        
        input; % 0=fits, 1=tiff-stack
        fnames; % cell with all filenames, only for tiff-stack
    end
    
    methods
        %constructor
        function obj = movie(pname, fname, first, last, sequence)
            obj.sequence = sequence;
            obj.pname = pname;
            obj.fname = fname;
            
            if strcmp(fname(end-2:end), 'tif')
                obj.input = 1; % read tiff data
            else
                obj.input = 0; % read fits data
            end
            
            if obj.input == 1 % tiff-stack
                tmp = dir([pname filesep '*.tif']);
                obj.fnames = {tmp.name};
                obj.info = 'Tiff-Stack info';
                obj.sizeX = size(imread([pname filesep obj.fnames{1}]),2);
                obj.sizeY = size(imread([pname filesep obj.fnames{1}]),1);
                obj.mov_length = length(obj.fnames);
            else % fits
                obj.info = fitsinfo([obj.pname obj.fname]);
                obj.sizeX = obj.info.PrimaryData.Size(1); 
                obj.sizeY = obj.info.PrimaryData.Size(2);
                obj.mov_length = obj.info.PrimaryData.Size(3);
            end
            
            
            obj.first = first;
            if last == -1
                obj.last = obj.mov_length;
            else
                obj.last = last;
            end
            obj.frames = obj.getFrames(obj.sequence, obj.first, obj.last);
            
        end
        
        %generate a list of images to be read from the movie
         function frames = getFrames(obj, sequence, first, last)
            frames = [];
            tmp = first:last;
            for i=1:length(tmp)
               if(  sequence(  mod(i-1,size(sequence,2))+1 )  )
                  frames = [frames tmp(i) ];       
               end    
            end
         end
         
         % initialize the counter         
         function  initRead(obj)
              obj.counter = 1;
         end
          
         %reads the next N_read frames
         function [tmp, frames_out, go_on] = readNext(obj)
            read_stop = obj.counter+obj.N_read-1; % last frame to read
            go_on = 1;
            if read_stop >= length(obj.frames)
                read_stop = length(obj.frames);
                go_on = 0;
            end
            
            
            tmp = zeros(obj.sizeX, obj.sizeY, read_stop-obj.counter+1);
            frames_out = obj.frames(obj.counter:read_stop);
            
            display(['Reading frame ' num2str(frames_out(1)) ' to ' num2str(frames_out(end))])
            for i=1:length(frames_out)
                if obj.input == 1 % tif
                    tmp(:,:,i) = double(imread([obj.pname filesep obj.fnames{frames_out(i)}]));
                else
                    tmp(:,:,i) = fitsread([obj.pname obj.fname],  'Info', obj.info, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [frames_out(i) frames_out(i)] });  % fits               
                end
            end
            
            obj.counter = obj.counter + length(frames_out); 
         end
         
         
         
         
         % trace the movie 
         function [ traces, itraces, avg_frame ] = trace_movie(obj, h_min, r_find, r_integrate, min_length )
            traces = cell(0,1);
            itraces = cell(0,1);
            avg_frame = zeros(obj.sizeX, obj.sizeY);
            
            go_on = 1;
            obj.initRead;
            N = 0;
                     
            display('Tracing movie... please wait')
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                [traces, itraces] = append_traces(movie, traces, itraces, frames, h_min, r_find, r_integrate, min_length);
            
                avg_frame = avg_frame + sum(movie,3);
                N = N + length(frames);
            end
            avg_frame = avg_frame ./ N; 
            display('Done tracing movie.')
         end
         
         % integrate specific reagions in movie
         function [itraces ] = traces_movie_position(obj, positions )
            itraces = cell(0,1);
            go_on = 1;
            obj.initRead;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                itraces = append_traces_to_position(movie, itraces, frames, positions, r_integrate);
            end
         end
         
         % determine peak-finding threshholds
         function [h_min ] = get_h_min(obj, r_find)
            if obj.input == 1 % tiff
                img = double(imread([obj.pname filesep obj.fnames{obj.frames(1)}]));
            else % fits
                img = fitsread([obj.pname obj.fname],  'Info', obj.info, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(1)] }); % read first frame                
            end
            p = find_peaks2d(img, r_find*2, min(img(:)), 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
            
            close all
            figure('units','normalized','outerposition',[0 0 1 1])
            img_mean = mean(img(:));
            img_std = std(img(:));            
            
            p_7std = p(find(p(:,4)>=7*img_std), :); % this is just an estimate, #of peaks found may vary since peak_find algorithm return height as int not double
            p_5std = p(find(p(:,4)>=5*img_std), :);
            p_3std = p(find(p(:,4)>=3*img_std), :);

            % plot
            subplot(1, 2, 1)
            imagesc(img), colorbar, axis image, colormap gray, hold on
            h(1) = plot(p_3std(:,1)+1, p_3std(:,2)+1, 'ro');
            h(2) = plot(p_5std(:,1)+1, p_5std(:,2)+1, 'go');
            h(3) = plot(p_7std(:,1)+1, p_7std(:,2)+1, 'bo');
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
                  
            
            subplot(1, 2, 2)
            xhist = min(p(:,4)):5:max(p(:,4));
            n = hist(p(:,4), xhist);
            semilogy(xhist, sum(n)-cumsum(n)), hold on
            h(1) = vline(3*img_std, 'r');
            h(2) = vline(5*img_std, 'g');
            h(3) = vline(7*img_std, 'b');
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
            set(gca, 'XLim', [xhist(1) xhist(end)])
            xlabel('Minimal height'), ylabel('# of peaks found')
            axis square
            
                       
            % promp
            options.WindowStyle='normal';
            prompt={'Enter min heigth (default=5*sigma):'};
            def={num2str(round(5*img_std))};
            threshold = inputdlg(prompt, strcat('Enter threshold:'), 1, def, options);
            h_min = str2double(threshold(1));
            close all
            
            obj.h_min = h_min;
            
         end


    end
    

    
    
    
end
        