classdef CalTimeRemain < handle
    %CalTimeRemain Calculate the remaining time to finish the for loop
    %   When first creating the class, a vector that indicating the total
    %   iteration times is required.
    %   Here is an example:
    %       Timer1 = CalTimeRemain([5,10]);
    %       % In this case, the vector [5,10] indicates that there are two
    %       % for-loops. Use method TimeRemain to estimate the remaining time
    %       for i = 1:5
    %           for j = 1:10
    %               do something;
    %               Timer1.timeRemain([i,j]); % estimate the remaing time
    %           end
    %       end
    %   
    %   By default, the notification pops on the command window will be
    %       "Time remaining to finish:  xxx min"
    %   Key words "Time remaining to finish:" can be changed by calling
    %   'setNotification'.
    %   The unit can be changed by calling 'setUnit'.
    %   
    %   After the loop finishes, please call "delete" to delete the class
    %
    %   By Ruizhi Cao, 02/11/2020, Biophotonics Lab, Caltech
    %   Modified on 03/08/2023
    
    properties(SetAccess = private, GetAccess = public)
        UseLastNTimer = []; 
        % use the time used among last N times that we call TimeRemain to
        % estimate the remaining time
        NotificationWord = 'Time remaining to finish: ';
        unitInUse = ' min';
    end
    
    properties(SetAccess = private, GetAccess = private)
        stopwatch;
        Units = {' hrs',' min',' sec'};
        scalef = 60;
        Record = []; % first row is the time used, 
                     % second row is the # iter finished
        Tidx = 1;
        charLast = 0;
        totalIterStep; % iterations for each for-loop
        totalIter;     % iterations in total
        IterCount      % iterations finished when the i-th idx changes
                       % e.g.  totalIterStep = [10,6,4] which means the
                       % first for-loop runs 10 iterations, the second
                       % for-loop runs 6 and the third runs 4.
                       % In this case, whenever the first index
                       % increases by 1, wev'e run 6*4 = 24 iterations.
                       % And IterCount(1) = 24. Similarly IterCount(2) = 4
    end
    
    methods
        function obj = CalTimeRemain(totalIter,UseNLoop)
            obj.totalIterStep = totalIter;
            obj.IterCount = zeros(1,length(totalIter)-1);
            tmp = 1;
            for idx = 1:(length(totalIter)-1)
                tmp = tmp*obj.totalIterStep(end-idx+1);
                obj.IterCount(end-idx+1) = tmp;
            end
            obj.totalIter = tmp*obj.totalIterStep(1);
            
            if nargin >1
                obj.UseLastNTimer = UseNLoop;
                obj.Record = zeros(2,UseNLoop);
            end
            obj.stopwatch = tic;
        end
        
        function setUnit(obj,unit)
            switch lower(unit)
                case {'hrs','hr'}
                    obj.unitInUse = obj.Units{1};
                    obj.scalef = 60*24;
                case 'min'
                    obj.unitInUse = obj.Units{2};
                    obj.scalef = 60;
                case 'sec'
                    obj.unitInUse = obj.Units{3};
                    obj.scalef = 1;
                otherwise
                    error('Supported uints are hrs, min and sec.');     
            end
        end
        
        function setNotification(obj,words)
            obj.NotificationWord = words;
        end
        
        function timeRemain(obj,StepNow)
            [xsize,~] = size(StepNow);
            if xsize == 1
                StepNow = StepNow.';
            end
            if ~isempty(obj.IterCount)
                iterNow = obj.IterCount*(StepNow(1:end-1,1) - 1) + StepNow(end);
            else
                iterNow = StepNow;
            end
            timeused = toc(obj.stopwatch);
            
            if ~isempty(obj.UseLastNTimer)
                idx2use = mod(obj.Tidx-1,obj.UseLastNTimer)+1;
                TimeBetween = timeused - obj.Record(1,idx2use);
                idxBetween = iterNow - obj.Record(2,idx2use);
                obj.Record(1,idx2use) = timeused;
                obj.Record(2,idx2use) = iterNow;
                timeremains = TimeBetween/idxBetween*(obj.totalIter - iterNow);
            else
                timeremains =  timeused/iterNow*(obj.totalIter - iterNow);
            end
            if obj.Tidx == 1
                charTime = num2str(timeremains/obj.scalef,'%.2f');
                fprintf([obj.NotificationWord,charTime,obj.unitInUse,'\n']);
                obj.charLast = length(charTime) + 5;
            else
                charTime = num2str(timeremains/obj.scalef,'%.2f');
                fprintf([repmat('\b',1,obj.charLast),charTime,obj.unitInUse,'\n']);
                obj.charLast = length(charTime) + 5;
            end
            obj.Tidx = obj.Tidx + 1;
        end
        
        function msgOut = readMsg(obj)
            totalTime = toc(obj.stopwatch);
            charTime = num2str(totalTime/obj.scalef,'%.2f');
            msgOut = ['Time used: ', charTime,obj.unitInUse];
        end
        
        function delete(obj) % varargout = delete(obj) not a valid destructor
            totalTime = toc(obj.stopwatch);
            charTime = num2str(totalTime/obj.scalef,'%.2f');
            lengthword = length(obj.NotificationWord);
            fprintf([repmat('\b',1,obj.charLast + lengthword),...
                'Finished. Total time used: ', charTime,obj.unitInUse,'\n']);
            clear obj.stopwatch
            delete(obj);
        end
    end
end

