classdef beadprocessor
    %%BEADPROCESSOR Process probe particle microrheology data
    %
    %%%%%%%%%%%%%%%%%%%%%%%% TO DO: %%%%%%%%%%%%%%%%%%%%%%%%
    % Write some documentation
    % Implement:
    %  Data organisation
    %   Store all centres in cell array, á la msdanalyzer's tracks,
    %    size: {nBeads, nReps}
    %   Apply processing to values and return to same property,
    %    updating the state properties; store state of data in list of
    %    properties, accessed through a display function
    %  Setup and configuration
    %   Write a constructor
    %   Write a configurator public method for each process and plot
    %    Processes need 
    %   Add private properties and state variable for each process
    %  Loading data
    %   dat and csv files. Configurable multi-D csv handling.
    %  Loading images
    %  Processing data
    %  Plotting data
    %  Tracking if stored data is up to date with configuration
    properties (Constant)
        % Constants go in here
        c = 3e8;
    end
    
    properties (SetAccess = private)
        % Private properties go in here
        dSize = 'big';
        
        % Example properties for msd calculation
        msdOpts = struct('upToDate', false);
        msdOpts.doNorm = false;
    end
    
    properties (SetAccess = private, GetAccess = private, Hidden = true)
        % Secret properties go in here
        password = hunter42;
        
    end
    
    
    %% Constructor
    methods
        
        function obj = beadprocessor(args)
            % et voila! You have a value
        end
    end
    
    %% Private methods
    
    methods (Access = private)
        % Call these methods internally with obj.methodName
        x = methodName(obj)
        
    end
    
    %% Static methods
    
    methods (Static, Access = private)
        % Short private functions locally
        
        function y = plusOne(z)
            y = z + 1;
        end
        

    end
    
    
end

