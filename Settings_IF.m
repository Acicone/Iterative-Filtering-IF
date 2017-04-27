function options = Settings_IF(varargin)
% 
% Settings_IF Constructs option structure for the algorithms
%
% EXAMPLES
%
% OPTIONS = Settings_IF with no input arguments returns
%   setting structure with default values
%
% OPTIONS = Settings_IF('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = Settings_IF(OLDOPTIONS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTIONS.
%
%
% Settings_IF PROPERTIES :
%
% GENERAL
% 
%   saveEnd          - (0) If value is 1 save workspace at the end
%
%   saveInter        - (0) If value is 1 save intermediate results every
%                      time a new IMF is generated.
%
%   verbose          - (1) 0 the method is silent, 1 normal, >1 loud, 
%
%   maxTime          - (Inf) Approximate max time in seconds one allows the
%                      algorithm to run. At the end of each iteration the  
%                      algorithm checks the elapsed time, if it has run for
%                      more than maxTime or maxTime - (time of last iteration)
%                      it stops iterating.
%                      Termination due to reaching of maxTime is signaled
%                      with a msg and in info.status. 
%       
%   plots            - (0) the algorithm does not produce plots,
%                      1 it produces them 
%   
%   saveplots        - (0) Set to 1 to save automatically each plot
%
%
% SPECIFIC PARAMETERS
%
%
%  IF.delta       (0.001) Stopping criterion
%  IF.ExtPoints   (3)     Number of extrema allowed in the remainder
%  IF.NIMFs       (1)     Number of IMFs we want to produce, not counting
%                           the remainder
%  IF.Xi          (1.6)   Parameter we use to tune the mask length
%  IF.extensionType 
%                  ('c') - constant
%                   'p'  - periodical
%                   'r'  - reflection
%  IF.alpha       ('ave') Parameter used for the mask length computation.
%                           Allowed values [0,1], 'ave' or 'Almost_min'.
%                           If set to 0 the mask length is proportinal to the 
%                           minimum distance between two subsequent extrema.
%                           If set to 1 then it is proportional to the 
%                           maximum distance. If set to 'ave' the mask length 
%                           equals round(2*Xi*(length of the signal)/
%                           (number of extrema)). Finally if set to
%                           'Almost_min' it is set to a value close to the
%                           minimum.
%  IF.MaxInner    (200)   Maximum number of inner steps allowed.
%
% ------------------------------------------------------
% EXAMPLE
%          
%   >> options = Settings_IF('IF.delta',0.08,'IF.NIMFs',5,'plots',1) 
%   >> IMF = IF_v2(x,options)
%              
%  Executes algorithm IF_v2 with delta = 0.08, stops after we have at most 5 IMFs and a trend, it produces plots.                            
% ------------------------------------------------------      
%
% See also IF_v2
%
% Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%      Signal Decomposition and Instantaneous Frequency analysis'. 
%      http://arxiv.org/abs/1411.6051
%

% (Ripped from sdpsettings.m by Johan Lufberg)


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help Settings_IF
    return;
end


Names = {
    % General   
    'saveEnd'
    'saveInter'
    'verbose'    
    'maxTime'
    'plots'
    'saveplots'   
             
    % IF
    'IF.delta'
    'IF.ExtPoints'
    'IF.NIMFs'
    'IF.Xi'
    'IF.extensionType'
    'IF.alpha'
    'IF.MaxInner'
};

obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    
    options.saveEnd = 0;
    options.saveInter = 0;
    options.verbose = 1;    
    options.maxTime = Inf;
    options.plots = 0.0;
    options.saveplots = 0;     
        
    % IF
    options.IF.delta = 0.001;
    options.IF.ExtPoints=3;
    options.IF.NIMFs=1;
    options.IF.Xi=1.6;
    options.IF.extensionType='c';
    options.IF.alpha='ave';
    options.IF.MaxInner=200;
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;       % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;    % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end

