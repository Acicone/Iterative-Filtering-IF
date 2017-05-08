function [IMF,logM] = IF_v6(f,options,M)


%
%  function IMF = IF_v6(f,options)
%
% It generates the decomposition of the signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(K, :)
%
% where the last row in the matrix IMF is the trend and the other rows
% are actual IMFs
%
%                                Inputs
%
%   f         Signal to be decomposed
%
%   options    Structure, generated using function decompSettings_IF, containing
%              all the parameters needed in the various algorithms
%
%   M         Mask length values for each Inner Loop
%
%                               Output
%
%   IMF       Matrices containg in row i the i-th IMF. The last row
%              contains the remainder-trend.
%
%   logM      Mask length values used for each IMF
%
%   See also SETTINGS_IF, GETMASK_V1, MAXMINS_v3, PLOT_IMF_V8.
%
%  Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051




%% deal with the input

if nargin < 1,  help IF_v6; return; end
if nargin < 2, options = Settings_IF; end
if nargin < 3, M = []; end

FigCol = 'ckmygr'; % Plot Colors
logM=zeros(1,options.IF.NIMFs);

N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end
if size(f,1)>1
    disp('Wrong dataset, the signal must be a single row vector')
end
IMF =zeros(options.IF.NIMFs,N);

nameFile=sprintf('%1.0d',sum(round(clock*1000)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Iterative Filtering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
load('prefixed_double_filter','MM');



% k = length(maxmins);
% diffMaxmins=diff(maxmins);
%% Create a signal without zero regions and compute the number of extrema
f_pp=f;
f_pp(abs(f)<=10^-18)=[];
maxmins_pp=Maxmins_v3(f_pp,options.IF.extensionType);
diffMaxmins_pp=diff(maxmins_pp);
N_pp=length(f_pp);
k_pp = length(maxmins_pp);

if options.plots>=1
    maxmins=Maxmins_v3(f,options.IF.extensionType);
    figMask=figure;
    title(['Maxima for IMF_{1}'])
    plot(f,'b','LineWidth',2)
    hold on
    plot(maxmins,f(maxmins),'kx','markersize',10,'linewidth',2)
    legend('Signal','Maxes and mins')
    set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    hold off
    
    figMask2=figure;
    title(['Mask length IMF_{1}'])
    plot(maxmins(1:end-1),diff(maxmins),'kx','markersize',10,'linewidth',2)
    hold on
    plot(spline(maxmins(1:end-1),diff(maxmins),1:N),'linewidth',2)
    set(figMask2,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
end



countIMFs=0;

while countIMFs < options.IF.NIMFs && k_pp>=options.IF.ExtPoints && toc<options.maxTime
    countIMFs=countIMFs+1;
    
    SD=1;
    h=f;
    
    if isempty(M) || length(M)<countIMFs
        
        if isa(options.IF.alpha,'char')
            if strcmp(options.IF.alpha,'ave') % Using an average mask length
                m = 2*round(N_pp/k_pp*options.IF.Xi);
            elseif strcmp(options.IF.alpha,'Almost_min') % Using an almost min mask length
                if 2*round(options.IF.Xi*prctile(diffMaxmins_pp,30))<2*round(N_pp/k_pp*options.IF.Xi)
                    m = 2*round(options.IF.Xi*prctile(diffMaxmins_pp,30));
                else
                    m = 2*round(N_pp/k_pp*options.IF.Xi);
                end
                if countIMFs>1
                    if m<=logM(countIMFs-1)
                        m=ceil(logM(countIMFs-1)*1.1);
                    end
                end
                
            else
                disp(' Value of alpha not recognized')
                return
            end
        else % using a fixed value alpha
            m = 2*round(options.IF.Xi*(max(diffMaxmins_pp)*options.IF.alpha+min(diffMaxmins_pp)*(1-options.IF.alpha)));
        end
    else
        m=M(countIMFs);
    end
    
    inStepN=0;
    if options.verbose>0
        fprintf('\n IMF # %1.0d   -   # Extreme points %5.0d\n',countIMFs,k_pp)
        fprintf('\n  step #            SD             Mask length        %% Negative Maxima     %% Positive Minima \n\n')
    end
    logM(countIMFs)=m;
    a = get_mask_v1(MM,m);
    if options.plots>0 %&& rem(inStepN,5)==0
        if gcf > 30
            close all
        end
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    while SD>options.IF.delta && inStepN < options.IF.MaxInner
        inStepN=inStepN+1;
        
        %maxmins = Maxmins_v3(f,options.IF.extensionType);
        %k = length(maxmins);
        %diffMaxmins=diff(maxmins);
        
        
        %             if m > N
        %                 ff =[h];
        %                 for i=1:floor(m/N)
        %                     ff=[h ff h];
        %                 end
        %                 len=rem(m,N);
        %                 ff=[h(N-len+1:N) ff h(1:len)];
        %             else
        %                 ff=[h(N-m+1:N) h h(1:m)];
        %             end
        if strcmp(options.IF.extensionType,'p')
            if m > N
                ff = [h(N-rem(m,N)+1:N), h, h, h, h(1:rem(m,N))];
            else
                ff = [h(N-m+1:N), h, h(1:m)]; %this could also just be ff = [h, h, h]
            end
        elseif strcmp(options.IF.extensionType,'r')
            if m > N
                ff = [h(N-rem(m,N)+1:N), fliplr(h), h, fliplr(h), h(1:rem(m,N))];
            else
                ff = [h(m:-1:1), h, h(N:-1:N-m+1)];
            end
        elseif strcmp(options.IF.extensionType,'c')
            ff = [h(1)*ones(1,m), h, h(end)*ones(1,m)];
        end
        
        
        
        h_ave=zeros(1,N);
        for i=1:N
            h_ave(i)=a*transpose(ff(i:i+2*m));
        end
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        SD=norm(h_ave)^2/norm(h)^2;
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        h=h-h_ave;
        
        if options.verbose>0 || options.plots>=1
            [Maxs,Mins]=Maxmins_v3(h,options.IF.extensionType);
            if not(isempty(Maxs))
                NnegMax=sum(h(Maxs)<0);
                NposMin=sum(h(Mins)>0);
                if options.verbose>0
                    
                    
                    fprintf('    %2.0d      %1.14f          %2.0d                   %2.4f               %2.4f\n',inStepN,SD,m,NnegMax/length(Maxs)*100,NposMin/length(Mins)*100)
                end
                if options.plots>=1
                    hMaxs=h(Maxs);
                    posNegMax=find(hMaxs<0);
                    hMins=h(Mins);
                    posPosMin=find(hMins>0);
                    figure(figMask2);
                    title(['Maxima for ' sprintf('IMF%1.0d  Step # %5.0d',countIMFs,inStepN)])
                    plot(h,'b','LineWidth',2)
                    hold on
                    plot(Maxs(posNegMax),hMaxs(posNegMax),'rx','markersize',10,'linewidth',2)
                    plot(Mins(posPosMin),hMins(posPosMin),'yx','markersize',10,'linewidth',2)
                    legend('h','Negative Maxes','Positive Mins')
                    set(figMask2,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                    hold off
                end
            else
                fprintf('    %2.0d      %1.14f          %2.0d\n',inStepN,SD,m)
            end
        end
        
        if options.plots>0 %&& rem(inStepN,5)==0
            if countIMFs>1
                plot_imf_v8([IMF(countIMFs-1,:);h;f-h],1:length(h),sprintf('IMF%1.0d  Step # %5.0d',countIMFs,inStepN),[],figN);
            else
                plot_imf_v8([h;f-h],1:length(h),sprintf('IMF%1.0d  Step # %5.0d',countIMFs,inStepN),[],figN);
                %plot_imf_v8([h;f-h],1:length(h),sprintf('IMF%1.0d  Step # %5.0d',countIMFs,inStepN),[],figN);
            end
            if options.saveplots>0
                saveas(figN,[nameFile '_IMF' num2str(countIMFs) '_fig_' num2str(inStepN)], 'fig')
                saveas(figN,[nameFile '_IMF' num2str(countIMFs) '_fig_' num2str(inStepN)], 'epsc')
                saveas(figN,[nameFile '_IMF' num2str(countIMFs) '_fig_' num2str(inStepN)], 'png')
            end
            pause(0.01)
        end
    end
    if inStepN >= options.IF.MaxInner
        disp('Max # of inner steps reached')
        %return
    end
    
    IMF(countIMFs,:) = h;
    f=f-h;
    
    %% Create a signal without zero regions and compute the number of extrema
        f_pp=f;
        f_pp(abs(f)<=10^-18)=[];
        maxmins_pp=Maxmins_v3(f_pp,options.IF.extensionType);
        diffMaxmins_pp=diff(maxmins_pp);
        N_pp=length(f_pp);
        k_pp = length(maxmins_pp);
        
    if options.plots>=1
%         figure(figMask);
%         title(['Mask length IMF_{' num2str(countIMFs) '}'])
%         set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%         plot(f ,'b','LineWidth',2)
%         hold on
%         plot(maxmins,f(maxmins),'kx','markersize',10,'linewidth',2)
%         legend('Signal','Maxes and mins')
%         set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        
        hMaxs=h(Maxs);
        posNegMax=find(hMaxs<0);
        hMins=h(Mins);
        posPosMin=find(hMins>0);
        figure(figMask2);
        title(['Maxima for ' sprintf('IMF%1.0d  Step # %5.0d',countIMFs,inStepN)])
        plot(h,'b','LineWidth',2)
        hold on
        plot(Maxs(posNegMax),hMaxs(posNegMax),'rx','markersize',10,'linewidth',2)
        plot(Mins(posPosMin),hMins(posPosMin),'yx','markersize',10,'linewidth',2)
        legend('h','Negative Maxes','Positive Mins')
        set(figMask2,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        hold off
    end
    if options.saveInter==1
        save([nameFile '_intermediate_IF_v6.mat'],'IMF','f','logM','-v7.3');
    end
end %end of while

IMF = [IMF(1:countIMFs,:); f];
logM = logM(1:countIMFs);


if options.plots>=1
    if gcf > 30
        close all
    end
    figN=plot_imf_v8(IMF,1:N);
    for i=1:length(figN)
        set(figN(i),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        if options.saveplots>0
            saveas(figN(i),[nameFile '_IMFs'], 'fig')
            saveas(figN(i),[nameFile '_IMFs'], 'epsc')
            saveas(figN(i),[nameFile '_IMFs'], 'png')
        end
        
    end
end

if options.saveEnd == 1
    save([ 'Final_' nameFile '_IF_v6.mat'],'IMF','logM','-v7.3');
end

end


%% Auxiliar functions


function a=get_mask_v1(y,k)
%
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>10^-14
        fprintf('\n\n Warning!\n\n')
        fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
        fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        a=a/norm(a,1);
    end
end

end


function varargout = Maxmins_v3(f,extensionType)

% Based on Version 2. 
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

if nargin == 1, extensionType = 'c'; end
N = length(f);
Maxs=zeros(1,N);
Mins=zeros(1,N);
df = diff(f);


h = 1;
cIn=0;
if strcmp(extensionType,'p') && df(1) == 0 && df(end) == 0
    while df(h)==0
        cIn=cIn+1;
        h=h+1;
    end
    if df(h) < 0
        Initial_df=-1;
    else
        Initial_df=+1;
    end
end

cmaxs=0;
cmins=0;

if strcmp(extensionType,'c') && df(1) == 0
    while df(h)==0
        h=h+1;
    end
    if df(h) < 0
        cmaxs=cmaxs+1;
        Maxs(cmaxs)=h;
    else
        cmins=cmins+1;
        Mins(cmins)=h;
    end
end

c = 0;

last_df=[];
for i=h:N-2
    if   df(i)*df(i+1) == 0
        if df(i) < 0
            last_df=-1;
            posc = i;
        elseif df(i) > 0
            last_df=+1;
            posc = i;
        end
        c = c + 1;
        if df(i+1) < 0
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=posc+floor((c-1)/2)+1;                
            end
            c=0;
        end
        if df(i+1) > 0
            if last_df==-1
                cmins=cmins+1;
                Mins(cmins)=posc+floor((c-1)/2)+1;                
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1) < 0
        if df(i) < df(i+1)
            cmins=cmins+1;
            Mins(cmins)=i+1;
            last_df=-1;
        else
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=i+1;
            last_df=+1;
        end
    end
end
if c > 0
    if strcmp(extensionType,'p')
        disp('Code to be completed!')
        if Initial_df < 0
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);                
            end            
        elseif Initial_df > 0
            if last_df==+1
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);                 
            end   
        else
            disp('Code missing!')
        end        
    end        
    if strcmp(extensionType,'c')        
        if last_df > 0
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=posc+1;
        else
            cmins=cmins+1;
            Mins(cmins)=posc+1;
        end
    end 
    if Mins(cmins)==0
        Mins(cmins)=N;
    end
    if Maxs(cmaxs)==0
        Maxs(cmaxs)=N;
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
maxmins=sort([Maxs Mins]);

if strcmp(extensionType,'p') % we deal with a periodic signal
    disp('Code to be completed')
    if isempty(maxmins)
        maxmins = 1;
    else
        if maxmins(1)~=1 && maxmins(end)~=N
            if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
                maxmins=[1 maxmins];
            end
        end
    end
elseif strcmp(extensionType,'c')
    if not(isempty(maxmins))
        if maxmins(1) ~= 1 && maxmins(end) ~= N && df(1)~=0 && df(end)~=0
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
            maxmins = [1, maxmins, N];
        elseif maxmins(1) ~= 1 && df(1)~=0
            maxmins = [1, maxmins];
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
        elseif  maxmins(end) ~= N && df(end)~=0
            maxmins = [maxmins, N];
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
        end
    end
elseif strcmp(extensionType,'r')
    disp('Code to be completed')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end