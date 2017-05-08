% Test Example
%
% Example 1 page 16
%
%  Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051

dt=0.001;

t=0:0.001:1;

x=(2*(t-0.5).^2+0.2).*sin(20*pi*t+0.2*cos(40*pi*t));

y=4*(t-0.5).^2;

z=x+y+1;

plot_imf_v8([x;y+ones(1,length(t))])

figure
plot(t,z,'b','linewidth',2)
set(gca,'fontsize', 20);


%%

opts=Settings_IF('IF.delta',10^-3,'IF.NIMFs',1,'plots',0,'IF.Xi',2,'IF.extensionType','p','IF.alpha','ave');

[IMF_2,logM] = IF_v6(z,opts);

plot_imf_v8(IMF_2)

plot_imf_v7(IMF_2,[x;y+ones(1,length(t))],t,'IMFs','Ground truth')



  