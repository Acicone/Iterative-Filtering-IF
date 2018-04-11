% Test Example
%
% Example 1 page 16
%
%  Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051
%
%  A. Cicone. 'Nonstationary signal decomposition for dummies'. 
%  To appear in the book Advances in Mechanics and Mathematics.
%  ArXiv https://arxiv.org/abs/1710.04844
%
%  A. Cicone, H. Zhou. 'Iterative Filtering algorithm numerical analysis 
%  with new efficient implementations based on FFT'
%  ArXiv http://arxiv.org/abs/1802.01359

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

[IMF0,logM] = IF_v8_2(z);

plot_imf_v8(IMF0)

%%

opts=Settings_IF_v1('IF.Xi',2,'IF.alpha','ave');

[IMF,logM] = IF_v8_2(z,opts);

plot_imf_v8(IMF)

%%
plot_imf_v7(IMF,[x;y+ones(1,length(t))],t,'IMFs','Ground truth')



  