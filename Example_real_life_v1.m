% Test Example
%
% Example 7 page 25 - Length of the day dataset
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
%
% dataset obtained from http://hpiers.obspm.fr/eoppc/eop/eopc04/eopc04.62-now

load LengthOftheDay_LOD_ALIF_paper

figure
plot(x,'b')
set(gca,'fontsize', 20);

%%

opts=Settings_IF_v1('IF.delta',10^-2,'IF.Xi',3);

[IMF_1,logM] = IF_v8_2(x,opts);

plot_imf_v8(IMF_1)


%%

opts=Settings_IF('IF.delta',10^-2,'IF.NIMFs',100,'plots',0,'IF.Xi',3,'IF.extensionType','c','IF.alpha','ave');

[IMF_2,logM] = IF_v6(x,opts);

plot_imf_v8(IMF_2)



  