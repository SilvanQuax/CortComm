%% Generate figures from paper

clear all
close all

%% Print panels figure 2

load('Coherence_abs_data.mat')

settings.Figure=2;
settings.firepattern=1;
settings.sth=1;
settings.pow=1;
settings.Gamma=1;

parameters.p1_t=5;
parameters.p2_t=1;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

load('Gamma_EI_data.mat')

settings.Figure=2;
settings.firepattern=1;
settings.sth=1;
settings.pow=1;
settings.EI_Gamma=1;

parameters.p1_t=5;
parameters.p2_t=1;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 3

load('Coherence_abs_data.mat')

settings.Figure=3;
settings.firepattern=1;
settings.sth=1;
settings.pow=1;

parameters.p1_t=5;
parameters.p2_t=4;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 4

load('Coherence_abs_data.mat')

settings.Figure=4;
settings.firepattern=1;
settings.sth=1;
settings.frate=1;
settings.pow=1;
settings.coherence=1;

parameters.p1_t=5;
parameters.p2_t=4;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 5

load('V1V2_data_restored.mat')

settings.Figure=5;
settings.firepattern=1;
settings.sth=1;
settings.frate=1;
settings.pow=1;
settings.coherence=0;
settings.V1V2=1;

parameters.p1_t=5;
parameters.p2_t=5;
parameters.tr_t=1;
data_A=[];
% [data_A]=analysis_paper(parameters,settings,firings_all);
   
[data_A]=wavelet_analysis_paper(settings,parameters,firings_all,data_A);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 6

load('Stim1_data2.mat')

settings.Figure=6;
settings.firepattern=1;
settings.sth=1;
settings.frate=1;
settings.pow=0;
settings.coherence=0;
settings.cohspecgram=0;
settings.pearson=0;
settings.granger=0;
settings.grangergram=0; 
settings.crossgram=0;
settings.mutual_info=0;
settings.stim_phase=0;
settings.phase_diff=0;
settings.V1V2=0;


parameters.p1_t=7;
parameters.p2_t=1;
parameters.tr_t=3;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 7

load('MI_data2.mat')
settings.Figure=7;
settings.firepattern=1;
settings.sth=1;
settings.frate=1;
settings.mutual_info=1;
settings.MI=1;

parameters.p1_t=1;
parameters.p2_t=1;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

load('MI_mul_stim_data2.mat')
settings.Figure=7;
settings.firepattern=1;
settings.sth=1;
settings.frate=1;
settings.mutual_info=1;
settings.MI=2;

parameters.p1_t=1;
parameters.p2_t=1;
parameters.tr_t=1;

[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

%% Print panels figure 8

load('Grang2_data.mat')
settings.Figure=8;
settings.firepattern=1;
settings.sth=1;
settings.granger=1;

parameters.p1_t=1;
parameters.p2_t=1;
parameters.tr_t=1;

fig_generator_paper(settings,parameters,data_A)

settings.granger=0;
settings.coherence=1;
[data_A]=analysis_paper(parameters,settings,firings_all);

fig_generator_paper(settings,parameters,data_A)

