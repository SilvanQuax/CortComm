function [data_A] = analysis_new(parameters,settings,firings_all)
% Analysis of generated data
addpath(genpath('D:\U_schijf_master\chronux'))
v2struct(parameters) % Convert parameters to variables
v2struct(settings) % Convert settings to variables
data_A=[];
stimwin=30;
%% To select certain range of parameter values (not used)
if exist('parrange1') 
    ch_par1=ch_par1(parrange1)';
    if isrow(ch_par1)
        ch_par1=ch_par1';
    end
    firings_all=firings_all(parrange1,:,:);
end
if exist('parrange2')
    ch_par2=ch_par2(parrange2);
    if isrow(ch_par2)
        ch_par2=ch_par2';
    end
    firings_all=firings_all(:,parrange2,:);
end
if exist('trrange')
    firings_all=firings_all(:,:,trrange);
end
%% Determine parameter range
szpar1=size(firings_all,1);
szpar2=size(firings_all,2);
sztr=size(firings_all,3);
% Function for plotting SEM
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

%% Fire patterns
if firepattern==1
    szpar1=size(firings_all,1)
    szpar2=size(firings_all,2);
    sztr=size(firings_all,3);
    for ii=1:szpar1
        for iii=1:szpar2
            for tr1=1:sztr
                firings=firings_all{ii,iii,tr1};
                if npop>=1
                firNe11=find(firings(:,2)<=Ne);
                firNi111=find(firings(:,2)> npop*Ne & firings(:,2)<=npop*Ne+Ni1);
                firNi211=find(firings(:,2)> npop*Ne+npop*Ni1 & firings(:,2)<=npop*Ne+npop*Ni1+Ni2);

                firNe1=[firings(firNe11,1),firings(firNe11,2)];
                firNi11=[firings(firNi111,1),firings(firNi111,2)];      
                firNi21=[firings(firNi211,1),firings(firNi211,2)];

                yNe1tr(:,tr1,ii,iii) =histc(firings(firNe11,1),p_in:rt);
                yNi1tr(:,tr1,ii,iii) =histc(firings(firNi111,1),p_in:rt)+histc(firings(firNi211,1),p_in:rt);
                if stim==1
                    if length(stimonset)==szpar1
                    amp1(tr1,ii,iii)=max(yNe1tr(round(stimonset(ii)):round(stimonset(ii))+stimwin,tr1,ii,iii));
                    frateNe1st(tr1,ii,iii)=length(find(firings(:,2)<=1*Ne & firings(:,1)>=stimonset(ii) & firings(:,1)<=stimonset(ii)+stimwin))/(Ne)*1000/stimwin;
                    else
                    frateNe1st(tr1,ii,iii)=length(find(firings(:,2)<=1*Ne & firings(:,1)>=stimonset & firings(:,1)<=stimonset+stimwin))/(Ne)*1000/stimwin;
                    amp1(tr1,ii,iii)=max(yNe1tr(round(stimonset):round(stimonset)+stimwin,tr1,ii,iii));
                    end
                end
                frateNe1(ii,iii,tr1)=length(find(firings(:,2)<=1*Ne))/(Ne)*1000/rt;
                frateNi1(ii,iii,tr1)=length(find(firings(:,2)<=Ne+Ni1 & firings(:,2)>Ne))/(Ni1)*1000/rt;
                frateNi2(ii,iii,tr1)=length(find(firings(:,2)<=Ne+Ni1+Ni2& firings(:,2)>Ne+Ni1))/(Ni2)*1000/rt;
                end

                if npop>=2
                stimwin2=20;
                firNe21=find(firings(:,2)> Ne & firings(:,2)<=2*Ne);%2*Ne
                firNi121=find(firings(:,2)> npop*Ne+Ni1 & firings(:,2)<=npop*Ne+2*Ni1);
                firNi221=find(firings(:,2)> npop*Ne+npop*Ni1+Ni2 & firings(:,2)<=npop*Ne+npop*Ni1+2*Ni2);

                firNe2=[firings(firNe21,1),firings(firNe21,2)];
                firNi12=[firings(firNi121,1),firings(firNi121,2)];
                firNi22=[firings(firNi221,1),firings(firNi221,2)];

                yNe2tr(:,tr1,ii,iii) =histc(firings(firNe21,1),p_in:rt);
                yNi2tr(:,tr1,ii,iii) =histc(firings(firNi121,1),p_in:rt)+histc(firings(firNi221,1),p_in:rt);
                if stim==1
                    if length(stimonset)==szpar1
                    amp2(tr1,ii,iii)=max(yNe2tr(round(stimonset(ii)):round(stimonset(ii))+stimwin2,tr1,ii,iii));
                    frateNe2st(tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset(ii) & firings(:,1)<=stimonset(ii)+stimwin2))/(Ne)*1000/stimwin2;
                    else
                    amp2(tr1,ii,iii)=max(yNe2tr(round(stimonset):round(stimonset)+stimwin2,tr1,ii,iii));
                    frateNe2st(tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset & firings(:,1)<=stimonset+stimwin2))/(Ne)*1000/stimwin2;
                    end  
                end
                frateNe1(ii,iii,tr1)=length(find(firings(:,2)<=1*Ne))/(Ne)*1000/rt;
                frateNe2(ii,iii,tr1)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne))/(Ne)*1000/rt;
                frateNi11(ii,iii,tr1)=length(find(firings(:,2)<=npop*Ne+Ni1 & firings(:,2)>npop*Ne))/(Ni1)*1000/rt;
                frateNi12(ii,iii,tr1)=length(find(firings(:,2)<=npop*Ne+npop*Ni1 & firings(:,2)>npop*Ne+Ni1))/(Ni1)*1000/rt;
                frateNi1(ii,iii,tr1)=length(find(firings(:,2)>=npop*Ne & firings(:,2)<=npop*Ne+Ni1 | firings(:,2)>=npop*Ne+npop*Ni1 & firings(:,2)<=npop*Ne+npop*Ni1+Ni2))/((Ni1+Ni2))*1000/rt;
                frateNi2(ii,iii,tr1)=length(find(firings(:,2)>=npop*Ne+Ni1 & firings(:,2)<=npop*Ne+npop*Ni1 | firings(:,2)>=npop*Ne+npop*Ni1+Ni2 & firings(:,2)<=npop*Ne+npop*Ni1+npop*Ni2))/((Ni1+Ni2))*1000/rt;

                % Analysed data to struct
                end

                if npop>=3
                firNe31=find(firings(:,2)> 2*Ne & firings(:,2)<=3*Ne);
                firNi131=find(firings(:,2)> npop*Ne+2*Ni1 & firings(:,2)<=npop*Ne+3*Ni1);
                firNi231=find(firings(:,2)> npop*Ne+npop*Ni1+2*Ni2 & firings(:,2)<=npop*Ne+npop*Ni1+3*Ni2);

                firNe3=[firings(firNe21,1),firings(firNe21,2)];
                firNi13=[firings(firNi121,1),firings(firNi121,2)];
                firNi23=[firings(firNi221,1),firings(firNi221,2)];

                yNe3tr(:,tr1,ii,iii) =histc(firings(firNe31,1),p_in:rt);
                yNi3tr(:,tr1,ii,iii) =histc(firings(firNi131,1),p_in:rt)+histc(firings(firNi231,1),p_in:rt);
                end


                
               
                if ii==p1_t && iii==p2_t && tr1==tr_t
                firNe1t=firNe1;
                firNi11t=firNi11;
                firNi21t=firNi21;
                if npop>=2;
                firNe2t=firNe2;
                firNi12t=firNi12;
                firNi22t=firNi22;
                end
               
                end
            end
        end
    end
    data_A=v2struct(yNe1tr,yNi1tr,firNe1t,firNi11t,firNi21t,frateNe1,frateNi1,frateNi2);
    
    if npop>=2
     data_A=v2struct(yNe1tr,yNi1tr,firNe1t,firNi11t,firNi21t,frateNe1,frateNi1,frateNi2,yNe2tr,yNi2tr,firNe2t,firNi12t,firNi22t,frateNe2,frateNi12);
%    data_A=v2struct(yNe1tr,yNi1tr,firNe1t,firNi11t,firNi21t,yNe2tr,yNi2tr,firNe2t,firNi12t,firNi22t);

    end
    
    if stim==1
        data_A.frateNe1st=frateNe1st;
        if npop==2
            data_A.frateNe2st=frateNe2st;
        end
    
    end
    
end

%% Power
if pow==1
        for ii=1:szpar1
            for iii=1:szpar2

                params.Fs=1000;
                params.fpass=5:100;
                params.tapers=[5 10];
%                 params.tapers=[1 2 1];

                params.trialave=1;
                params.err=[2,0.05];

                [pxx f Serr]=mtspectrumc(yNe1tr(:,:,ii,iii),params);
                pxx1(:,:,ii,iii)=pxx;
                f1(:,:,ii,iii)=f';
                S1err(:,:,ii,iii)=Serr;
                lb_gam=find(f1(:,:,ii,iii)>20,1);
                pxx_sm=smooth(pxx1(:,:,ii,iii),5,'moving');
                [P_gam,I]=max(pxx1(lb_gam:end,:,ii,iii));

                F_gam=f1(I+lb_gam-1,:,ii,iii);
                P_gam_all(ii,iii)=P_gam;
                F_gam_all(ii,iii)=F_gam;
                lb_a=find(f1>7,1);
                ub_a=find(f1<13,1,'last');
                P_alp_all(ii,iii)=mean(pxx1(lb_a:ub_a));

                if npop>=2

                    [pxx2 f2]=mtspectrumc(yNe2tr(:,:,ii,iii),params);
                end
                if npop>=3
                    yNe3(:,:)=yNe3tr(:,1,ii,iii);
                    [pxx3 f3]=mtspectrumc(yNe3(100:rt),params);
                end

            end
        end
data_A.pxx1=pxx1;
data_A.f1=f1;
data_A.S1err=S1err;

end
%% Coherence
if coherence==1

    if npop>=2
        for ii=1:szpar1
            for iii=1:szpar2
                params.Fs=1000;
                params.fpass=5:100;
                params.tapers=[5,(rt-p_in)/1000,1];
                params.trialave=1;
                params.err=[2,0.05];


                [C1(:,:,ii,iii),phi1(:,:,ii,iii),S121,S11,S21,fc1(:,:,ii,iii),~,~,Cerr(:,:,:,ii,iii)]=coherencyc(yNe1tr(:,:,ii,iii),yNe2tr(:,:,ii,iii),params);
                [C2(:,:,ii,iii),phi2,S122,S12,S22,fc2(:,:,ii,iii),~,~,Cerr2(:,:,:,ii,iii)]=coherencyc_phase(yNe1tr(:,:,ii,iii),yNe2tr(:,:,ii,iii),params);
                [C3(:,:,ii,iii),phi3,S123,S13,S23,fc3(:,:,ii,iii),~,~,Cerr3(:,:,:,ii,iii)]=coherencyc_amp(yNe1tr(:,:,ii,iii),yNe2tr(:,:,ii,iii),params);
                [C4(:,:,ii,iii),phi4,S124,S14,S24,fc4(:,:,ii,iii),~,~,Cerr4(:,:,:,ii,iii)]=coherencyc(yNe1tr(:,:,ii,iii),circshift(yNe2tr(:,:,ii,iii),[0,1]),params);
                [C5(:,:,ii,iii),phi5,S125,S15,S25,fc5(:,:,ii,iii),~,~,Cerr5(:,:,:,ii,iii)]=coherencyc_phase(yNe1tr(:,:,ii,iii),circshift(yNe2tr(:,:,ii,iii),[0,1]),params);
                [C6(:,:,ii,iii),phi6,S126,S16,S26,fc6(:,:,ii,iii),~,~,Cerr6(:,:,:,ii,iii)]=coherencyc_amp(yNe1tr(:,:,ii,iii),circshift(yNe2tr(:,:,ii,iii),[0,1]),params);

                Imin=find(fc1(:,:,ii,iii)>30,1);
                Imax=find(fc1(:,:,ii,iii)<50,1,'last');
                
                C_gam=mean(C1(Imin:Imax,:,ii,iii));
                C_gam2=mean(C2(Imin:Imax,:,ii,iii));
                C_gam3=mean(C3(Imin:Imax,:,ii,iii));
                C_gam4=mean(C4(Imin:Imax,:,ii,iii));
                C_gam5=mean(C5(Imin:Imax,:,ii,iii));
                C_gam6=mean(C6(Imin:Imax,:,ii,iii));
                
%                 phigam(ii,iii)=mean(phi1(Imin:Imax,:,ii,iii));

                Imean=round((Imin+Imax)/2);
%                 if iii==[7]
%                     phi1a=phi1(Imean,:,ii,3:7);
%                     phi1a=phi1a(:);
%                     phase_range=[-pi:0.1*pi:pi];
%                     phigam(ii,:)=histc(phi1a,phase_range)
% 
%                     figure(112)
%                     phibin=[-19/20*pi:pi/10:19/20*pi,-19/20*pi];   
%                     phigam(:,end)=phigam(:,1);
%                     stralp={'-0.5\pi/2','-0.4\pi','-0.3\pi','-0.2\pi','0.1\pi','0','0.1\pi','0.2\pi','0.3\pi','0.4\pi','0.5\pi'};
%                     figure(9)
%                     set(gcf,'Position',[300 300 1500 220])
%                     ha = tight_subplot(1,9,[.01 .03],[.1 .01],[.01 .01]);
% 
%                 end
                
                if params.trialave==1
                C_gam_err(ii,iii)=mean(Cerr(2,:,:,ii,iii)-Cerr(1,:,:,ii,iii))/2;
                C_gam_err2(ii,iii)=mean(Cerr2(2,:,:,ii,iii)-Cerr2(1,:,:,ii,iii))/2;
                C_gam_err3(ii,iii)=mean(Cerr3(2,:,:,ii,iii)-Cerr3(1,:,:,ii,iii))/2;
                C_gam_err4(ii,iii)=mean(Cerr4(2,:,:,ii,iii)-Cerr4(1,:,:,ii,iii))/2;
                C_gam_err5(ii,iii)=mean(Cerr5(2,:,:,ii,iii)-Cerr5(1,:,:,ii,iii))/2;
                C_gam_err6(ii,iii)=mean(Cerr6(2,:,:,ii,iii)-Cerr6(1,:,:,ii,iii))/2;
                else
                C_gam_err(ii,iii)=0;
                C_gam_err2(ii,iii)=0;
                C_gam_err3(ii,iii)=0;
                C_gam_err4(ii,iii)=0;
                C_gam_err5(ii,iii)=0;
                C_gam_err6(ii,iii)=0;
                end
                C_gam_all(ii,iii)=mean(C_gam);
                C_gam_all2(ii,iii)=mean(C_gam2);
                C_gam_all3(ii,iii)=mean(C_gam3);
                C_gam_all4(ii,iii)=mean(C_gam4);
                C_gam_all5(ii,iii)=mean(C_gam5);
                C_gam_all6(ii,iii)=mean(C_gam6);

            end
        end
%         for ii=1:9
%             axes(ha(ii));
%             polar(1,1) 
%             hold on
%             h(ii)=polar(phibin,phigam(ii,:)/norm(phigam(ii,:)));
%             strpha=['\Delta\phi \alpha = ' stralp{ii}];
%             xlabel(strpha)
%             title('\gamma Phase')
%             set(h(ii),'LineWidth',2)
%             hold off
%          end
        
        
        
        data_A.C1=C1;
        data_A.fc1=fc1;
        data_A.Cerr=Cerr;
        data_A.C2=C2;
        data_A.fc2=fc2;
        data_A.Cerr2=Cerr2;
        data_A.C3=C3;
        data_A.fc3=fc3;
        data_A.Cerr3=Cerr3;
        data_A.C4=C4;
        data_A.fc4=fc4;
        data_A.Cerr4=Cerr4;
        data_A.C5=C5;
        data_A.fc5=fc5;
        data_A.Cerr5=Cerr5;
        data_A.C6=C6;
        data_A.fc6=fc6;
        data_A.Cerr6=Cerr6;
        
    end
end
%% Spectogram/Cohgram
if cohspecgram==1
    figure(9); clf;
    title('Coh gram')
    % ColorSet = [linspace(0,1,ch_par1)',linspace(1,0,ch_par1)',linspace(0,1,ch_par1)'];
    % set(gca,'ColorOrder',ColorSet);

    window1=[0.3 0.1];
    window2=[0.5 0.1];
    
    params.Fs=1000;
    params.fpass=[25 100];
    params.tapers=[10 window1(1) 1];
    params.trialave=1;

        for ii=1:szpar1

            [C1,phi1,S121,S11,S21,t1,f1]=cohgramc(yNe1tr(100:rt,:,ii,p2_t),yNe2tr(100:rt,:,ii,p2_t),window1,params);
        h(ii)=subplot(2,length(ch_par1),(ii));
        pcolor(t1,f1,C1')
        shading flat
        titlestr=[str_ch_par1,'=',num2str(ch_par1(ii))];
        title(h(ii),titlestr)

        end
    %     colorbar



    params.fpass=[5 30];
    params.tapers=[5 window2(1) 1];
    params.trialave=1;

        for ii=1:szpar1

            [C2,phi2,S121,S11,S21,t,f2]=cohgramc(yNe1tr(100:rt,:,ii,p2_t),yNe2tr(100:rt,:,ii,p2_t),window2,params);
        subplot(2,length(ch_par1),(length(ch_par1)+ii));
        pcolor(t,f2,C2')
        shading flat
        end
    %     colorbar



    figure(10)
        pcolor(t1,f1,phi1')
        shading flat
    %     colorbar



    params.fpass=[25 100];
    params.tapers=[10 window1(1) 1];
    params.trialave=1;

    figure(11)
    title('Spec gram E2')  
        for ii=1:szpar1

            [S11,t,ff1]=mtspecgramc(yNe1tr(100:rt,:,ii,p2_t),window1,params);
            [S12,t,ff1]=mtspecgramc(yNe2tr(100:rt,:,ii,p2_t),window1,params);
            h(ii)=subplot(2,length(ch_par1),ii);
        %     I=find(S11<0.025&C1>0.9);
        %     S11(ind2sub(size(S11),I))=1;

            pcolor(t,ff1,S11')
            shading interp
            titlestr=[str_ch_par1,'=',num2str(ch_par1(ii))];
            title(h(ii),titlestr)
%             colorbar
            end

    params.fpass=[5 30];
    params.tapers=[5 window2(1) 1];
    params.trialave=1;
        for ii=1:szpar1

            [S21,t,ff2]=mtspecgramc(yNe1tr(100:rt,:,ii,p2_t),window2,params);
            [S22,t,ff2]=mtspecgramc(yNe2tr(100:rt,:,ii,p2_t),window2,params);
        subplot(2,length(ch_par1),length(ch_par1)+ii);
        pcolor(t,ff2,S21')
        shading flat
        end
        
end
%% Pearson correlation    
if pearson==1
    ii=3:5;  
    vC1=reshape(C1(:,ii),[numel(C1(:,ii)),1]);
    vphi1=reshape(phi1(:,ii),[numel(phi1(:,ii)),1]);
    vS11=reshape(S11(:,ii),[numel(S11(:,ii)),1]);
    vS12=reshape(S12(:,ii),[numel(S12(:,ii)),1]);


    % scatter(vC1,vS11);
    % lsline
    figure(12)
    distC1=histc(vC1,0:0.05:1);
    distS11=histc(vS11,0:0.02:0.4);
    hold all

    scatter(vC1,vS11);
    plot((0:0.05:1)+0.025,distC1/max(distC1)*max(vS11))
    plot(distS11/max(distS11)*max(vC1),(0:0.02:0.4)+0.01)
    lsline
    xlabel('Coherence')
    ylabel('Gamma power')
    hold

    figure(13)
    distphi1=histc(vphi1,-3.14:0.1:3.14);

    hold all
    II=find(vphi1>-1.5&vphi1<3);

    scatter(vphi1,vC1);
    plot(-3.14:0.1:3.14,distphi1/max(distphi1))
    plot(distC1/max(distC1)*max(2*pi)-pi,(0:0.05:1)+0.025)
    xlabel('Coherence')
    ylabel('Gamma phase difference')
    lsline

    [r1,p1]=corr(vC1,vS11)
    [r2,p2]=corr(vC1,vS12)
    [r12,p12]=corr(vC1,vS12.*vS11)
    hold

    [rphi1,pphi1]=corr(vC1(II),vphi1(II))
end
%% Granger causality
if granger==1
    for p1_t=1:szpar1;
   
    
    ntrials=size(yNe1tr,2);
    nvars=2;
    regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
    icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

    morder    = 'BIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
    momax     = 20;%12     % maximum model order for model order estimation
    nsamps    = 10;     % number of bootstrap samples
    acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)

    tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
    alpha     = 0.05;   % significance level for significance test
    mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
    
    data_red  = 2
    nobs      = length(yNe1tr)/data_red;
    fs        = 1000/data_red;    % sample rate (Hz)
    fres      = [];     % frequency resolution (empty for automatic calculation)


    %%
    X(1,:,:)=yNe1tr(1:data_red:end,:,p1_t,p2_t);
    X(2,:,:)=yNe2tr(1:data_red:end,:,p1_t,p2_t);
% %     X(1,:,:)=reshape(sum(reshape(yNe1tr(101:end,:,phi_test,alpha_test),4,250,tr)),250,tr);
% %     X(2,:,:)=reshape(sum(reshape(yNe2tr(101:end,:,phi_test,alpha_test),4,250,tr)),250,tr);

%     X1=(3*sin((1:(rt-p_in+1))*2*pi*10/1000-pi/2))'; % Activate for conditional GC
%     if exist('A11_all','var')
%         X(3,:,:)=A11_all(p_in:end,:,p1_t,p2_t)+0.1*randn(rt-p_in+1,tr);
%     else
%         X(3,:,:)=repmat(X1,1,tr)+0.5*randn(rt-p_in+1,tr);
%     end
%     size(X)
% %     X=X(1:2,:,:);
    %% Calculate information criteria up to specified maximum model order.

    ptic('\n*** tsdata_to_infocrit\n');
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
    ptoc('*** tsdata_to_infocrit took ');

    % Plot information criteria.

    figure(20); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');



    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);


    % Select model order.


    if strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
    end
    %% Estimate VAR model of selected order from data.

    ptic('\n*** tsdata_to_var... ');
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    ptoc;

    % Check for failed regression

    assert(~isbad(A),'VAR estimation failed');
    %% Estimate autocovariance
    ptic('*** var_to_autocov... ');
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    ptoc;

    % The above routine does a LOT of error checking and issues useful diagnostics.
    % If there are problems with your data (e.g. non-stationarity, colinearity,
    % etc.) there's a good chance it'll show up at this point - and the diagnostics
    % may supply useful information as to what went wrong. It is thus essential to
    % report and check for errors here.

    var_info(info,true); % report results (and bail out on error)

    %% Calculate time-domain pairwise-conditional causalities - this just requires
    % the autocovariance sequence.

    ptic('*** autocov_to_pwcgc... ');
    F = autocov_to_pwcgc(G)
    ptoc;

    % Check for failed GC calculation

    assert(~isbad(F,false),'GC calculation failed');

    % Significance test using theoretical null distribution, adjusting for multiple
    % hypotheses.

    pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat) % take careful note of arguments!
    sig  = significance(pval,alpha,mhtc)

    % Plot time-domain causal graph, p-values and significance.

    figure(21); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])

    % For good measure we calculate Seth's causal density (cd) measure - the mean
    % pairwise-conditional causality. We don't have a theoretical sampling
    % distribution for this.

    cd = mean(F(~isnan(F)));

    fprintf('\ncausal density = %f\n',cd);

    %% Calculate spectral pairwise-conditional causalities at given frequency
    % resolution - again, this only requires the autocovariance sequence.

    ptic('\n*** autocov_to_spwcgc... ');
    % [f,fres] = autocov_to_spwcgc(G,fres);
    [f21,fres]=autocov_to_smvgc(G,1,2);
    [f12]=autocov_to_smvgc(G,2,1,fres);
    ptoc;

    %%
    ptic('\n*** bootstrap_tsdata_to_pwcgc\n');

    % FSAMP = bootstrap_tsdata_to_spwcgc(X,morder,fres,nsamps);
    FSAMP21 = bootstrap_tsdata_to_smvgc(X,1,2,morder,fres,nsamps);
    FSAMP12 = bootstrap_tsdata_to_smvgc(X,2,1,morder,fres,nsamps);
    ptoc('*** bootstrap_tsdata_to_pwcgc took ',[],1);

    % (We should really check for failed bootstrap estimates here.)

    % Bootstrap (empirical) confidence intervals.

    [FSUP21,FSLO21] = empirical_confint(alpha,FSAMP21);
    [FSUP12,FSLO12] = empirical_confint(alpha,FSAMP12);


    %%
    % Check for failed spectral GC calculation
    % f21=shiftdim(f(1,2,:),2)';
    % f12=shiftdim(f(2,1,:),2)';
    % f32=shiftdim(f(2,3,:),2)';
    % f31=shiftdim(f(1,3,:),2)';
    % FSUP21=shiftdim(FSUP(1,2,:),2)';
    % FSUP12=shiftdim(FSUP(2,1,:),2)';
    % FSLO21=shiftdim(FSLO(1,2,:),2)';
    % FSLO12=shiftdim(FSLO(2,1,:),2)';
    assert(~isbad(f21,false),'spectral GC calculation failed');
    assert(~isbad(f12,false),'spectral GC calculation failed');
    % Plot spectral causal graph.
    h=length(f21);
    fres = h-1;
    lam = sfreqs(fres,fs)';
    I40=find(lam>40,1);
    mean(FSAMP12(:,I40))
    std(FSAMP12(:,I40))
    mean(FSAMP21(:,I40))
    std(FSAMP21(:,I40))
    [h,p]=ttest2(FSAMP12(:,I40),FSAMP21(:,I40))
    

    plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

    figure(22); clf;
    % plot_spw(f,fs);
    hold all
    ha(1)=plot(lam,f12,'b');
    ha(2)=plot(lam,f21,'r');
    legend(ha([1,2]),'1->2','2->1')
    hold off
    xlim([0 100])
    xlabel('F (Hz)');ylabel('Granger causality');

    figure(23); clf;
    % plot_spw(f,fs);
    hold all
    plot_variance(lam,FSLO12,FSUP12,[0.4,0.4,1])
    ha1(1)=plot(lam,f12,'b');
    plot_variance(lam,FSLO21,FSUP21,[1,0.4,0.4])
    ha1(2)=plot(lam,f21,'r');
    legend(ha1([1,2]),'1->2','2->1')
    hold off
    xlim([0 100])
    xlabel('F (Hz)');ylabel('Granger causality');
    data_A.lam{p1_t}=lam;
    data_A.FSLO12{p1_t}=FSLO12;
    data_A.FSUP12{p1_t}=FSUP12;
    data_A.FSLO21{p1_t}=FSLO21;
    data_A.FSUP21{p1_t}=FSUP21;
    data_A.f12{p1_t}=f12;
    data_A.f21{p1_t}=f21;
    end
%     save('GC_FFI_2','dat2')
end

%% Grangergram
if grangergram==1
    for p1_t=1:9
    for tr_t=1:tr
    regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
    icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
    nsamps    = 5;
    acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)

    tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
    alpha     = 0.05;   % significance level for significance test
    mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

    fs        = 1000;    % sample rate (Hz)
    fres      = [1910];     % frequency resolution (empty for automatic calculation)
    morder=7;
    window=[0.2 0.1];
    rt
    nsteps=floor((rt-p_in+1-window(1)*fs)/(window(2)*fs))+1
    xvec=0.5*window(1)*fs:fs*window(2):(nsteps-1)*fs*window(2)+0.5*window(1)*fs
    f21_all=[];
    f12_all=[];
%     for tr_t=1:tr
    for ii=1:nsteps
        %%
        startstep=1+(ii-1)*window(2)*fs
        endstep=startstep+window(1)*fs-1
        X(1,:,:)=yNe1tr(startstep:endstep,tr_t,p1_t,p2_t);
        X(2,:,:)=yNe2tr(startstep:endstep,tr_t,p1_t,p2_t);
        % X(1,:,:)=reshape(sum(reshape(yNe1tr(101:end,:,phi_test,alpha_test),4,250,tr)),250,tr);
        % X(2,:,:)=reshape(sum(reshape(yNe2tr(101:end,:,phi_test,alpha_test),4,250,tr)),250,tr);
        X1=(3*sin((startstep:endstep)*2*pi*10/1000-pi/2))';
        size(X1)
        X(3,:,:)=repmat(X1,1,1)+0.5*randn(window(1)*fs,1);
        size(X)

        %% Estimate VAR model of selected order from data.

        [A,SIG,~,SIG_err] = tsdata_to_var(X,morder,regmode);
        SIG_all(ii)=SIG(1,1);
        SIG_err_all(ii)=SIG_err(1,1)
        assert(~isbad(A),'VAR estimation failed');
        %% Estimate autocovariance

        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        var_info(info,true); % report results (and bail out on error)

        %% Calculate spectral pairwise-conditional causalities at given frequency
        % resolution - again, this only requires the autocovariance sequence.

        ptic('\n*** autocov_to_spwcgc... ');
        % [f,fres] = autocov_to_spwcgc(G,fres);
        [f21(:,ii,tr_t)]=autocov_to_smvgc(G,1,2,fres);
 
        [f12(:,ii,tr_t)]=autocov_to_smvgc(G,2,1,fres);
        ptoc;

        %%
        % Check for failed spectral GC calculation
        assert(~isbad(f21,false),'spectral GC calculation failed');
        assert(~isbad(f12,false),'spectral GC calculation failed');

%         ptic('\n*** bootstrap_tsdata_to_pwcgc\n');
% 
%         % FSAMP = bootstrap_tsdata_to_spwcgc(X,morder,fres,nsamps);
%         FSAMP21 = bootstrap_tsdata_to_smvgc(X,1,2,morder,fres,nsamps);
%         FSAMP12 = bootstrap_tsdata_to_smvgc(X,2,1,morder,fres,nsamps);
%         ptoc('*** bootstrap_tsdata_to_pwcgc took ',[],1);
% 
%         % (We should really check for failed bootstrap estimates here.)
%         
%         % Bootstrap (empirical) confidence intervals.
%         h=size(f21,1);
%         fres = h-1;
%         lam = sfreqs(fres,fs)';
%         I40=find(lam>40,1);
%         [FSUP21(:,ii),FSLO21(:,ii)] = empirical_confint(alpha,FSAMP21(:,I40));
%         [FSUP12(:,ii),FSLO12(:,ii)] = empirical_confint(alpha,FSAMP12(:,I40));
    end

    
%     f21_all=[f21_all,f21];
%     f12_all=[f12_all,f12];
%     end
        %%
        % Plot spectral causal graph.
        plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);

        h=size(f21,1);
        fres = h-1;
        lam = sfreqs(fres,fs)';
        I100=find(lam<100,1,'last');
        I30=find(lam>30,1);
        I50=find(lam<50,1,'last');
        I40=find(lam<40,1,'last');
        f12m=mean(f12(I30:I50,:,tr_t))
        f21m=mean(f21(I30:I50,:,tr_t))
        var12(p1_t,tr_t)=std(f12m,[],2);
        var21(p1_t,tr_t)=std(f21m,[],2);
    end
end
        figure(24); clf;
        pcolor(xvec,lam(1:I100),f12(1:I100,:,tr_t))
        shading flat
        colorbar
        figure(25); clf;
        pcolor(xvec,lam(1:I100),f21(1:I100,:,tr_t))
        shading flat
        colorbar
        
        I30=find(lam>30,1);
        I50=find(lam<50,1,'last');
        I40=find(lam<40,1,'last');
        f12m=mean(f12(I30:I50,:,tr_t));
        f21m=mean(f21(I30:I50,:,tr_t));
        
        figure(26); clf;
        hold on
%         plot(xvec,f12(I30,:,tr_t),xvec,f21(I40,:,tr_t))
        plot(xvec,f12m,xvec,f21m)

%         plot_variance(xvec,FSLO12,FSUP12,[0.4,0.4,1])
%         plot_variance(xvec,FSLO21,FSUP21,[1,0.4,0.4])
        hold off
        data_A.f12g=f12;
        data_A.f21g=f21;
        data_A.xvec=xvec;
        data_A.I40=I40;
        data_A.var12=var12;
        data_A.var21=var21;
%         p=polyfit(f12_all(I40,:),f21_all(I40,:),1);
%         x1 = linspace(0,0.8);
%         y1 = polyval(p,x1);
%         figure(27); clf;
%         hold on
%         scatter(f12_all(I40,:),f21_all(I40,:))
%         plot(x1,y1)
%         hold off
        
        figure(28);clf;
%         for tr_t=1:tr
        CorrABt(:,tr_t)=xcorr(f12(I40,:,tr_t),f21(I40,:,tr_t))
%         end
        CorrAB=mean(CorrABt,2);
        plot(CorrAB)
        
        figure(29)
        errorbar(xvec,SIG_all,SIG_err_all)
end
%%
if crossgram==1
%UNTITLED Summary of this function goes here
window=100;
step=10;
x=yNe1tr(101:end,:,p1_t,p2_t);
y=yNe2tr(101:end,:,p1_t,p2_t);
rt=length(x);

size(y)

nsteps=(rt-window)/step;
for ii=7
for i=1:nsteps
    xi=x(i*step:i*step+window,ii);
    yi=y(i*step:i*step+window,ii);
	r(i,:,ii)=xcorr(xi,yi);
end
end
rm=mean(r,3);

windowp=[0.1 0.01];
params.fpass=[30 100];
params.Fs=1000;
params.tapers=[10 windowp(1) 1];
params.trialave=1;

[S1,t,ff1]=mtspecgrampb(y,windowp,params);

figure(29)
subplot(2,1,1)
pcolor(1+window/2:step:(rt-window)+window/2,-window:window,rm')
shading flat
colorbar
title('Time resolved cross-correlation')
subplot(2,1,2)
pcolor(t,ff1,S1')
shading flat
colorbar
title('Gamma amplitude')
end
%%
if mutual_info==1
    fr=1;
    spk_time=0;
    win=30;
    szwin=length(win);
    nbootstrap=100;
%     NN_all=round(logspace(2.3,0,10));
NN_all=20;
if MI==1
    NN_all=round(logspace(2.3,0,10));
end

% stimonset=810:20:890
stimonset

    for NN=1:length(NN_all)

    if stim==1
                szpar1=szpar2;
        szpar2=length(stimonset);
        for ii=1:szpar1
            for iii=1:szwin
                for tr1=1:tr
                    firings=firings_all{ii,1,tr1};
                    if length(stimonset)==szpar1
                    resp_data(:,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset(ii) & firings(:,1)<=stimonset(ii)+win(iii)))/(Ne)*1000/win(iii);
                    else
                    resp_data(:,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset & firings(:,1)<=stimonset+win(iii)))/(Ne)*1000/win(iii);
%                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset & firings(:,1)<=stimonset+win(iii),2),401:800)*1000/(win(iii));

                    end
                end
            end
        end
        elseif stim==2
            if MI==1
            %         tr=200;

        
            for ii=1:szpar1
                for iii=1:szpar2
                    for tr1=1:tr
                    ii
                    iii
                    firings=firings_all{ii,iii,tr1};
%                     if stimonset==0
%                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=200,2),401:800)*1000/(rt-200);
%                     else
%                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset-00 & firings(:,1)<=stimonset-00+win(iii),2),401:800)*1000/(win(iii));
%                     stimonset=stimonset_new(tr1,ii)
                    resp_data1(1,tr1,ii,iii)=length(find(firings(:,2)<=Ne+NN_all(NN) & firings(:,2)>Ne & firings(:,1)>=stimonset(iii) & firings(:,1)<=stimonset(iii)+win))/(Ne)*1000/win;
                    resp_data1(2,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne-200+NN_all(NN) & firings(:,2)>Ne+200 & firings(:,1)>=stimonset(iii) & firings(:,1)<=stimonset(iii)+win))/(Ne)*1000/win;
                
                    end 
                   
                end
            end

            elseif MI==2
            szpar1=szpar2;
            szpar2=length(stimonset);
                for ii=1:szpar1
                    for iii=1:szpar2
                        for tr1=1:tr
                        ii
                        iii
                        firings=firings_all{iii,ii,tr1};
    %                     if stimonset==0
    %                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=200,2),401:800)*1000/(rt-200);
    %                     else
    %                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset-00 & firings(:,1)<=stimonset-00+win(iii),2),401:800)*1000/(win(iii));
    %                     stimonset=stimonset_new(tr1,ii)
                        resp_data1(1,tr1,ii,iii)=length(find(firings(:,2)<=Ne+NN_all(NN) & firings(:,2)>Ne & firings(:,1)>=stimonset(iii) & firings(:,1)<=stimonset(iii)+win))/(Ne)*1000/win;
                        resp_data1(2,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne-200+NN_all(NN) & firings(:,2)>Ne+200 & firings(:,1)>=stimonset(iii) & firings(:,1)<=stimonset(iii)+win))/(Ne)*1000/win;

                        end 

                    end
                end
            end
        elseif stim==3
            szpar1=szpar2;
szpar2=length(stimonset);
        for ii=1:szpar1
            for iii=1:szwin
                for tr1=1:tr
                    firings=firings_all{ii,iii,tr1};
                    if length(stimonset)==szpar1
                    resp_data(:,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset(ii) & firings(:,1)<=stimonset(ii)+win(iii)))/(Ne)*1000/win(iii);
                    else
                    resp_data(:,tr1,ii,iii)=length(find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset_new(tr1,ii) & firings(:,1)<=stimonset_new(tr1,ii)+win(iii)))/(Ne)*1000/win(iii);
%                     resp_data(:,tr1,ii,iii)=histc(firings(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=stimonset & firings(:,1)<=stimonset+win(iii),2),401:800)*1000/(win(iii));

                    end
                end
            end
        end
    end
resp_data=(resp_data1(1,:,:,:)-resp_data1(2,:,:,:))%./max(resp_data1(:,:,:),[],1);
size(resp_data)
resp_data(1,1:20,:)
% mean(resp_data(1,1:2:end,3,6),2)
% mean(resp_data(1,2:2:end,3,6),2)
% mean(resp_data(1,1:2:end,7,6),2)
% mean(resp_data(1,2:2:end,7,6),2)

% for ii=1:szpar1
%     figure(ii); clf;
%     size(stimonset_new(:,ii))
% %     size(resp_data1(1,:,ii))
%     hold on
%     scatter(stimonset_new(1:2:end,ii),resp_data1(1,1:2:end,ii)','r')
%     scatter(stimonset_new(2:2:end,ii),resp_data1(2,2:2:end,ii)','r')
%     scatter(stimonset_new(2:2:end,ii),resp_data1(1,2:2:end,ii)','b')
%     scatter(stimonset_new(1:2:end,ii),resp_data1(2,1:2:end,ii)','b')
%     scatter(stimonset_new(1:2:end,ii),resp_data(1,1:2:end,ii)','g')
%     scatter(stimonset_new(2:2:end,ii),resp_data(1,2:2:end,ii)','k')
% 
%     hold off
% 
% end
%     inD=find(firings(:,2)<=2*Ne & firings(:,2)>Ne & firings(:,1)>=200)
%     firings(inD,2)
%     histc(firings(inD,2),401:800)
size(resp_data)
    Itrain=ceil(tr/2);
    Itest=floor(tr/2);
    Y=zeros(tr,1);
    Iodd=1:2:tr;
    Ieven=2:2:tr;
    Y(Iodd)=1;  
%     Y1=Y(1:Itrain);%Train input
%     Y2=Y(Itrain+1:end);%Test input
%     Ix0=find(Y2==0);
%     Ix1=find(Y2==1);

%     mpop1s=[mean(resp_data(1,Iodd,1,2),2),std(resp_data(1,Iodd,1,2),[],2)]
%     mpop1n=[mean(resp_data(1,Ieven,1,2),2),std(resp_data(1,Ieven,1,2),[],2)]
%     mpop2s=[mean(resp_data(2,Ieven,1,2),2),std(resp_data(2,Ieven,1,2),[],2)]
%     mpop2n=[mean(resp_data(2,Iodd,1,2),2),std(resp_data(2,Iodd,1,2),[],2)]
%     mpop1s2=[mean(resp_data(1,Iodd,2,2),2),std(resp_data(1,Iodd,2,2),[],2)]
%     mpop1n2=[mean(resp_data(1,Ieven,2,2),2),std(resp_data(1,Ieven,2,2),[],2)]
%     mpop2s2=[mean(resp_data(2,Ieven,2,2),2),std(resp_data(2,Ieven,2,2),[],2)]
%     mpop2n2=[mean(resp_data(2,Iodd,2,2),2),std(resp_data(2,Iodd,2,2),[],2)]

    for iii=1:szpar2
        for ii=1:szpar1
            for n=1:nbootstrap
                Iboot=[randi(Itrain/2,1,Itrain/2)*2,randi(Itrain/2,1,Itrain/2)*2-1];

                Y1=Y(Iboot);
                Y2=Y(Itrain+Iboot);
                Ix0=find(Y2==0);
                Ix1=find(Y2==1);
            X1=resp_data(:,Iboot,ii,iii)';%Train response
            X2=resp_data(:,Itrain+Iboot,ii,iii)';%Test response
            svmstruct=svmtrain(X1,Y1);%Train SVM
            Y2t=svmclassify(svmstruct,X2);%Classify test responses
            ClassError1(n)=length(find((Y2-Y2t)~=0))/length(Y2t);%Degree misclassified
            if iii==1 && ii==3 && n==1
            Imissclass=find((Y2-Y2t)~=0);
            Iclass=find((Y2-Y2t)==0);
            Imiss=Iboot(Imissclass)+Itrain
            Icl=Iboot(Iclass)+Itrain
            end
            %Mutual information
            Px(1,1)=-1*sum(Y2-1)/length(Y2);
            Px(1,2)=sum(Y2)/length(Y2);
            Py(1,1)=-1*sum(Y2t-1)/length(Y2t);
            Py(1,2)=sum(Y2t)/length(Y2t);
            Pyx(1,1)=-1*sum(Y2t(Ix0)-1)/length(Ix0);
            Pyx(2,1)=1-Pyx(1,1);
            Pyx(2,2)=sum(Y2t(Ix1))/length(Ix1);
            Pyx(1,2)=1-Pyx(2,2);
            if Pyx(1,1)==1 && Pyx(2,2)==1
                logp=[0 0;
                      0 0];
            else
            logp=log2(Pyx);
            end
            H=dot(Pyx',logp',2);
            Info1(n,ii,iii,NN)=Px*H-log2(Py)*Py';
            
%             if iii==1 && ii==2
%                 Px
%                 Py
%                 Pyx
%                 Inf
%                 dot(Pyx',log2(Pyx),2)
%             end
            
            end
            sum(isnan(Info1(:,ii,iii,NN)))
            ClassError(ii,iii)=mean(ClassError1)
            ClassError_err(ii,iii)=std(ClassError1);
            Info(ii,iii,NN)=nanmean(Info1(:,ii,iii,NN),1)
            Info_err(ii,iii,NN)=nanstd(Info1(:,ii,iii,NN),1);
        end
    end
    end
    data_A.Info=Info;
    data_A.Info_err=Info_err;
    data_A.ClassError=ClassError;
    data_A.ClassError_err=ClassError_err;



end
end


