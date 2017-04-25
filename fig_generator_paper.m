function []=fig_generator_paper(settings,parameters,data_A)
% Generates figures for paper
v2struct(settings)
v2struct(parameters)
v2struct(data_A)
szpar1=size(yNe1tr,3);
szpar2=size(yNe1tr,4);

if alpha_abs==1
    ch_par1=ch_par1*2;
    x_data=ch_par2/2;
else
    x_data=ch_par2;
end
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);


if Figure==2
    if Gamma==1
        % Plot firing patterns
        fig1=figure(1); clf;
        hold on
        plot(firNe1t(:,1)-200,firNe1t(:,2),'r.');
        plot(firNi11t(:,1)-200,firNi11t(:,2)-400,'b.');
        plot(firNi21t(:,1)-200,firNi21t(:,2)-475,'b.');
        hold off

        xlim([0 1000])
        title('2A')
        xlabel('Time (ms)')
        ylabel('Neuron #')
        set(fig1,'PaperUnits','inches','PaperPosition',[0 0 4 3])

        % Plot STH
        fig2=figure(2); clf;
        hold on
        area(smooth(yNe1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','r')
        area(-smooth(yNi1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','b')
        hold off
        legstr={'Exc.','Inh.'};

        xlim([0 1000])
        title('2B')
        xlabel('Time (ms)')
        ylabel('# Neurons firing')
        set(fig2,'PaperUnits','inches','PaperPosition',[0 0 4 3])

        fig3=figure(3); clf;
        hold on
        plot_variance(f1(:,:,p1_t,p2_t)',10*log10(S1err(1,:,p1_t,p2_t,1)),10*log10(S1err(2,:,p1_t,p2_t,1)),[1,0.5,0.5])
        hp1=plot(f1(:,:,p1_t,p2_t),10*log10(pxx1(:,:,p1_t,p2_t)),'r','Linewidth',2);
        hold off

        title('2C')
        xlabel('Frequency (Hz)')
        ylabel('Power (dB)')
        ylim([-30,0]);
        set(fig3,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end
    

    if EI_Gamma==1
        
        % Plot power/frequency
        for ii=1:szpar1
            for iii=1:szpar2
                lb_gam=find(f1(:,:,ii,iii)>20,1)
                pxx_sm=smooth(pxx1(:,:,ii,iii),5,'moving')
                [P_gam,I]=max(pxx1(lb_gam:end,:,ii,iii))
                F_gam=f1(I+lb_gam-1,:,ii,iii);
                P_gam_all(ii,iii)=P_gam;
                F_gam_all(ii,iii)=F_gam;
            end
        end
    
        fig4=figure(4); clf;
            imagesc(15*ch_par1,15*ch_par2,10*log10(P_gam_all))

            title('2D')
            xlabel('I (pA)','FontName','Times New Roman')
            ylabel('E (pA)','FontName','Times New Roman')
            set(gca,'YDir','normal','FontName','Times New Roman') 
            b1=colorbar;
            yl1=ylabel(b1,'Gamma power (dB)');
            set(yl1,'FontName','Times New Roman');
            set(fig4,'PaperUnits','inches','PaperPosition',[0 0 4 3])

        fig5=figure(5); clf;
            imagesc(15*ch_par1,15*ch_par2,F_gam_all)

            title('2E')
            xlabel('I (pA)','FontName','Times New Roman')
            ylabel('E (pA)','FontName','Times New Roman')
            set(gca,'YDir','normal','FontName','Times New Roman') 
            b2=colorbar;
            yl2=ylabel(b2,'Gamma frequency (Hz)');
            set(yl2)
            set(yl2,'FontName','Times New Roman');
            set(fig5,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end

end

if Figure==3
% Plot firing patterns
    fig6=figure(6); clf;
    hold on
    plot(firNe1t(:,1)-200,firNe1t(:,2),'r.');
    plot(firNi11t(:,1)-200,firNi11t(:,2)-400,'b.');
    plot(firNi21t(:,1)-200,firNi21t(:,2)-475,'b.');
    hold off

    xlim([0 1000])
    title('3A')
    xlabel('Time (ms)')
    ylabel('Neuron #')
    set(fig6,'PaperUnits','inches','PaperPosition',[0 0 4 3])

% Plot STH
    fig7=figure(7); clf;
    hold on
    area(smooth(yNe1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','r')
    area(-smooth(yNi1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','b')
    hold off
    legstr={'Exc.','Inh.'};

    xlim([0 1000])
    title('3B')
    xlabel('Time (ms)')
    ylabel('# Neurons firing')
    set(fig7,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    
    fig8=figure(8); clf;
    p_t=[1,7]
    hold on
    plot_variance(f1(:,1,5,p_t(1))',10*log10(S1err(1,:,5,p_t(1))),10*log10(S1err(2,:,5,p_t(1))),[1,0.5,0.5])
    plot(f1(:,1,5,p_t(1)),10*log10(pxx1(:,1,5,p_t(1))),'r','Linewidth',2);
    plot_variance(f1(:,1,5,p_t(2))',10*log10(S1err(1,:,5,p_t(2))),10*log10(S1err(2,:,5,p_t(2))),[0.5,1,0.5])
    plot(f1(:,1,5,p_t(2)),10*log10(pxx1(:,1,5,p_t(2))),'g','Linewidth',2);

    hold off
    title('3C')
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    ylim([-30,0]);
    set(fig8,'PaperUnits','inches','PaperPosition',[0 0 4 3])

end

% Plot coherence
if Figure==4
    for ii=1:szpar1
        for iii=1:szpar2
            Imin=find(fc1(:,:,ii,iii)>30,1);
            Imax=find(fc1(:,:,ii,iii)<50,1,'last');

                C_gam=mean(C1(Imin:Imax,:,ii,iii));
                C_gam2=mean(C2(Imin:Imax,:,ii,iii));
                C_gam3=mean(C3(Imin:Imax,:,ii,iii));
                C_gam4=mean(C4(Imin:Imax,:,ii,iii));
                C_gam5=mean(C5(Imin:Imax,:,ii,iii));
                C_gam6=mean(C6(Imin:Imax,:,ii,iii));


                size(Cerr)
                C_gam_err(ii,iii)=mean(Cerr(2,:,:,ii,iii)-Cerr(1,:,:,ii,iii))/2;
                C_gam_err2(ii,iii)=mean(Cerr2(2,:,:,ii,iii)-Cerr2(1,:,:,ii,iii))/2;
                C_gam_err3(ii,iii)=mean(Cerr3(2,:,:,ii,iii)-Cerr3(1,:,:,ii,iii))/2;
                C_gam_err4(ii,iii)=mean(Cerr4(2,:,:,ii,iii)-Cerr4(1,:,:,ii,iii))/2;
                C_gam_err5(ii,iii)=mean(Cerr5(2,:,:,ii,iii)-Cerr5(1,:,:,ii,iii))/2;
                C_gam_err6(ii,iii)=mean(Cerr6(2,:,:,ii,iii)-Cerr6(1,:,:,ii,iii))/2;

                C_gam_all(ii,iii)=mean(C_gam);
                C_gam_all2(ii,iii)=mean(C_gam2);
                C_gam_all3(ii,iii)=mean(C_gam3);
                C_gam_all4(ii,iii)=mean(C_gam4);
                C_gam_all5(ii,iii)=mean(C_gam5);
                C_gam_all6(ii,iii)=mean(C_gam6);
        end
    end
    
fig11=figure(11);clf;
    if strcmp(str_ch_par1,'AlphaI1')
        for ii=p2_t
            errorbar(ch_par1,C1(:,ii),Cerr(:,ii)/sqrt(tr),'k','LineWidth',1);
            xlabel('Alpha amplitude')
            ylabel('Gamma coherence')
        end
    else
        set(gca,'XTick',[-180,0,180])
        ColorSet = [linspace(0.9,0,length(ch_par2))',linspace(0.9,0,length(ch_par2))',linspace(0.9,0,length(ch_par2))'];
        set(gca, 'ColorOrder', ColorSet);
        colormap gray
        colormap(flipud(colormap));
        hold all
        for ii=p2_t
        errorbar(ch_par1/pi*180,C_gam_all(:,ii),C_gam_err(:,ii)/sqrt(tr),'k','LineWidth',2);
        errorbar(ch_par1/pi*180,C_gam_all4(:,ii),C_gam_err4(:,ii)/sqrt(tr),'k--','LineWidth',2);

        end
        hold off
        title('4A')
        ylim([0,1])
        xlabel('Alpha phase difference (deg.)')
        ylabel('Gamma coherence')
        set(fig11,'PaperUnits','inches','PaperPosition',[0 0 4 3])

fig12=figure(12);clf;           
        hold all
        errorbar(ch_par1/pi*180,mean(frateNe2(:,p2_t,:),3),std(frateNe2(:,p2_t,:),[],3),'r','LineWidth',2)
        errorbar(ch_par1/pi*180,mean(frateNi2(:,p2_t,:),3),std(frateNi2(:,p2_t,:),[],3),'b','LineWidth',2)
        hold off
        ylim([0 35])
        set(gca,'XTick',[-180,0,180])
        set(fig12,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('4B')
        xlabel('Alpha phase difference')
        ylabel('Firing rate')

fig13=figure(13);clf;   
        set(0,'defaultAxesColorOrder',[0 0 1;0 1 0;1 0 0])
        if alpha_abs==1
            x_data=ch_par2/2;
        else
            x_data=ch_par2;
        end
        hold all
        for ii=[3,7]
            errorbar(x_data*15,C_gam_all(ii,:),C_gam_err(ii,:)/sqrt(tr),'LineWidth',2);
        end
        hold off

        ylim([0 1])
        xlim([-5 50])
        if CohE==1
            xlim([-5 20])
        end
        title('4F')
        xlabel('Alpha amplitude (pA)')
        ylabel('Gamma coherence')
        set(fig13,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        
fig14=figure(14);clf;
        hold all
        errorbar(ch_par1/pi*180,C_gam_all2(:,p2_t),C_gam_err2(:,p2_t)/sqrt(tr),'k','LineWidth',2);
        errorbar(ch_par1/pi*180,C_gam_all5(:,p2_t),C_gam_err5(:,p2_t)/sqrt(tr),'k--','LineWidth',2);
        hold off

        set(gca,'XTick',[-180,0,180])
        title('4D')
        xlabel('Alpha phase difference (deg.)')
        ylabel('Gamma phase coherence')
        ylim([0,1])
        set(fig14,'PaperUnits','inches','PaperPosition',[0 0 4 3])

fig15=figure(15);clf;
        hold all
        errorbar(ch_par1/pi*180,C_gam_all3(:,p2_t),C_gam_err3(:,p2_t)/sqrt(tr),'k','LineWidth',2);
        errorbar(ch_par1/pi*180,C_gam_all6(:,p2_t),C_gam_err6(:,p2_t)/sqrt(tr),'k--','LineWidth',2);
        hold off

        set(gca,'XTick',[-180,0,180])
        title('4C')
        xlabel('Alpha phase difference (deg.)')
        ylabel('Gamma amplitude coherence')
        ylim([0,1])

        set(fig15,'PaperUnits','inches','PaperPosition',[0 0 4 3])
            
fig16=figure(16);clf;
        hold all
        plot_variance(fc1(:,:,3,p2_t),Cerr(1,:,:,3,p2_t),Cerr(2,:,:,3,p2_t),[0.5,0.5,1])
        plot(fc1(:,:,3,p2_t),C1(:,:,3,p2_t),'b','Linewidth',2)
        plot_variance(fc1(:,:,7,p2_t),Cerr(1,:,:,7,p2_t),Cerr(2,:,:,7,p2_t),[0.5,1,0.5])
        plot(fc1(:,:,7,p2_t),C1(:,:,7,p2_t),'g','Linewidth',2)
        title('4E')
        xlabel('Frequency')
        ylabel('Coherence')
        set(fig16,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end
end 

if Figure==5
    fig17=figure(17); clf;
    errorbar(V11m(:,:,1),V11err(:,:,1),'LineWidth',2)
    xlim([0,9])

    ylim([0,0.3])
    title('5A')
    xlabel('Area 1 gamma phase (deg.)')
    ylabel('Normalized area 1 spike counts')
    set(gca,'XTick',[1,4.5,8])
    set(gca,'XTickLabel',[-360,-180,0])
    set(fig17,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    
    fig18=figure(18); clf;
    errorbar(V22m(:,:,1),V22err(:,:,1),'LineWidth',2)
    xlim([0,9])
    ylim([0,0.3])
	title('5B')
    xlabel('Area 2 gamma phase (deg.)')
    ylabel('Normalized area 2 spike counts')
    set(gca,'XTick',[1,4.5,8])
    set(gca,'XTickLabel',[-360,-180,0])

    set(fig18,'PaperUnits','inches','PaperPosition',[0 0 4 3])

fig19=figure(19); clf;
    errorbar(V12m(:,:,1),V12err(:,:,1),'LineWidth',2)
    xlim([0,9])
    ylim([0,0.2])
	title('5C')
    xlabel('Area 2 gamma phase (deg.)')
    ylabel('Normalized area 1 spike counts')
    set(fig19,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    set(gca,'XTick',[1,4.5,8])
    set(gca,'XTickLabel',[-360,-180,0])
    
fig20=figure(20); clf;
    errorbar(probm(:,:,1),proberr(:,:,1),'LineWidth',2)
    xlim([0,9])
    ylim([0,0.25])
    set(gca,'XTick',[1,4.5,8])
    set(gca,'XTickLabel',[-360,-180,0])
    title('5D')
    xlabel('Area 1 gamma phase (deg.)')
    ylabel(sprintf('Normalized probability of area 2 spike \n following area 1 spike'))
    set(fig20,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    legend('1','2','3','4','5','6','7','8','9')
            
fig21=figure(21); clf;
    errorbar(probm(:,:,2),proberr(:,:,2),'LineWidth',2)
    xlim([0,9])
    ylim([0,0.25])
    title('5E')
    xlabel('Area 2 gamma phase (deg.)')
    ylabel(sprintf('Normalized probability of area 2 spike \n following area 1 spike'))
    set(fig21,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    set(gca,'XTick',[1,4.5,8])
    set(gca,'XTickLabel',[-360,-180,0])
end
 

% Plot stimulus figures
if Figure==6
    if stim_phase==0
        [~,Ish]=min(mean(frateNe1st(11:20,:,p2_t),1),[],2);
        frateNe1st=circshift(frateNe1st,[0,-Ish,0]);

        % Plot STH
        fig22=figure(22);clf;
        ha=tight_subplot(2,1,0,0.1,0.1);
        axes(ha(1));
        hold on

        area(smooth(yNe1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','r')
        area(-smooth(yNi1tr(:,tr_t,p1_t,p2_t),10),'FaceColor','b')
        hold off
        ylabel('# Neurons firing')
        set(gca,'YTick',[0,10,20])
        set(gca,'YTickLabel',[0,10,20])
        set(gca,'Position',[0.12,0.31,0.78,0.65])
        box off

        axes(ha(2));
        stim_curr=ones(1,length(yNe1tr))*Emns;
        stim_curr((stimonset(p1_t)-p_in):end)=Ems;
        plot(10*stim_curr,'k','LineWidth',2)
        ylim([20,60])
        title('6A')
        ylabel('Input current (pA)')
        xlabel('Time (ms)')
        box off
        set(gca,'Position',[0.12,0.12,0.78,0.18])

        set(fig22,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    
    fig23=figure(23);clf;
        hold all
        errorbar(-180:36:180,mean(frateNe1st(1:10,:,p2_t),1),std(frateNe1st(1:10,:,p2_t),[],1)/sqrt(tr/2),'k','Linewidth',2)
        errorbar(-180:36:180,mean(frateNe1st(11:20,:,p2_t),1),std(frateNe1st(11:20,:,p2_t),[],1)/sqrt(tr/2),'k--','Linewidth',2)
        hold off
        set(gca,'XTick',[-180,0,180])
        set(fig23,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('6B')
        xlabel('Alpha phase (deg.)')
        ylabel('Firing rate')
            
    fig24=figure(24);clf;
        errorbar(-180:36:180,mean(frateNe1st(1:10,:,p2_t)-frateNe1st(11:20,:,p2_t),1),std(frateNe1st(1:10,:,p2_t)-frateNe1st(11:20,:,p2_t),[],1)/sqrt(tr/2),'k','Linewidth',2)
        set(gca,'XTick',[-180,0,180])
        set(fig24,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        box off
        title('6C')
        xlabel('Alpha phase (deg.)')
        ylabel('Firing rate difference')

        frateNe2stc=circshift(frateNe2st,[0,-Ish,0]);
        p1_t=6;            
    fig25=figure(25); clf;
        hold all
        errorbar(ch_par2/pi*180,mean(frateNe2stc(1:10,p1_t,:),1),std(frateNe2stc(1:10,p1_t,:),[],1)/sqrt(tr/2),'k','Linewidth',2)
        errorbar(ch_par2/pi*180,mean(frateNe2stc(11:20,p1_t,:),1),std(frateNe2stc(11:20,p1_t,:),[],1)/sqrt(tr/2),'k--','Linewidth',2)
        hold off

        set(gca,'XTick',[-180,0,180])
        set(fig25,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('6E')
        xlabel('Alpha phase difference (deg.)')
        ylabel('Firing rate')

    fig26=figure(26); clf;
        errorbar(ch_par2/pi*180,mean(frateNe2stc(1:10,p1_t,:)-frateNe2stc(11:20,p1_t,:),1),std(frateNe2stc(1:10,p1_t,:)-frateNe2stc(11:20,p1_t,:),[],1)/sqrt(tr/2),'k','Linewidth',2)
        set(gca,'XTick',[-180,0,180])
        ylim([0,25])
        box off

        set(fig26,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('6F')
        xlabel('Alpha phase difference (deg.)')
        ylabel('Firing rate difference')

    fig27=figure(27);
        stim_resp_all=squeeze(mean(frateNe2st(1:10,:,:),1)-mean(frateNe2st(11:20,:,:),1));
        imagesc(-180:45:180,-180:36:180,stim_resp_all)
        set(gca,'XTick',[-180,0,180])
        set(gca,'YTick',[-180,0,180])
        title('6D')
        ylabel('Stimulus w.r.t. alpha phase area 1')
        xlabel('Alpha phase difference (deg.)')
        h1=colorbar;
        ylabel(h1,'Firings rate difference area 2')
        set(gcf,'renderer','painters')
        set(fig27,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end
end

% Plot mutual information figures
if Figure==7
    if MI==1
        fig28=figure(28);clf;
            colormap gray
            colormap(flipud(colormap));
            errorbar(ch_par1/pi*180,Info(:,1,6),Info_err(:,1,6),'k','LineWidth',2);
            set(gca,'XTick',[-180,0,180])

            set(fig28,'PaperUnits','inches','PaperPosition',[0 0 4 3])
            ylim([0,1])
            box off
            title('7C')
            xlabel('Alpha phase difference (deg.)')
            ylabel('Mutual information')

        fig29=figure(29);
            NN_all=round(logspace(2.3,0,10));
            minMI=min(squeeze(Info),[],1);
            maxMI=(max(squeeze(Info),[],1));
            C=bsxfun(@rdivide,squeeze(Info)-repmat(minMI,9,1),maxMI-minMI);

            C(isnan(C))=1;

            surf(1:10,-180:45:180,squeeze(Info),C)
            view(130,30)
            set(gca, 'YTick', [-180,0,180])
            set(gca, 'XTick', 1:2:10)
            set(gca, 'XTickLabel', NN_all(1:2:10))
            
            title('7D')
            ylabel('Alpha phase difference (deg.)')
            xlabel('# of neurons used')
            h1=get(gca,'xlabel');
            set(h1,'rotation',25)
            h2=get(gca,'ylabel');
            set(h2,'rotation',-17)
            get(h1,'Position')
            set(h1,'Position',get(h1,'Position') + [.3 0 -0.05])
            get(h1,'Position')

            set(h2,'Position',get(h2,'Position') - [0 0 0])

            zlabel('Mutual information')
            set(fig29,'PaperUnits','inches','PaperPosition',[0 0 4 3])
            
    elseif MI==2         
        fig30=figure(30);
            imagesc(-180:45:180,-180:36:180,Info(:,[2:11,1])')
            set(gca,'XTick',[-180,0,180],'FontName','Times New Roman')
            set(gca,'YTick',[-180,0,180])
            title('7B')
            ylabel('Stimulus w.r.t. alpha phase area 1','FontName','Times New Roman')
            xlabel('Alpha phase difference (deg.)','FontName','Times New Roman')
            h1=colorbar;
            ylabel(h1,'Mutual information','FontName','Times New Roman')
            set(gcf,'renderer','painters')
            set(fig30,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end
            
end
    
% Plot Granger causality figures
if Figure==8
    if coherence==1
        for ii=1:szpar1
            for iii=1:szpar2
                Imin=find(fc1(:,:,ii,iii)>30,1);
                Imax=find(fc1(:,:,ii,iii)<50,1,'last');

                C_gam=mean(C1(Imin:Imax,:,ii,iii));
                C_gam2=mean(C2(Imin:Imax,:,ii,iii));
                C_gam3=mean(C3(Imin:Imax,:,ii,iii));
                C_gam4=mean(C4(Imin:Imax,:,ii,iii));
                C_gam5=mean(C5(Imin:Imax,:,ii,iii));
                C_gam6=mean(C6(Imin:Imax,:,ii,iii));


                size(Cerr)
                C_gam_err(ii,iii)=mean(Cerr(2,:,:,ii,iii)-Cerr(1,:,:,ii,iii))/2;
                C_gam_err2(ii,iii)=mean(Cerr2(2,:,:,ii,iii)-Cerr2(1,:,:,ii,iii))/2;
                C_gam_err3(ii,iii)=mean(Cerr3(2,:,:,ii,iii)-Cerr3(1,:,:,ii,iii))/2;
                C_gam_err4(ii,iii)=mean(Cerr4(2,:,:,ii,iii)-Cerr4(1,:,:,ii,iii))/2;
                C_gam_err5(ii,iii)=mean(Cerr5(2,:,:,ii,iii)-Cerr5(1,:,:,ii,iii))/2;
                C_gam_err6(ii,iii)=mean(Cerr6(2,:,:,ii,iii)-Cerr6(1,:,:,ii,iii))/2;

                C_gam_all(ii,iii)=mean(C_gam);
                C_gam_all2(ii,iii)=mean(C_gam2);
                C_gam_all3(ii,iii)=mean(C_gam3);
                C_gam_all4(ii,iii)=mean(C_gam4);
                C_gam_all5(ii,iii)=mean(C_gam5);
                C_gam_all6(ii,iii)=mean(C_gam6);
            end
        end

        fig31=figure(31);clf;
            set(gca,'XTick',[-180,0,180])
            ColorSet = [linspace(0.9,0,length(ch_par2))',linspace(0.9,0,length(ch_par2))',linspace(0.9,0,length(ch_par2))'];
            set(gca, 'ColorOrder', ColorSet);
            colormap gray
            colormap(flipud(colormap));
            hold all
            
            for ii=p2_t
                errorbar(ch_par1/pi*180,C_gam_all(:,ii),C_gam_err(:,ii)/sqrt(tr),'k','LineWidth',2);
                errorbar(ch_par1/pi*180,C_gam_all4(:,ii),C_gam_err4(:,ii)/sqrt(tr),'k--','LineWidth',2);
            end
            
            hold off
            ylim([0,1])
            title('8B')
            xlabel('Alpha phase difference (deg.)')
            ylabel('Gamma coherence')
            set(fig31,'PaperUnits','inches','PaperPosition',[0 0 4 3])
    end
    % Granger causality
    if granger==1

    fig32=figure(32);clf;
        p1_tt=3;
        hold all
        plot_variance(lam{p1_tt},FSLO12{p1_tt},FSUP12{p1_tt},[0.4,0.4,1])
        ha1(1)=plot(lam{p1_tt},f12{p1_tt},'b');
        plot_variance(lam{p1_tt},FSLO21{p1_tt},FSUP21{p1_tt},[1,0.4,0.4])
        ha1(2)=plot(lam{p1_tt},f21{p1_tt},'r');
        hold off

        set(fig32,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('8C')
        xlim([0 100])
        xlabel('f (Hz)');
        ylabel('Granger causality');

    fig33=figure(33);clf;
        p1_tt=7;
        hold all
        plot_variance(lam{p1_tt},FSLO12{p1_tt},FSUP12{p1_tt},[0.4,0.4,1])
        ha1(1)=plot(lam{p1_tt},f12{p1_tt},'b');
        plot_variance(lam{p1_tt},FSLO21{p1_tt},FSUP21{p1_tt},[1,0.4,0.4])
        ha1(2)=plot(lam{p1_tt},f21{p1_tt},'r');
        hold off

        set(fig33,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('8D')
        xlim([0 100])
        xlabel('f (Hz)');
        ylabel('Granger causality');

    fig34=figure(34);clf;
        p1_tt=5;
        hold all
        plot_variance(lam{p1_tt},FSLO12{p1_tt},FSUP12{p1_tt},[0.4,0.4,1])
        ha1(1)=plot(lam{p1_tt},f12{p1_tt},'b');
        plot_variance(lam{p1_tt},FSLO21{p1_tt},FSUP21{p1_tt},[1,0.4,0.4])
        ha1(2)=plot(lam{p1_tt},f21{p1_tt},'r');
        hold off

        set(fig34,'PaperUnits','inches','PaperPosition',[0 0 4 3])
        title('8E')
        xlim([0 100])
        xlabel('f (Hz)');
        ylabel('Granger causality');
    end
end
end


