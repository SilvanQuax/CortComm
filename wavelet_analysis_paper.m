function [data_A]=wavelet_analysis_new(settings,parameters,firings_all,data_A)

v2struct(settings);
v2struct(parameters);
v2struct(data_A);
szpar1=size(firings_all,1);
szpar2=size(firings_all,2);
scales=logspace(1,log10(150));
F = scal2frq(scales,'cmor1-1',1/1000)
t_plot=length(p_in:rt);
% if exist('coefs1','var')==0
figure(401);clf;
for ii=1:szpar1
    ii
         iii=1;
        for tr1=1:tr
             
            firings=firings_all{ii,p2_t,tr1};
            if npop>=1
                yNe1tr(:,tr1,ii,iii) =histc(firings(firings(:,2)<=Ne,1),p_in:rt);
                yNi1tr(:,tr1,ii,iii) =histc(firings(firings(:,2)> npop*Ne & firings(:,2)<=npop*Ne+Ni1,1),p_in:rt)...
                +histc(firings(firings(:,2)> npop*Ne+npop*Ni1 & firings(:,2)<=npop*Ne+npop*Ni1+Ni2,1),p_in:rt);

                x1=yNe1tr(:,tr1,ii,iii);
                x2=yNi1tr(:,tr1,ii,iii);
                coefs1(:,:,tr1,ii,iii) = cwt(x1,scales,'cmor1-1','plot');
                coefs2(:,:,tr1,ii,iii) = cwt(x2,scales,'cmor1-1','plot');
                SC(:,:,tr1,ii,iii)=wscalogram('image',coefs1(:,:,tr1,ii,iii),'scales',scales);
                SC2(:,:,tr1,ii,iii)=wscalogram('image',coefs2(:,:,tr1,ii,iii),'scales',scales);

                phases1(:,:,tr1,ii,iii)=-angle(coefs1(:,:,tr1,ii,iii));
            end

            if npop>=2 
                yNe2tr(:,tr1,ii,iii) =histc(firings(firings(:,2)> Ne & firings(:,2)<=2*Ne,1),p_in:rt);
                yNi2tr(:,tr1,ii,iii) =histc(firings(firings(:,2)> npop*Ne+Ni1 & firings(:,2)<=npop*Ne+2*Ni1,1),p_in:rt)...
                +histc(firings(firings(:,2)> npop*Ne+npop*Ni1+Ni2 & firings(:,2)<=npop*Ne+npop*Ni1+2*Ni2,1),p_in:rt);
            end
        
    end
end


    if npop>=2 

        for ii=1:szpar1
            iii=1;
                for tr1=1:tr
        x21=yNe2tr(:,tr1,ii,iii);
        coefs21(:,:,tr1,ii,iii)= cwt(x21,scales,'cmor1-1','plot');
        phases21(:,:,tr1,ii,iii)=-angle(coefs21(:,:,tr1,ii,iii));
                end
            
        end
    data_A.coefs21=coefs21;
    data_A.phases21=phases21;
    data_A.yNe2tr=yNe2tr;
    end
    data_A.coefs1=coefs1;
    data_A.phases1=phases1;
    data_A.yNe1tr=yNe1tr;
    data_A.SC=SC;
    data_A.SC2=SC2;

if V1V2==1
% V1-V2 spikes locked to gamma phase
nbins=9;
iii=0;
for p1_t=[3,7]
    iii=iii+1;
    for phas=[1,2]
for tr1=1:tr
    [~,bin]=histc(phases1(19,1:(end-50),tr1,p1_t,1),linspace(-pi,pi,nbins));
    [~,bin2]=histc(phases21(19,1:(end-50),tr1,p1_t,1),linspace(-pi,pi,nbins));
  
    for n=1:(nbins-1)
        Ind=find(bin==n);
        Ind2=find(bin2==n);        
        V11(n,tr1)=mean(yNe1tr(Ind,tr1,p1_t,1));
        V22(n,tr1)=mean(yNe2tr(Ind2,tr1,p1_t,1));
        V12(n,tr1)=mean(yNe1tr(Ind2,tr1,p1_t,1));
%         for ii=1:50
%         AA=yNe2tr(Ind+ii,tr1,p1_t,1)./yNe1tr(Ind,tr1,p1_t,1);
%         V12p1(n,ii,tr1)=mean(AA(isfinite(AA)));
%         AA2=yNe2tr(Ind2+ii,tr1,p1_t,1)./yNe1tr(Ind2,tr1,p1_t,1);
%         V12p2(n,ii,tr1)=mean(AA2(isfinite(AA2)));
%         end
    end
end
% Chances
swin=4;
kk_all=zeros(swin,8,10);
% [row,col]=find(S(401:450,1:400)>0.25);
% AA=size(row)
for tr1=1:tr
    firings1=firings_all{p1_t,p2_t,tr1};
    firings1(:,1)=firings1(:,1)-p_in+1;
    firings=[firings1(firings1(:,1)>0,1),firings1(firings1(:,1)>0,2)];
    [~,bin1]=histc(phases1(19,1:(end),tr1,p1_t,1),linspace(-pi,pi,nbins));
    [~,bin2]=histc(phases21(19,1:(end),tr1,p1_t,1),linspace(-pi,pi,nbins));

    I1=find(firings(:,2)<=Ne); % Firings pop1
    timings_V1=firings(I1,1); % Spike times pop1
    
    fV2=ismember(firings(:,2),401:800); % Firings pop2
    firingsV2=firings(fV2,:); % 
    MemI=[];
    timings_V2a=[];
    for n=1:(nbins-1)
        if phas==1
        Ind=find(bin2==n); % Gamma phase bin indices pop1
        elseif phas==2
        Ind=find(bin2==n); % Gamma phase bin indices pop2
        end
    MemI(:,n)=ismember(firingsV2(:,1),Ind); % Timings spikes pop2 in bin n
    Nfir(n)=sum(MemI(:,n),1); % # pop2 firings per bin
    end
    
    len=min(Nfir); % For normalizing pop2 spikes
    for n=1:(nbins-1) % Determining timings V2a (normalizing for V2 spikes)
        if phas==1
        Ind=find(bin2==n);
        elseif phas==2
        Ind=find(bin2==n);
        end


        tfir2=firingsV2.*repmat(MemI(:,n),1,2); % pop2 firings in bin n
        IInd=find(tfir2(:,1)); % Timings pop2 in bin n
        s1=length(IInd); % # of firings pop2 in bin n
        IndV2=randperm(s1,len); % Reduce number of firings to len
        firV2=tfir2(IInd(IndV2),:); % Select firings at timings IInd(IndV2)
        timings_V2a=[timings_V2a;firV2];
        size(firV2)
        frac(n)=s1/len;
    end
    frax=mean(frac);
    len_all(iii,phas,tr1)=len;
    
    for n=1:(nbins-1)
        if phas==1
        Ind=find(bin1==n);
        elseif phas==2
        Ind=find(bin2==n);
        end
       
        kk=0;
        kkk=0;
        n
        tfir=0;
        pp=0;
        kk11=zeros(1,swin);
        iall=0;
% timings_V2a=firingsV2; % Disable to compensate for different amount of V2 spikes for different gamma phases
        for t=101:140
            for tt=501:540
                I1=find(firings(:,2)==t);
                I2=find(timings_V2a(:,2)==tt);
                
                timings_V1n=firings(I1,1); % Spike timings neuron t
                timings_V2n=timings_V2a(I2,1);% Spike timings neuron tt
                pp=pp+length(I2);
                btim=timings_V1n(ismember(timings_V1n,Ind)); % Spike timings neuron t in bin n  
                btim2=timings_V2n(ismember(timings_V2n,Ind)); % Spike timings neuron tt in bin n
                
                if isempty(btim)==0
                    kk2=0;
                    for ii=1:swin
                        % Testing whether V1 stike is followed by V2 spike 1-4
                        % ms afterwards
                        III=ismember(timings_V1n,timings_V2n+ii);
                        iall=iall+sum(III);
                        kk1(ii)=sum(ismember(timings_V1n(III),Ind));
%                         timings_V2n(III)=0;
                    end
               
                    % Normalise for amount of V1 spikes 
                    
                    kk=kk+sum(kk1)/length(btim);
                    kk11=kk11+kk1./length(btim);
                end
                
%                 if phas==1
                         tfir=tfir+length(btim);   
%                 elseif phas==2
%                          tfir=tfir+length(btim2);   
%                 end

            end
        end
        ialll(n,tr1)=iall;
        kk_all(:,n,tr1)=kk11/1600;
        prob(n,tr1)=kk/1600;
        kka(n,tr1)=kk;
        tfira(n,tr1)=tfir;
    end
end
iallm(:,iii,phas)=mean(ialll,2)
kk_allm=mean(kk_all,3);
% Shift gamma phase so highest spike phase is on right
V11mm=mean(V11,2);
[~,Im1]=max(V11mm);
V22mm=mean(V22,2);
[~,Im2]=max(V22mm);

        if phas==1
        probm(:,iii,phas)=circshift(mean(prob,2),nbins-Im1-1)/sum(mean(prob,2));
        proba(:,iii,phas)=circshift(mean(prob,2),nbins-Im1-1);

        proberr(:,iii,phas)=circshift(std(prob,[],2),nbins-Im1-1)/sum(mean(prob,2))/sqrt(tr);
        probaerr(:,iii,phas)=circshift(std(prob,[],2),nbins-Im1-1)/sqrt(tr);

        tfirm(:,iii,phas)=circshift(mean(tfira,2),nbins-Im1-1);
        kkm(:,iii,phas)=circshift(mean(kka,2),nbins-Im1-1);

        elseif phas==2
        probm(:,iii,phas)=circshift(mean(prob,2),nbins-Im2-1)/sum(mean(prob,2));
        proba(:,iii,phas)=circshift(mean(prob,2),nbins-Im2-1);

        proberr(:,iii,phas)=circshift(std(prob,[],2),nbins-Im2-1)/sum(mean(prob,2))/sqrt(tr);
        probaerr(:,iii,phas)=circshift(std(prob,[],2),nbins-Im2-1)/sqrt(tr);

        tfirm(:,iii,phas)=circshift(mean(tfira,2),nbins-Im2-1);
        kkm(:,iii,phas)=circshift(mean(kka,2),nbins-Im2-1);

        end
    end
V11m(:,iii)=circshift(V11mm,nbins-Im1-1)/sum(V11mm);
V11err(:,iii)=circshift(std(V11,[],2),nbins-Im1-1)/sum(V11mm)/sqrt(tr);
V22m(:,iii)=circshift(V22mm,nbins-Im2-1)/sum(V11mm);
V22err(:,iii)=circshift(std(V22,[],2),nbins-Im2-1)/sum(V11mm)/sqrt(tr);
V12m(:,iii)=circshift(mean(V12,2),nbins-Im2-1)/sum(V11mm);
V12a(:,iii)=circshift(mean(V12,2),nbins-Im2-1);

V12err(:,iii)=circshift(std(V12,[],2),nbins-Im2-1)/sum(V11mm)/sqrt(tr);
end
data_A.V11m=V11m;
data_A.V11err=V11err;
data_A.V22m=V22m;
data_A.V22err=V22err;
data_A.V12m=V12m;
data_A.V12err=V12err;
data_A.probm=probm;
data_A.proba=proba;
data_A.proberr=proberr;
data_A.probaerr=probaerr;

data_A.v12=V12;
mns=[mean(tfirm(:,1,1)),mean(tfirm(:,1,2)),mean(kkm(:,1,1)),mean(kkm(:,1,2))]
len_all

end