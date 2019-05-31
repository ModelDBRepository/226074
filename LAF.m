%% This code runs based on the previously published 
%% paper: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003787
%% and code: https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=155565
%% note that all files in lib essentially stayed the same. So if you want to run a bigger sheet, please see previous code for the generation of the connectivity matrix.
%% contact yujiang.wang@ncl.ac.uk for more details

clear all
close all
addpath('lib')
%load('Conns_n150.mat') I recommend generating the connectivity matrix once
%and saving it and loading it in future...

tic;
n=150;

%calculate Connectivity Matrix for Py->Py 
CeLoc=GaussianLocConnFunc(n,@distTorus,5);%r_loc=500 micrometre, hence the standard deviation of the gaussian is 250 micrometre, which corresponds to 5 units.
CeLoc(CeLoc>0)=1;
%calculate Connectivity Matrix for Py->In
CeLocI=GaussianLocConnFunc(n,@distTorus,5);%r_loc=500 micrometer
CeLocI(CeLocI>0)=1;
toc;
%remote conn
tic;
nOut=round(mean(sum(CeLoc,1))*4/6);%number of outgoing connections per mini column
%patchSize*numPatches should be > nOut!!!
remRad=75;%3750 micrometers
nM=10;
patchSize=round(5^2*3.14/2);
numPatches=6;
nOverlap=3;
CeRem=ConnPatchyRemOverlap(n,nM,patchSize,numPatches,remRad,nOut,nOverlap,@distTorus,@makeCellClusterToroidal);
toc;



%% parameters
n=150;
tinterp=1;

parameters=getParam(n,CeRem,CeLoc,CeLocI);


%% prescan to get initial state
        %prescan===============================
        tend=1.5;
        
        T=0:parameters.h*tinterp:tend;%for plotting
        nIt=(tend)/parameters.h+1;
        tic
            %using all 0 initial condition for Py
            initCond=zeros(2*n^2,1);
            parameters.NValue=getNoise(nIt,n);

            Y=runSheet(initCond,parameters);
            bgInitc=Y(end,:);
        toc

        
        %% real run
        
        parameters=getParam(n,CeRem,CeLoc,CeLocI);
        tend=15;
        T=0:parameters.h*tinterp:tend;%for plotting
        nIt=(tend)/parameters.h+1;
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ClusterNum=15; %change this number to change the number of subclusters
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        PercentHetP=0.05;
        P=-2.5;
        parameters.PyInput=P*(ones(n^2,nIt));

        CellLocVAll=[];
        CellLocAll=zeros(n,n);
        r=25;%onset gap in units of 2ms
        Pt=2;%target P
        
        for cn=1:ClusterNum
            [CellLoc,CellLocV] = makeCellCluster(1,PercentHetP/ClusterNum,n);%to scan: clustering coeff and percent of bad cells
            CellLocVAll=[CellLocVAll, CellLocV];
            CellLocAll=CellLocAll+CellLoc;
            

            Ramp=[P*ones(1,1000), P:abs((P-Pt))/(1500+r):Pt, Pt*ones(1,5000-r)]; 
           
            parameters.PyInput(CellLocV,:)=repmat(Ramp,length(CellLocV),1);
            
        end
        CellLocVAll=unique(CellLocVAll);
        
        imagesc(CellLocAll)
        
        %Ramp=[P*ones(1,250), P:abs((P-1))/500:1, 1*ones(1,1750)]; 
        


        
 %%       
        %real runs===============================
        


        %using initial condition all 1
        tic
        initCond=bgInitc;
        
        parameters.NValue=getNoise(nIt,n);
        
        
        Y=runSheetPRamp(initCond,parameters);
        toc
        %%
        Py=Y(1:tinterp:end,1:n^2); 
        
        Inh=Y(1:tinterp:end,n^2+1:end);
        SCInput=parameters.NValue(:,1:tinterp:end);
        PyInput=parameters.PyInput(:,1:tinterp:end);
        LFP=parameters.Py2Py*double(Py') - parameters.Inh2Py*double(Inh') + 0.1*SCInput;% + PyInput + SCInput;
        LFP=single(LFP');

        mPy=mean(Py,2);sPy=std(Py,0,2);
        mLFP=mean(LFP,2);

        [mMacroCol,mMacroColLFP]=meanMacroCol(n,10,Py,LFP);


        
%%  plot
load('MayColourMap')


figure(1)
plot(T,mLFP)
title('Raw model output LFP')

samplingF=1/T(2);
filtmLFP=FilterEEG(mLFP, 1, samplingF, 'high', 3);
figure(2)
plot(T,filtmLFP)
title('Filtered model output LFP - with DC shift')



figure(5)
TPts=[1501:1000:5501];
for k=1:length(TPts)
    subplot(1,length(TPts),k)
    imagesc(reshape(Py(TPts(k),:),n,n) )
    %colormap(mycmap)
    caxis([0 0.5])
    title(sprintf('T=%g',T(TPts(k))))
end

