function parameters=getParamDelay(n,CeRem,CeLoc,CeLocI,h,propspeed)


parameters.n=n;%sheet is nxn


parameters.Inh2Py=25*speye(n^2);%in this parameter set the increase in this prolongs the onset & prevents too much going to the upper fixed point


parameters.PyInput=-2.5*(ones(n^2,1));
parameters.InhInput=-5*ones(n^2,1);
    


parameters.tauPy=1*ones(n^2,1)/25;%time scale of the populations
parameters.tauInh=0.5*ones(n^2,1)/25;
%parameters.tauInFast=0.1*ones(n^2,1);

parameters.SigThresh=4*ones(n^2,1);%sigmoid parameters
parameters.SigSteepness=1*ones(n^2,1);

parameters.h=h;

hdtime=h*propspeed*1000*1000; %stepsize times delaytime

CeLocDelay=getDelayMatrix(CeLoc,n,@distTorus,hdtime);
CeLocIDelay=getDelayMatrix(CeLocI,n,@distTorus,hdtime);
CeRemDelay=getDelayMatrix(CeRem,n,@distTorus,hdtime);


maxTS=max([length(CeLocDelay),length(CeRemDelay)]);
minTS=min([length(CeLocDelay),length(CeRemDelay)]);


parameters.Py2Py{1}=10*speye(n^2)+.15*CeLocDelay{1}+.05*CeRemDelay{1};%conn matrix Py to Py at time step 1 distance
for k=2:maxTS
    if k>minTS
        parameters.Py2Py{k}=.05*CeRemDelay{k};%conn matrix Py to Py at time step minTS until maxTS distance
    else
        parameters.Py2Py{k}=.15*CeLocDelay{k}+.05*CeRemDelay{k};%conn matrix Py to Py  at time step 2 until minTS distance
    end
end

parameters.Py2Inh{1}=0.1*CeLocIDelay{1}+15*speye(n^2);%
for k=2:length(CeLocIDelay)
    parameters.Py2Inh{k}=0.1*CeLocIDelay{k};%
end

parameters.delaysteps=maxTS;
end