function parameters=getParam(n,CeRem,CeLoc,CeLocI)


parameters.n=n;%sheet is nxn

parameters.Py2Py=10*speye(n^2)+.15*CeLoc+.05*CeRem;%conn matrix Py to Py
%parameters.Py2Py=10*speye(n^2)+.1*CeLoc+.1*CeRem;%conn matrix Py to Py
parameters.Inh2Py=25*speye(n^2);%in this parameter set the increase in this prolongs the onset & prevents too much going to the upper fixed point
parameters.Py2Inh=0.1*CeLocI+15*speye(n^2);%
parameters.Inh2Inh=0.0*speye(n^2);

parameters.PyInput=-2.5*(ones(n^2,1));
parameters.InhInput=-5*ones(n^2,1);
    


parameters.tauPy=1*ones(n^2,1)/25;%time scale of the populations
parameters.tauInh=0.5*ones(n^2,1)/25;
%parameters.tauInFast=0.1*ones(n^2,1);

parameters.SigThresh=4*ones(n^2,1);%sigmoid parameters
parameters.SigSteepness=1*ones(n^2,1);

parameters.h=1/(500);

end