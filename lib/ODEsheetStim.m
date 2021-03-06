function dYdt=ODEsheetStim(t,y,parameters)

nsq=parameters.n^2;


SigThresh=parameters.SigThresh;
SigSteepness=parameters.SigSteepness;
tauPy=parameters.tauPy;
tauInh=parameters.tauInh;


Py=y(1:nsq);
Inh=y(nsq+1:2*nsq);
Stim=parameters.stimampl*(sin(2*pi*parameters.stimfreq*t)+1)/2;
Py(parameters.stimlocs)=Stim;

dPydt     =(-Py     + Sigm(parameters.Py2Py*Py      - parameters.Inh2Py*Inh + parameters.PyInput ,  SigThresh,SigSteepness))./tauPy;
dInhdt    =(-Inh    + Sigm(parameters.Py2Inh*Py     - parameters.Inh2Inh*Inh             + parameters.InhInput,                SigThresh,SigSteepness))./tauInh;

dPydt(parameters.stimlocs)=0;

dYdt = [dPydt;dInhdt];
end