function Y=runSheetPRamp(y0,parameters)

%this is an euler solver with a noise term as subcortical input, which
%takes a ramp in P.


nsq=parameters.n^2;
h=parameters.h;


SigThresh=parameters.SigThresh;
SigSteepness=parameters.SigSteepness;
tauPy=parameters.tauPy;
tauInh=parameters.tauInh;


neq = length(y0);
Nl = size(parameters.NValue,2);
Y = zeros(neq,Nl);

Y(:,1) = y0;
for i = 1:Nl-1 
    Py=Y(1:nsq,i);
    Inh=Y(nsq+1:2*nsq,i);
    
    
    dPydt     =(-Py     + Sigm(parameters.Py2Py*Py      - parameters.Inh2Py*Inh + parameters.PyInput(:,i) + parameters.NValue(:,i),  SigThresh,SigSteepness))./tauPy;
    dInhdt    =(-Inh    + Sigm(parameters.Py2Inh*Py                  + parameters.InhInput,                SigThresh,SigSteepness))./tauInh;

    Y(:,i+1) = Y(:,i) + h*[dPydt;dInhdt];
  
  
end
Y = single(Y)';

end