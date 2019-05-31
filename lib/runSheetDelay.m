function Y=runSheetDelay(y0,parameters)

%this is an fixed step solver with a noise term as subcortical input and
%delay


nsq=parameters.n^2;
h=parameters.h;
ds=parameters.delaysteps;

SigThresh=parameters.SigThresh;
SigSteepness=parameters.SigSteepness;
tauPy=parameters.tauPy;
tauInh=parameters.tauInh;


neq = 2*parameters.n^2;
Nl = size(parameters.NValue,2);
Y = zeros(neq,Nl+ds);

Y(:,1:ds) = y0;
for i = ds:ds+Nl-1 
    
    Py=Y(1:nsq,i-ds+1:i);
    Inh=Y(nsq+1:2*nsq,i);
    
    PyPyInp=zeros(size(Py,1),1);
    PyInhInp=zeros(size(Inh,1),1);
    for k=1:ds
        dsk=ds-k+1;
        PyPyInp=PyPyInp+parameters.Py2Py{k}*Py(:,dsk);
        if k<=length(parameters.Py2Inh)
            PyInhInp=PyInhInp+parameters.Py2Inh{k}*Py(:,dsk);
        end
    end
    
    
    
    dPydt     =(-Py(:,end)     + Sigm(PyPyInp      - parameters.Inh2Py*Inh + parameters.PyInput + parameters.NValue(:,i-ds+1),  SigThresh,SigSteepness))./tauPy;
    dInhdt    =(-Inh    + Sigm(PyInhInp                  + parameters.InhInput,                SigThresh,SigSteepness))./tauInh;

    Y(:,i+1) = Y(:,i) + h*[dPydt;dInhdt];
  
  
end
Y = single(Y)';

end