function sout=Sigm(x,thresh,steepness)
    sout=1./(1+exp(-steepness.*(x-thresh)));
end