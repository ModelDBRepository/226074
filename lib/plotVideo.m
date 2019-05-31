function plotVideo(n,T,mTS,TS,mycmap,caxisBounds,titles,outputFN,varargin)

%needed: T, mPy, mLFP, 
% n=100;
% mTS=mLFP;
% TS=LFP;
% mycmap;
%load('MayColourMap')
writerObj = VideoWriter(outputFN);
nVarargs = length(varargin);
if nVarargs>=1
    writerObj.FrameRate = varargin{1};
else
    writerObj.FrameRate = 50;
end
open(writerObj);
for k=1:length(T)
    hFig=figure(1);
    set(hFig, 'Position', [100 300 700 650]);
    subplot('position',[0.2 0.1 0.6 0.15])%[left bottom width height]
    minTS=min(mTS);
    maxTS=max(mTS);
    ampl=abs(minTS-maxTS);
    plot(T,mTS,'-k','Linewidth',2)
    hold on
    plot([T(k) T(k)],[minTS-0.1*ampl maxTS+0.1*ampl],'-r')
    hold off
    title(titles.TS)
    xlabel('Time [sec]')
    axis([T(1) T(end) minTS-0.1*ampl maxTS+0.1*ampl])
    
    subplot('position',[0.2 0.35 0.7 0.6])
    sidel=50e-3:50e-3:150*50e-3;
    imagesc(sidel,sidel,reshape(TS(k,:),n,n))
    colormap(mycmap)
    colorbar
    caxis(caxisBounds)
    xlabel('[mm]')
    ylabel('[mm]')
    title(titles.image)
    
    frame = getframe(hFig);
    writeVideo(writerObj,frame);
end
close(writerObj);

end