function ClearSignal=PeakRemover(file,n)
[z,zfs]=audioread(file);
ThreeP = AutoPeak(file);
m=MidFinder(ThreeP,n);

m2=round(m(2)*zfs);
mend=round(m(end)*zfs);

ClearSignal=[z(1:mend); zeros((m2-mend),1); z(m2:end-1)];
end