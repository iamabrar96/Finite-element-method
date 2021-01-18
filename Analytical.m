function Analyticalsol=Analytical(v,p,E,r1,r2,rnodes)
Analyticalsol=zeros(length(rnodes),1);
Constant=(1+v)*(p/E)*((r1*r1)/((r2*r2)-(r1*r1)));
for i=1:length(rnodes)
    Analyticalsol(i)=Constant*(((1-2*v)*rnodes(i))+((r2*r2)/rnodes(i)));
end
end
