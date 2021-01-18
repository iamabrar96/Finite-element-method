% Generate list of position of nodes according to a geometric series
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2020
%    lecturer in charge: Dr. Geralf H�tter
%
function rnodes=meshGenerator(a,b,nelem)
    if nelem==1
        rnodes=[a;b];
    else
        nelem=10;               %number of elements
        meshrefinementfactor=2; %ratio of element sizes at outer and inner radius

%ratio between element sizes of subsequent elements for a geometric series
        q=meshrefinementfactor^(1./(nelem-1));
%size of first interval
        dr=(b-a)*(1-q)/(1-meshrefinementfactor*q);
        rnode=a;
        rnodes=[a];
%loop over all elements
        for i=1:nelem
            rnode=rnode+dr;
            rnodes=[rnodes;rnode];
            dr=dr*q;
        end
    end
end


