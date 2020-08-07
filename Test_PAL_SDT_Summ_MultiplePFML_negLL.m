%
%SDT_Summ_MultiplePFML_negLL     (negative) Log Likelihood associated with fit 
%   of summation psychometric function 
%
%Internal Function
%
%Introduced: Palamedes version 1.8.0 (FK & NP)

function negLL = Test_SDT_Summ_MultiplePFML_negLL(paramsIn,StimLevels,NumPos,OutOfNum,SummModel,M,Q,CM)

[nrows, ncols, nnums]=size(StimLevels);
PC=zeros(nrows,nnums);
StimVec=zeros(1,ncols); % vector for corresponding stimulus levels for summation equation

params=paramsIn*CM;
num=length(params);
gParams=params(1:num/2);
pParams=params(num/2+1:num);


% Calculate PCs for the components 
if strcmp(SummModel,'PS')==1 
    PC(1,:)=PAL_SDT_PS_SLtoPC(StimLevels(1,1,:),gParams(1),pParams(1),M,Q(1),1);
    PC(nrows,:)=PAL_SDT_PS_SLtoPC(StimLevels(nrows,2,:),gParams(2),pParams(2),M,Q(1),1);
    SummFuncCompound=@PAL_SDT_PS_uneqSLtoPC;
end
if strcmp(SummModel,'AS')==1 
    PC(1,:)=PAL_SDT_AS_SLtoPC(StimLevels(1,1,:),gParams(1),pParams(1),M,Q(1),1);
    PC(nrows,:)=PAL_SDT_AS_SLtoPC(StimLevels(nrows,2,:),gParams(2),pParams(2),M,Q(1),1);
    SummFuncCompound=@PAL_SDT_AS_uneqSLtoPC;
end

for j=2:nrows-1
    for i=1:nnums
        for k=1:ncols
            StimVec(1,k)=StimLevels(j,k,i);
        end
        PC(j,i) = SummFuncCompound(StimVec,gParams,pParams,M,Q(2));
    end
end


negLL = -sum(PAL_nansum(NumPos.*log(PC))+PAL_nansum((OutOfNum-NumPos).*log(1 - PC)));