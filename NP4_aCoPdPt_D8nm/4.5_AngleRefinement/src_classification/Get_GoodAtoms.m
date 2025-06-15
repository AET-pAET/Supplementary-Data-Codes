function good_atom_pos = Get_GoodAtoms(Data_tracing,N,ourputstring_prefix,saveflag)
CoeffArr=Data_tracing.CoeffArr(:,1:end-N);
Orders=Data_tracing.Orders;
PosArr=Data_tracing.PosArr(1:end-N,:);
SearchRad=Data_tracing.SearchRad;
TotPosArr=Data_tracing.TotPosArr(1:end-N,:);
exitFlagArr=Data_tracing.exitFlagArr(:,1:end-N);
maxXYZ=Data_tracing.maxXYZ(1:end-N,:);

good_atom_pos = TotPosArr(exitFlagArr==0,:); %

if saveflag==1
save(sprintf('%s_result.mat',ourputstring_prefix),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr','SearchRad','maxXYZ');
end
end
    