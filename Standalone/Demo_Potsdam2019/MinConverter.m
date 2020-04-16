%% MinConverter
% Mineral Converter 1.0 (Last modification 13.05.19 - ED)
%
% converts oxides [wt.%] to atoms [mol]

function [WorkVariXMap] = MinConverter(WorkVariXMap)

Converter.OxyNames={'SI','TI','AL','F3','FE','MN','MG','CA','NA','K'};
Converter.OxyMass=[60.085,79.9,101.963,159.68,71.84,70.938,40.305,56.08,61.97954,94.196];
Converter.NbCat=[1,1,2,2,1,1,1,1,2,2];
Converter.NbOx=[2,2,3,3,1,1,1,1,1,1];


WorkVariXMap.COMP=zeros(size(WorkVariXMap.COMPoxy));
summe_oxy=zeros(WorkVariXMap.NbPhases,1);
for i=1:WorkVariXMap.NbPhases
    if ismember(WorkVariXMap.Format(i),'O')
        for j=2:WorkVariXMap.NbEl-2
            [IDX.IsIn,IDX.Indices] = ismember(WorkVariXMap.Els(j),Converter.OxyNames);
            WorkVariXMap.COMP(i,j)=Converter.NbCat(IDX.Indices).*WorkVariXMap.COMPoxy(i,j)./Converter.OxyMass(IDX.Indices);
            summe_oxy(i)=summe_oxy(i) + Converter.NbOx(IDX.Indices).*WorkVariXMap.COMPoxy(i,j)./Converter.OxyMass(IDX.Indices);
        end 
        WorkVariXMap.COMP(i,:)=WorkVariXMap.COMP(i,:)./(summe_oxy(i)./WorkVariXMap.COMPoxy(i,1));
    else
        WorkVariXMap.COMP(i,:)=WorkVariXMap.COMPoxy(i,:);
    end
end
WorkVariXMap.COMP(:,1)=WorkVariXMap.COMPoxy(:,1);
%%
end