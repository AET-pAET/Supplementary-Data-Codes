function [GoodRings_No,Atom_type_Ring,Atom_pos_Ring] = Cal_GoodRings(cycles_sort_ref,cycles_num_ref,Pos_all,Type_all,n)

for j=1:length(n)
ID_ring=find(cycles_num_ref==n(j));
Num_tot=length(ID_ring);
count=0;
for i=1:Num_tot
ID_temp=ID_ring(i);
Atom_ID_Ring=cycles_sort_ref{ID_temp};
Atom_pos_Ring=Pos_all(:,Atom_ID_Ring);
Atom_type_Ring=Type_all(:,Atom_ID_Ring);

Atom_type_Ref=zeros(1,n(j));
Atom_type=unique(Atom_type_Ring);

if length(Atom_type)==2

Atom_type_Ref(1,1:2:n(j)-1)=Atom_type(1);
Atom_type_Ref(1,2:2:n(j))=Atom_type(2);
else
    Atom_type_Ref(1,1:2:n(j)-1)=Atom_type(1);
Atom_type_Ref(1,2:2:n(j))=0;
end
%Atom_type_Ring
%Atom_type_Ref

if Atom_type_Ring==Atom_type_Ref
    count=count+1;
    GoodRings(ID_temp)=1;
else
    GoodRings(ID_temp)=0;
end
end
GoodRings_No(j,:)=[n(j),count,count/Num_tot];
end

end

    