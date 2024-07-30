function CBC = set_cbc(Hbc,bc_cell);

CBC = zeros(size(Hbc,1),6);

nbcid = length(bc_cell);
for i=1:nbcid
   id_list = bc_cell{i};
   for j=1:length(id_list)
      is_bc = find(Hbc==id_list(j));
      CBC(is_bc) = i;
   end
end
