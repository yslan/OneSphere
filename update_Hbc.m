function Hbc = update_Hbc(Hbc,Hbc_new);

bc_shift = max(abs(Hbc(:))); 
bc_shift = (floor(bc_shift/10)+1)*10;
Hbc_new(Hbc_new~=0)=Hbc_new(Hbc_new~=0)+bc_shift; 
Hbc=[Hbc;Hbc_new];

