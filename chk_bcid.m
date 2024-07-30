function set_new = bcid(set_old,Hbc,tag,imode);

set_new = unique(Hbc(Hbc~=0));
set_diff = setdiff(set_new,set_old);

if imode==1
   fprintf('  %s diff BC id:',tag);
   for i = set_diff
      fprintf(' %d',i);
   end
   fprintf('\n');
end


