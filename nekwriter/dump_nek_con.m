function dump_nek_con(fname,Hexes,ifco2,verbose)

t0=tic; ext='.con'; if (ifco2); ext='.co2'; end; fname=[fname ext];
if (verbose>0); fprintf('Dump con: %s ... ',fname); end

if(ifco2); write_co2(fname,Hexes); 
else;      write_con(fname,Hexes); end

[osize,otype]=comp_fsize(fname);
if(verbose>0); fprintf(['   done! (%3.1f %s %2.4e sec)\n'],osize,otype,toc(t0)); end


function write_con(fname,Hexes)
  ivord = [1,2,4,3,5,6,8,7];
  [nH,nv]=size(Hexes); fmt='';for i=1:nv+1;fmt=[fmt '%12d'];end;fmt=[fmt '\n'];
  map=[(1:nH)',Hexes(:,ivord(1:nv))];

  fid=fopen(fname,'w');
  fprintf(fid,'#v001%12d%12d%12d\n',nH,nH,nv);
  fprintf(fid,fmt,map');
  fclose(fid);

function write_co2(fname,Hexes)
  ivord = [1,2,4,3,5,6,8,7];
  [nH,nv]=size(Hexes); map=[(1:nH)',Hexes(:,ivord(1:nv))];
  etag=654321; etag=etag*1e-5; emode = 'le';
  [fid,message] = fopen(fname,'w',['ieee-' emode]);
  if fid == -1, disp(message), status = -1; return, end

  header=sprintf('#v001%12d%12d%12d',nH,nH,nv);header(end+1:132) = ' ';
  fwrite(fid,header,'char');
  fwrite(fid,etag,'float32');
  fwrite(fid,map','int32');
  fclose(fid);


