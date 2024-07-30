function draw_quad_vtk(X,QuadsAll,Qdata,fldr,str,imode)

% imode:
%   =0 default (ascii)
%    1 (ascii) plot all facets 
%   -1 (binary) plot all facets 

t0=tic; nQ=size(QuadsAll,1); Qdata=Qdata(:);

fname = 'QuadsAll.vtk'; if(~isempty(str));fname=['QuadsAll_' str '.vtk'];end; 
if(~isempty(fldr)); fname=[fldr '/' fname]; end
format='ascii'; if(imode<0);format='binary';end

if(isinf(X(end,1))==1); X=X(1:end-1,:);end
nX=size(X,1); 

if (abs(imode)==1 || imode==0);
  idq=1:nQ; nq=nQ;
else
  error('un-support imode in draw_Ufacets_vtk');
end

vtk_title='Quads'; vtk_title = vtk_title(1:min(length(vtk_title),256));

switch format
  case 'ascii'
    fid=fopen(fname,'wt');
    assert(fid>0,['fail to open file: ' fname]);
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,[vtk_title '\n']);
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'\n');

    fprintf(fid,'POINTS %d float\n',nX);
    s='%f %f %f \n'; fprintf(fid,s,X');fprintf(fid,'\n');

    fprintf(fid,'CELLS %d %d\n',nq,nq*5);
    dat=[4*ones(nq,1),QuadsAll(idq,:)-1];
    fprintf(fid,'%d %d %d %d %d \n',dat');
    fprintf(fid,'\n');

    fprintf(fid,'CELL_TYPES %d\n',nq);
    fprintf(fid,'%d\n',7*ones(nq,1));
    fprintf(fid,'\n');

    fprintf(fid,'CELL_DATA %d\n',nq);
    fprintf(fid,'SCALARS cell_id int 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%d\n',1:nQ);
    fprintf(fid,'SCALARS data int 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%d\n',Qdata);
    fprintf(fid,'\n');

  case 'binary'
    write_a = @(fid,str) [fprintf(fid,str); fprintf(fid,'\n');];
    write_b = @(fid,dat,prec) [fwrite(fid,dat,prec); fprintf(fid,'\n');];
    write_b2= @(fid,dat,prec) [fwrite(fid,dat,prec)];

    fid=fopen(fname,'wb','ieee-be');
    assert(fid>0,['fail to open file: ' fname]);

    write_a(fid,'# vtk DataFile Version 2.0');

    write_a(fid,vtk_title);
    write_a(fid,'BINARY\n');
    write_a(fid,'DATASET UNSTRUCTURED_GRID');

    write_a(fid,['POINTS ' num2str(nX) ' float']);
    write_b(fid,X','float32');
    write_a(fid,'');
          
    write_a(fid,['CELLS ' num2str(nQ) ' ' num2str(5*nQ)]);
    dat = uint32([4*ones(nq,1),QuadsAll(idq,:)-1]);
    write_b2(fid,dat','uint32'); % no break line
    write_a(fid,''); write_a(fid,'');

    write_a(fid,['CELL_TYPES ' num2str(nq)]);
    write_b(fid,uint32(7*ones(1,nq)),'uint32');
    write_a(fid,'');

    write_a(fid,['CELL_DATA ' num2str(nq)]);
    write_a(fid,'SCALARS cell_id int 1');
    write_a(fid,'LOOKUP_TABLE default');
    write_b2(fid,uint32(1:nQ),'uint32');
    write_a(fid,'SCALARS data int 1');
    write_a(fid,'LOOKUP_TABLE default');
    write_b2(fid,uint32(Qdata),'uint32');
    write_a(fid,'');
    write_a(fid,'');

  otherwise
    error('wrong input dummy :P');
end


fclose(fid);
[osize,otype]=comp_fsize(fname);
fprintf(['    vtk: dump %d Quads into ' fname ' (%3.1f %s, %2.4e sec)\n'],nQ,osize,otype,toc(t0))
