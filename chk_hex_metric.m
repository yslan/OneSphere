function dat=chk_hex_metric(X,Hexes,str,verbose,varargin)

t0=tic; 
dim = size(X,2); nv=2^dim;
E=size(Hexes,1); N=8;%distribution

if(length(varargin)>0);N=varargin{1};end
for d=1:dim; Xtmp=X(Hexes,d); Xl(:,:,d)=reshape(Xtmp,E,nv); end % GtoL 

mma=@(v)[min(v(:)),max(v(:)),sum(v(:))/numel(v)]; % [min,max,ave]

% spacing
if (dim==3)
   iv1=[1,1,1]; iv2=[2,4,5];
else
   iv1=[1,1]; iv2=[2,4];
end
Xdif=Xl(:,iv1,:)-Xl(:,iv2,:); Elen=sqrt(sum(Xdif.^2,3)); dxmin=min(Elen,[],2);
dat.dxmin=mma(dxmin);

% scaled Jac
jacm=comp_Jacobian_v3_2(Xl); jacm=reshape(jacm,nv,E)'; sc_jac=min(jacm,[],2)./max(jacm,[],2);
dat.sc_jac=mma(sc_jac);

% aspect ratio
if (dim==3); % 12 edges
   iv1=[1,2,3,4,5,6,7,8,1,2,3,4];
   iv2=[2,3,4,1,6,7,8,5,5,6,7,8];
else; % 4 edges
   iv1=[1,2,3,4];
   iv2=[2,3,4,1];
end
Xdif=Xl(:,iv1,:)-Xl(:,iv2,:); Elen=sqrt(sum(Xdif.^2,3)); aratio=max(Elen,[],2)./min(Elen,[],2);
dat.aratio=mma(aratio);

% max multiplicity
dat.max_mult=max(accumarray(Hexes(:),1)); % only works for water tight mesh.

% check and print
if (verbose>0)
if (sum(isnan(sc_jac(:)))>0); % nanJac
  nanJac=max(isnan(sc_jac),[],2); nanJac=find(nanJac>0); nerr=length(nanJac);
  nprt=min(3,nerr);eid = nanJac(1:nprt);
  fprintf('WARN: Found %d nanJac elements! eid=',nerr); for i=1:nprt;fprintf(' %3d',eid(i));end;fprintf('\n');
end
if (dat.sc_jac(1)<=0); % negJac
  minJac=min(jacm,[],2); negJac=find(minJac<0); nerr=length(negJac); 
  nprt=min(3,nerr);eid = negJac(1:nprt);
  fprintf('WARN: Found %d negJac elements! eid=',nerr); for i=1:nprt;fprintf(' %3d',eid(i));end;fprintf('\n');
end
if (min(dxmin)<=0); % zero / negative edges
  badEdge=find(dxmin<=0); nerr=length(badEdge);
  nprt=min(3,nerr);eid = badEdge(1:nprt);
  fprintf('WARN: Found %d elements with badEdge (<=0)! eid=',nerr); for i=1:nprt;fprintf(' %3d',eid(i));end;fprintf('\n');
end

ifhist=0; if(exist('histcounts')==2); ifhist=1; end
if(ifhist==1);
   [c1,e1,n1]=get_bin(dxmin,N); c1=int8(c1/E*100);
   [c2,e2,n2]=get_bin(sc_jac,N);c2=int8(c2/E*100);
   [c3,e3,n3]=get_bin(aratio,N);c3=int8(c3/E*100);
else
   N=0;
end

stype='Hex'; if (dim==2); stype='Quad'; end
fprintf('  %s metrics: %s (min/max/ave | distr) (%2.4e sec)\n',stype,str,toc(t0));
fprintf('    GLL grid spacing %9.2e %2.2e %2.2e |',dat.dxmin); for i=1:n1-1;fprintf(' %3d',c1(i));end;fprintf('\n');
fprintf('    scaled Jacobian  %9.2e %2.2e %2.2e |',dat.sc_jac);for i=1:n2-1;fprintf(' %3d',c2(i));end;fprintf('\n');
fprintf('    aspect ratio     %9.2e %2.2e %2.2e |',dat.aratio);for i=1:n3-1;fprintf(' %3d',c3(i));end;fprintf('\n');
fprintf('    max multiplicity %d\n',dat.max_mult);

end

if (verbose>1)
  dat.jacm = jacm;
end

function [c,e,n] = get_bin(v,N)
  assert(N>0,'get_bin');
  n = N;
  if sum(isnan(v(:)))>0 
    N = 0;
    fprintf('WARN: detect NaN\n');
    c=1;e=[0,1];n=0;
  else 
    [c,e] = histcounts(v,linspace(min(v(:)),max(v(:))+eps,N));
  end


