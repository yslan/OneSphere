function [Xnew,Hnew,Hbcnew,bout] = genbox(dat,verbose)

nelx=dat.nelx; xx=dat.xx;
nely=dat.nely; yy=dat.yy;
nelz=dat.nelz; zz=dat.zz;
fprintf('Genbox : %3d x %3d x %3d',abs(nelx),abs(nely),abs(nelz)); t0=tic;

chk_input('x',nelx,xx);
chk_input('y',nely,yy);
chk_input('z',nelz,zz);

xx = get_endpoints(nelx,xx,verbose-2);
yy = get_endpoints(nely,yy,verbose-2);
zz = get_endpoints(nelz,zz,verbose-2);

nelx = abs(nelx);
nely = abs(nely);
nelz = abs(nelz);
nelxy = nelx*nely;

bout.nelx = nelx;
bout.nely = nely;
bout.nelz = nelz;
bout.xx=xx;
bout.yy=yy;
bout.zz=zz;

nptx = nelx+1;
npty = nely+1;
nptz = nelz+1;
nptxy = nptx*npty;

[X,Y,Z] = ndgrid(xx,yy,zz);
Xnew = [X(:),Y(:),Z(:)];

Hnew = zeros(nelx*nely*nelz,8);
Hbcnew = zeros(nelx*nely*nelz,6);
for ez=1:nelz
for ey=1:nely
for ex=1:nelx
   e = (ez-1)*nelxy + (ey-1)*nelx + ex;
   iv = (ez-1)*nptxy + (ey-1)*nptx + ex;
   face5 = [iv,iv+1,iv+1+nptx,iv+nptx];
   face6 = face5 + nptxy;
   Hnew(e,:) = [face5,face6];

   if (ey==1);    Hbcnew(e,1) = 1; end
   if (ex==nelx); Hbcnew(e,2) = 2; end
   if (ey==nely); Hbcnew(e,3) = 3; end
   if (ex==1);    Hbcnew(e,4) = 4; end
   if (ez==1);    Hbcnew(e,5) = 5; end
   if (ez==nelz); Hbcnew(e,6) = 6; end
end
end
end

fprintf(' done!  (%2.4e sec)\n',toc(t0));

function chk_input(str,nelx,xx)
if nelx>0
   assert(length(xx)==nelx+1,['genbox: ' str ' size mismatched'])
elseif nelx<0
   assert(length(xx)==3,['genbox: ' str ' size mismatched'])
   x0=xx(1); x1=xx(2);
   assert(x1>x0, ['genbox: ' str ' x1>x0 is required'])
else
   error('genbox: nel=0 detected')
end

function xout = get_endpoints(nelx,xx,verbose)
if nelx>0
   xout = xx;
else
   nelx=abs(nelx);
   
   x0 = xx(1); x1 = xx(2); ratio = xx(3);

   del = zeros(nelx,1);
   del(1) = 1.0;
   for e=2:nelx
      del(e) = del(e-1) * ratio;
   end
   del = del * (x1-x0) / sum(del);

   xout = zeros(nelx+1,1);

   xout(1) = x0; 
   xout(2) = xout(1) + del(1);
   for e=2:nelx
      xout(e+1) = xout(e) + del(e);
   end
end

if (verbose>0)
   fprintf('\n  genbox endpt:');
   for i=1:length(xout)
      fprintf(' %f',xout(i));
   end
   fprintf('\n');
end

