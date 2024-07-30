function [X,Quad,face_front] = gen_fbox_qsph(nbx)
% Generate surface quad mesh for cubed sphere
% radius = 1
fprintf('CubeSph: gen Quad on sph (%d x %d / face)... ',nbx,nbx); t0 = tic;
nelf = nbx*nbx;
nptf = (nbx+1)*(nbx+1);

X = zeros(nptf*6,3);
Quad = zeros(nelf*6,4);

v=sqrt(3)/3;

% lexico graphical vertices
xv = [-v,v,-v,v,-v,v,-v,v];
yv = [-v,-v,v,v,-v,-v,v,v];
zv = [-v,-v,-v,-v,v,v,v,v];
findx = [1,2,6,5; 2,4,8,6; 4,3,7,8; 3,1,5,7; 3,4,2,1; 5,6,8,7];

nX=0;
for f=1:6
   vids=findx(f,:);
   [Xnew,Qnew] = gen_quad_on_surface(xv(vids),yv(vids),zv(vids),nbx);

   iq0 = (f-1)*nelf+1;
   iq1 = f*nelf;
   Quad(iq0:iq1,:)=Qnew + nX;

   ix0 = (f-1)*nptf+1;
   ix1 = f*nptf;
   X(ix0:ix1,:) = Xnew;
   nX = nX + nptf;
end

face_front = zeros(nelf,6);
face_front(:,1) = (1:nelf);
face_front(:,2) = (1:nelf) + 1*nelf;
face_front(:,3) = (1:nelf) + 2*nelf;
face_front(:,4) = (1:nelf) + 3*nelf;
face_front(:,5) = (1:nelf) + 4*nelf;
face_front(:,6) = (1:nelf) + 5*nelf;

fprintf(' done! nQ=%d (%2.4e sec)\n',size(Quad,1),toc(t0));

function [Xnew,Qnew]=gen_quad_on_surface(xv,yv,zv,nbx)
v1 = [xv(1),yv(1),zv(1)];
v2 = [xv(2),yv(2),zv(2)];
v4 = [xv(4),yv(4),zv(4)];
[Xnew,Qnew]=gen_fbox_geom(v1,v2,v4,nbx,2,1);

