function [Xnew,Qnew]=gen_fbox_geom(v1,v2,v4,nbx,gtype,ifproj)

% equi-space
wt_x = (0:nbx) / nbx;

% equi-angle
theta = pi/2 / nbx;
r = @(i) sin(i* theta) ./ sin(3*pi/4-i*theta);
wt_theta = r(0:nbx) / sqrt(2);

switch gtype
   case 1 % plane
      rwt = wt_x; swt = wt_x;
   case 2 % sph
      rwt = wt_theta; swt = wt_theta;
   case 3 % cyl
      rwt = wt_theta;
      swt = wt_x; % always assume y = z
% TODO
   otherwise
      error('gen_fbox_qsph: unsupport geom');
end
rwt = reshape(rwt,nbx+1,1); rwt(1) = 0.0; rwt(end) = 1.0;
swt = reshape(swt,nbx+1,1); swt(1) = 0.0; swt(end) = 1.0;

%v1 = [xv(1),yv(1),zv(1)];
%v2 = [xv(2),yv(2),zv(2)];
%v4 = [xv(4),yv(4),zv(4)];
rd = v2-v1;
sd = v4-v1;

Xnew = zeros((nbx+1)*(nbx+1),3);
for j=1:(nbx+1)
for i=1:(nbx+1)
   ix = (j-1)*(nbx+1)+i;
   Xnew(ix,:) = v1 + rd*rwt(i) + sd*swt(j);
end
end
%[ii,jj] = ndgrid(1:(nbx+1),1:(nbx+1)); ii=ii(:); jj=jj(:);
%Xnew = [v1(1) + rd(1)*wt(ii) + sd(1)*wt(jj),...
%        v1(2) + rd(2)*wt(ii) + sd(2)*wt(jj),...
%        v1(3) + rd(3)*wt(ii) + sd(3)*wt(jj)];

if(ifproj>0) 
   switch gtype
      case 1
         Xnew = Xnew / max(abs(Xnew(:)));
      case 2
         d = 1.0 ./ sqrt(sum(Xnew.^2,2));
         Xnew(:,1) = Xnew(:,1) .* d;
         Xnew(:,2) = Xnew(:,2) .* d;
         Xnew(:,3) = Xnew(:,3) .* d;
      case 3
         r = 1.0 ./ sqrt(Xnew(:,1).^2+Xnew(:,2).^2);
         Xnew(:,1) = Xnew(:,1) .* r;
         Xnew(:,2) = Xnew(:,2) .* r;
         h2 = max(abs(Xnew(:,3)));
         Xnew(:,3) = Xnew(:,3) / h2;
   end
end

Qnew = zeros(nbx*nbx,4);
iq=0;
for j=1:nbx
for i=1:nbx
   iv = (j-1)*(nbx+1)+i;
   iq = iq+1;
   Qnew(iq,:) = [iv,iv+1,iv+1+nbx+1,iv+nbx+1];
end
end
%[ii,jj] = ndgrid(1:nbx,1:nbx);
%iv = (jj(:)-1)*(nbx+1) + ii(:);
%Qnew = [iv,iv+1,iv+1+nbx+1,iv+nbx+1];



