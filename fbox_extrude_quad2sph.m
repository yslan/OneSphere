function [X,Hexes,Hface_front,Hbc,Hcurve] = fbox_extrude_quad2sphi...
         (X, Quad,Qface_front,nbx,R0,Rs,nlayers,ratio)

fprintf('Extrude: sph to sph, %d lyr',nlayers); t0 = tic;

% determine the thickness (in terms of ratio) per layer
wt = zeros(nlayers,1);
wt(1) = 1.0;
for i=2:nlayers
   wt(i) = ratio * wt(i-1);
end
wt = wt / sum(wt) * (Rs-R0);

fprintf(', radius: %2.2f',wt(1));
rad_new = zeros(nlayers,1);
rad_new(1) = R0 + wt(1);
for i=2:nlayers
   rad_new(i) = rad_new(i-1) + wt(i);
   fprintf(' %2.2f',rad_new(i));
end



nelf = nbx*nbx;
nptf = (nbx+1)*(nbx+1);

Xnew = zeros(6*nptf*nlayers,3);
Hnew = zeros(6*nelf*nlayers,8);
Hbc = zeros(6*nelf*nlayers,6);
Hcurve = zeros(6,12,nelf*nlayers);

nXnow = 0;
nHnow = 0;
nXtot = size(X,1);
for f=1:6

   Qface = Quad(Qface_front(:,f),:);

   v1 = X(Qface(1,1),:);
   v2 = X(Qface(nbx,2),:);
   v4 = X(Qface((nbx-1)*nbx+1,4),:);

   Qfab0 = Qface;
   rad0 = R0;
   [Xtmp,Qtmp] = gen_fbox_geom(v1,v2,v4,nbx,2,1);
   nXnew = size(Xtmp,1);
   nQnew = size(Qtmp,1);
   for i=1:nlayers
      rad = rad_new(i);
      rad1 = rad;

      Qfab1 = Qtmp + nXtot;

      Hnew((nHnow+1):(nHnow+nQnew),1:4) = Qfab0;
      Hnew((nHnow+1):(nHnow+nQnew),5:8) = Qfab1;

      if i==1
         Hbc((nHnow+1):(nHnow+nQnew),5) = f;
      elseif i==nlayers
         Hbc((nHnow+1):(nHnow+nQnew),6) = f+10;
      end
      Xnew((nXnow+1):(nXnow+nXnew),:) = Xtmp * rad;

      Hcurve(1,5,(nHnow+1):(nHnow+nQnew)) = 2;
      Hcurve(2:4,5,(nHnow+1):(nHnow+nQnew)) = 0;
      Hcurve(5,5,(nHnow+1):(nHnow+nQnew)) = rad0;
      Hcurve(1,6,(nHnow+1):(nHnow+nQnew)) = 2;
      Hcurve(2:4,6,(nHnow+1):(nHnow+nQnew)) = 0;
      Hcurve(5,6,(nHnow+1):(nHnow+nQnew)) = rad1;

      nHnow = nHnow + nQnew;
      nXnow = nXnow + nXnew;
      nXtot = nXtot + nXnew;
      Qfab0 = Qfab1;
      rad0 = rad1;
   end

   Hface_front(:,f) = (nHnow-nQnew+1):nHnow;
end

Hexes = Hnew;
X = [X;Xnew];
fprintf(' done! nH=%d (%2.4e sec)\n',size(Hexes,1),toc(t0));
