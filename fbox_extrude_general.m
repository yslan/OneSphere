function [X,Hexes,Hfront,Hbc,Hcurve] = fbox_extrude_general...
         (X,Hexes,Hfront,Hbc,Hcurve,nbx,geom1,geom2,nlayers,ratio)

s1 = gtype_to_str(geom1.type);
s2 = gtype_to_str(geom2.type);
fprintf('Extrude: %s to %s, %d lyr',s1,s2,nlayers); t0 = tic;

nface = 6;
if (geom2.type==3 || geom2.type==3); nface=4; end% cylinder: only do side faces 1-4

% determine the thickness (in terms of ratio) per layer
if (length(ratio)==1) % geometric series
   wt = zeros(nlayers,1);
   wt(1) = 1.0;
   for i=2:nlayers
      wt(i) = ratio * wt(i-1);
   end
   wt = wt / sum(wt);
   
   fprintf(', interp: %2.2f',wt(1));
   wt_acc = zeros(nlayers,1);
   wt_acc(1) = wt(1);
   for i=2:nlayers
      wt_acc(i) = wt_acc(i-1) + wt(i);
      fprintf(' %2.2f',wt_acc(i));
   end
else % custom ratio
   wt_acc = ratio;
   assert(nlayers==length(wt_acc), 'fbox_extrude_general: invalid custom ratio');
   fprintf(', interp c: %2.2f',wt_acc(1));
   for i=2:nlayers
      assert(wt_acc(i) > wt_acc(i-1), 'fbox_extrude_general: invalid custom ratio');
      fprintf(' %2.2f',wt_acc(i));
   end
end
wt_acc(end) = 1.0;

nelf = nbx*nbx;
nptf = (nbx+1)*(nbx+1);

Xnew = zeros(nface*nptf*nlayers,3);
Hnew = zeros(nface*nelf*nlayers,8);
Hcurve_new = zeros(6,12,nface*nelf*nlayers);
Hcurve = cat(3,Hcurve,Hcurve_new);
f_to_eid_3d = [1,5;2,6;3,7;4,8;0,0;0,0];
fff1 = @(X) 1.0/sqrt(X(:,1).^2+X(:,2).^2);
fff2 = @(X,R) [X(:,1).*R, X(:,2).*R, X(:,3)]; 

nXnow = 0;
nHnow = 0;
nXtot = size(X,1);
nHtot = size(Hexes,1); nHtot0=nHtot;
bc_replace_list=[];
for f=1:nface

   Qface = Hexes(Hfront(:,f),5:8); % face 6
   bc_replace_list = [bc_replace_list,Hbc(Hfront(:,f),6)];

   v1 = X(Qface(1,1),:);
   v2 = X(Qface(nbx,2),:);
   v4 = X(Qface((nbx-1)*nbx+1,4),:);

   [Xbot,Qbot] = gen_fbox_geom(v1,v2,v4,nbx,geom1.type,1);
   Xbot = Xbot * geom1.coef;
%   v1 = v1/sqrt(sum(v1.^2))*sqrt(3);
%   v2 = v2/sqrt(sum(v2.^2))*sqrt(3);
%   v4 = v4/sqrt(sum(v4.^2))*sqrt(3);
   [Xtop,Qtop] = gen_fbox_geom(v1,v2,v4,nbx,geom2.type,1);
   if (geom2.type==3)
      Xtop(:,1) = Xtop(:,1) * geom2.coef;
      Xtop(:,2) = Xtop(:,2) * geom2.coef;
      Xtop(:,3) = Xtop(:,3) * geom1.coef;
   else
      Xtop = Xtop * geom2.coef;
   end

   Qfab0 = Qface;
   nXnew = size(Xbot,1);
   nQnew = size(Qbot,1);
   for i=1:nlayers
      Xtmp = Xtop * wt_acc(i) + Xbot * (1-wt_acc(i));
      Qfab1 = Qtop + nXtot;

      Hnew((nHnow+1):(nHnow+nQnew),1:4) = Qfab0;
      Hnew((nHnow+1):(nHnow+nQnew),5:8) = Qfab1;

      if (i==1) % face 5
         if (geom1.type==2) 
            Hcurve(1,5,(nHnow+1):(nHnow+nQnew)) = 2;
            Hcurve(2:4,5,(nHnow+1):(nHnow+nQnew)) = 0;
            Hcurve(5,5,(nHnow+1):(nHnow+nQnew)) = geom1.coef;
         elseif (geom1.type==3) 
            % hollow cyl (outside)
            %Hcurve(1:2,f_to_eid_3d(5,1),(nHnow+1):(nHnow+nQnew)) = [3,-geom1.coef];
            %Hcurve(1:2,f_to_eid_3d(5,2),(nHnow+1):(nHnow+nQnew)) = [3,-geom1.coef];
            % midpt
            ex1 = (Xbot(Qtop(:,1),:)+Xbot(Qtop(:,2),:))/2; rr=fff1(ex1);ex1=fff2(ex1,rr); 
            ex2 = (Xbot(Qtop(:,2),:)+Xbot(Qtop(:,3),:))/2; rr=fff1(ex2);ex2=fff2(ex2,rr);
            ex3 = (Xbot(Qtop(:,3),:)+Xbot(Qtop(:,4),:))/2; rr=fff1(ex3);ex3=fff2(ex3,rr);
            ex4 = (Xbot(Qtop(:,4),:)+Xbot(Qtop(:,1),:))/2; rr=fff1(ex4);ex4=fff2(ex4,rr);
            Hcurve(1,[1,2,3,4],(nHnow+1):(nHnow+nQnew)) = 1;
            Hcurve(2:4,1,(nHnow+1):(nHnow+nQnew)) = ex1;
            Hcurve(2:4,2,(nHnow+1):(nHnow+nQnew)) = ex2;
            Hcurve(2:4,3,(nHnow+1):(nHnow+nQnew)) = ex3;
            Hcurve(2:4,4,(nHnow+1):(nHnow+nQnew)) = ex4;
         end
      end
      if (i==nlayers) % face 6
         if (geom2.type==2) 
            Hcurve(1,6,(nHnow+1):(nHnow+nQnew)) = 2;
            Hcurve(2:4,6,(nHnow+1):(nHnow+nQnew)) = 0;
            Hcurve(5,6,(nHnow+1):(nHnow+nQnew)) = geom2.coef;
         elseif (geom1.type==3) 
            %Hcurve(1:2,f_to_eid_3d(6,1),(nHnow+1):(nHnow+nQnew)) = [3,-geom1.coef];
            %Hcurve(1:2,f_to_eid_3d(6,2),(nHnow+1):(nHnow+nQnew)) = [3,-geom1.coef];
            % midpt
            ex1 = (Xtop(Qtop(:,1),:)+Xtop(Qtop(:,2),:))/2; rr=fff1(ex1);ex1=fff2(ex1,rr);
            ex2 = (Xtop(Qtop(:,2),:)+Xtop(Qtop(:,3),:))/2; rr=fff1(ex2);ex2=fff2(ex2,rr);
            ex3 = (Xtop(Qtop(:,3),:)+Xtop(Qtop(:,4),:))/2; rr=fff1(ex3);ex3=fff2(ex3,rr);
            ex4 = (Xtop(Qtop(:,4),:)+Xtop(Qtop(:,1),:))/2; rr=fff1(ex4);ex4=fff2(ex4,rr);
            Hcurve(1,[5,6,7,8],(nHnow+1):(nHnow+nQnew)) = 1;
            Hcurve(2:4,5,(nHnow+1):(nHnow+nQnew)) = ex1;
            Hcurve(2:4,6,(nHnow+1):(nHnow+nQnew)) = ex2;
            Hcurve(2:4,7,(nHnow+1):(nHnow+nQnew)) = ex3;
            Hcurve(2:4,8,(nHnow+1):(nHnow+nQnew)) = ex4;
         end
      end

      Xnew((nXnow+1):(nXnow+nXnew),:) = Xtmp;

      nHnow = nHnow + nQnew;
      nHtot = nHtot + nQnew;
      nXnow = nXnow + nXnew;
      nXtot = nXtot + nXnew;
      Qfab0 = Qfab1;
   end

   Hfront(:,f) = (nHtot-nQnew+1):nHtot;
end

Hexes=[Hexes;Hnew];
X = [X;Xnew];
bc_shift = max(abs(Hbc(:))); bc_shift = (floor(bc_shift/10)+1)*10;
Hbc = [Hbc;zeros(nHnow,6)];
for f=1:nface
   Hbc(Hfront(:,f),6) = f + bc_shift;
end 

% flip BCid that is no longer boundary
bc_replace_list = unique(bc_replace_list);
for i = 1:length(bc_replace_list)
   bc = bc_replace_list(i);
   id_bc = find(Hbc==bc);
   Hbc(id_bc) = -abs(Hbc(id_bc));
end

% mesh is not water tight so it relies on tol to set top/bot BC for cylinder
iftoiv=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
tol=1e-6;
Zming = min(Xnew(:,3));
Zmaxg = max(Xnew(:,3));
if (geom2.type==3)
   for e=(nHtot0+1):nHtot
   for f=1:4
      quad = Hexes(e,iftoiv(f,:));
      zmin = min(X(quad,3));
      zmax = max(X(quad,3));
      if (zmax-zmin<tol);
         zmid = (zmin+zmax)/2;
         if (abs(zmid-Zming)<tol); Hbc(e,f) = 5 + bc_shift; end
         if (abs(zmid-Zmaxg)<tol); Hbc(e,f) = 6 + bc_shift; end
      end
   end
   end
end

fprintf(' done! nH=%d (%2.4e sec)\n',size(Hexes,1),toc(t0));

function s=gtype_to_str(gtype)
switch gtype
   case 1
      s='box';
   case 2
      s='sph';
   case 3
      s='cyl';
   otherwise
      error('fbox_extrude_general: invalid gtype')
end
