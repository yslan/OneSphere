function [X,Hexes,Hbc] = extrude_onto_geom(X,Hexes,Hbc,bcid_list,func_geom_pj,nlayers,ratio)

fprintf('Extrude: geom to geom, %d lyr',nlayers); t0 = tic;

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

iftoiv=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;4 3 2 1;5 6 7 8];
iftofnei=[5 2 6 4;5 3 6 1;5 4 6 2;5 1 6 3;3 2 1 4;1 2 3 4];

Nhex0 = size(Hexes,1);
nbc = length(bcid_list);
elist=zeros(Nhex0,1);
flist=zeros(Nhex0,1);
neibclist=zeros(Nhex0,4);
Quad=zeros(Nhex0,4); iq0=1;iq1=0;
nQuad=zeros(nbc,1);
fprintf(', bcid:');
for i=1:nbc
   bc = bcid_list(i); fprintf(' %3d',bc);
   [ihex,ifac]=find(Hbc==bc);

   idf=sub2ind([Nhex0,6],ihex,ifac);
   id1=sub2ind([Nhex0,8],ihex,iftoiv(ifac,1));
   id2=sub2ind([Nhex0,8],ihex,iftoiv(ifac,2));
   id3=sub2ind([Nhex0,8],ihex,iftoiv(ifac,3));
   id4=sub2ind([Nhex0,8],ihex,iftoiv(ifac,4));

   idneif1=sub2ind([Nhex0,6],ihex,iftofnei(ifac,1));
   idneif2=sub2ind([Nhex0,6],ihex,iftofnei(ifac,2));
   idneif3=sub2ind([Nhex0,6],ihex,iftofnei(ifac,3));
   idneif4=sub2ind([Nhex0,6],ihex,iftofnei(ifac,4));

   Qtmp=[Hexes(id1),Hexes(id2),Hexes(id3),Hexes(id4)];nQtmp=size(Qtmp,1);
   iq1=iq1+nQtmp;
   Quad(iq0:iq1,:)=Qtmp;
   elist(iq0:iq1) = ihex;
   flist(iq0:iq1) = ifac;
   neibclist(iq0:iq1,:) = [Hbc(idneif1),Hbc(idneif2),Hbc(idneif3),Hbc(idneif4)];
   iq0=iq0+nQtmp;
   nQuad(i)=length(ihex);

   Hbc(idf)=-abs(bc);
end
Quad=Quad(1:iq1,:);elist=elist(1:iq1,:);flist=flist(1:iq1,:);neibclist=neibclist(1:iq1,:);

Xquad = X(Quad,:);
Xtarget = func_geom_pj(Xquad);




nXnow = 0;
nHnow = 0;
nXtot = size(X,1);
nHtot = size(Hexes,1); nHtot0=nHtot;
nQnew = size(Quad,1);

Quad_orig = Quad;

% only generates uniq points(well, this is extra...)
[ia,ib,ic] = unique(Quad);
Quad_uniq = reshape(ic,nQnew,4);
Xtarget = Xtarget(ib,:);
Xquad = Xquad(ib,:);

nXnew = size(Xtarget,1);
Xnew = zeros(nXnew*nlayers,3);
Hnew = zeros(nQnew*nlayers,8);
Hbcnew = zeros(nQnew*nlayers,6);
Qfab0 = Quad_orig;
bc_shift = max(abs(Hbc(:))); bc_shift = (floor(bc_shift/10)+1)*10;
for i=1:nlayers

   Xtmp = Xtarget*wt_acc(i) + Xquad*(1.0-wt_acc(i));
   Qfab1 = Quad_uniq + nXtot;

   Hnew((nHnow+1):(nHnow+nQnew),1:4) = Qfab0;
   Hnew((nHnow+1):(nHnow+nQnew),5:8) = Qfab1;

   Hbcnew((nHnow+1):(nHnow+nQnew),1) = abs(neibclist(:,1));
   Hbcnew((nHnow+1):(nHnow+nQnew),2) = abs(neibclist(:,2));
   Hbcnew((nHnow+1):(nHnow+nQnew),3) = abs(neibclist(:,3));
   Hbcnew((nHnow+1):(nHnow+nQnew),4) = abs(neibclist(:,4));
   if (i==nlayers)
      Hbcnew((nHnow+1):(nHnow+nQnew),6) = bc_shift + 6;
   end

   Xnew((nXnow+1):(nXnow+nXnew),:) = Xtmp;

   nHnow = nHnow + nQnew;
   nHtot = nHtot + nQnew;
   nXnow = nXnow + nXnew;
   nXtot = nXtot + nXnew;
   Qfab0 = Qfab1;

end

Hexes=[Hexes;Hnew];
X = [X;Xnew];
Hbc = [Hbc;Hbcnew];

fprintf(' done! nH=%d (%2.4e sec)\n',size(Hexes,1),toc(t0));
