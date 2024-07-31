function Hcurve = fix_cyl_curve_via_BC(X,Hexes,Hbc,Hcurve,bcid_list,Rcyl);

fprintf('Curve:   fix cyl curves Rcyl=%g',Rcyl); t0 = tic;
nHex = size(Hexes,1);
nHcurve = size(Hcurve,3);
nbc = length(bcid_list);

ncurve0 = length(find(Hcurve(1,:,:)));
if (nHcurve<nHex);
   Hcurve = cat(3,Hcurve,zeros(6,12,nHex-nHcurve));
end
%return % TODO

iftoiedge=[1 10 2 9;2 11 6 10;3 12 7 11;4 9 8 12;1 2 3 4;5 6 7 8];
iedgetoiv=[1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5;1 5;2 6;3 7;4 8];
iftoiv=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;4 3 2 1;5 6 7 8];

fprintf(', bcid:');
for i=1:nbc
   bc = bcid_list(i); fprintf(' %3d',bc);
   [ihex,ifac]=find(Hbc==bc);

   for j=1:4 % edge
      iedge = iftoiedge(ifac,j);
      v1 = iedgetoiv(iedge,1);
      v2 = iedgetoiv(iedge,2);
      id1=sub2ind([nHex,8],ihex,v1);ix1=Hexes(id1);
      id2=sub2ind([nHex,8],ihex,v2);ix2=Hexes(id2);
      xmid = (X(ix1,:)+X(ix2,:))/2;

      % proj
      rmid = sqrt(xmid(:,1).^2+xmid(:,2).^2);
      xmid(:,1) = xmid(:,1)./rmid*Rcyl; 
      xmid(:,2) = xmid(:,2)./rmid*Rcyl; 

      itmp=sub2ind([6,12,nHex],1+0*ihex,iedge,ihex); Hcurve(itmp)=1;
      itmp=sub2ind([6,12,nHex],2+0*ihex,iedge,ihex); Hcurve(itmp)=xmid(:,1);
      itmp=sub2ind([6,12,nHex],3+0*ihex,iedge,ihex); Hcurve(itmp)=xmid(:,2);
      itmp=sub2ind([6,12,nHex],4+0*ihex,iedge,ihex); Hcurve(itmp)=xmid(:,3);
   end
end
ncurve1 = length(find(Hcurve(1,:,:)));

fprintf(' done! ncurv: bfr=%d aft=%d (%2.4e sec)\n',ncurve0,ncurve1,toc(t0));
