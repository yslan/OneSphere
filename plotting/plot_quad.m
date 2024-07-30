function plot_quad(ifig,X,Quad)
if(ifig==0);return;end; if(ifig<0); ifig=abs(ifig); hold off; end; figure(ifig);

   nX=size(X,1);
   id=min((Quad>0) & (Quad<=nX),[],2);
   Quad = Quad(id,:);

   nQ=size(Quad,1);

   Xq = reshape(X(Quad(:,[1,2,3,4,1]),1),nQ,5);
   Yq = reshape(X(Quad(:,[1,2,3,4,1]),2),nQ,5);
   Zq = reshape(X(Quad(:,[1,2,3,4,1]),3),nQ,5);

   fill3(Xq',Yq',Zq','w'); hold on
   plot3(Xq',Yq',Zq','k-'); hold on
   
%   for i=1:size(Quad,1)
%      Xq = X(Quad(i,[1,2,3,4,1]),:);
%      plot3(Xq(:,1),Xq(:,2),Xq(:,3),'k-'); hold on
%   end

   axis equal
