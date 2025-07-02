 %Polygon file
 %2nd column V: number of vertices for each polygon
 %3rd column x: x coordinates of the vertices; 4th column y: y coordinates
 %Shape file
 %1st column X: X coordinates; 2nd column Y: Y coordinates

ntotal=length(V); %Total number of coordinates
NN0=sum(1.0./V);
NN=ceil(NN0-0.01); % Number of polygons
S=zeros(NN,1); C=zeros(NN,1); Roundness=zeros(NN,1); vertices=zeros(NN,1); % polygon area, perimeter, roundness and vertices number
Results1=zeros(NN,8); Results2=0; xc_normalized=zeros(NN,1); yc_normalized=zeros(NN,1); xcenter=zeros(NN,1);ycenter=zeros(NN,1);concave=zeros(NN,1);


% Calculate vertices number of each polygon
iv=1.0./V;
for ii=1:ntotal
       Nvertices(ii)=ceil(sum(iv(1:ii))-0.01);
end

for j=1:NN
   countfirst=0;
   for i=1:ntotal
      if (Nvertices(i)==j)
          countfirst=countfirst+1;
          if (countfirst<2)    %prevent reading typos
            vertices(j)=V(i);
          end
      end
   end
end


%Calculate geometrical entropy
En=0;
for I=min(vertices):max(vertices)
    if (sum(vertices==I)>0)
    CC(I)=sum(vertices==I)/NN;
    En=En-log(sum(vertices==I)/NN)*sum(vertices==I)/NN;
    end 
end


% Calculate area, perimeter, geometrical center and roundness of each polygon
for k=1:NN
    if k==1
        sumnumber=0;
    else
        sumnumber=sum(vertices(1:k-1));
    end 
    for p=sumnumber+1:sumnumber+vertices(k)-1
       S(k)=S(k)+0.5*(x(p)*y(p+1)-x(p+1)*y(p));
       C(k)=C(k)+sqrt((x(p)-x(p+1))^2+(y(p)-y(p+1))^2);
    end
    S(k)=S(k)+0.5*(x(sumnumber+vertices(k))*y(sumnumber+1)-x(sumnumber+1)*y(sumnumber+vertices(k)));  %area
    S(k)=abs(S(k));
    C(k)=C(k)+sqrt((x(sumnumber+vertices(k))-x(sumnumber+1))^2+(y(sumnumber+vertices(k))-y(sumnumber+1))^2);  %perimeter
    Roundness(k)=4*pi*S(k)/(C(k)^2);
    xcenter(k)=sum(x(sumnumber+1:sumnumber+vertices(k)))/vertices(k);  %geometrical center of each polygon
    ycenter(k)=sum(y(sumnumber+1:sumnumber+vertices(k)))/vertices(k); 

    % Judge whether it is a concave polygon
    if vertices(k)==3   % All triangles are convex polygon!
        concave(k)=0;
    else
        sumsig=0;
        for l=sumnumber+1:sumnumber+vertices(k)-2         
            sumsig=sumsig+sign((x(l+1)-x(l))*(y(l+2)-y(l+1))-(x(l+2)-x(l+1))*(y(l+1)-y(l)));
        end
        sumsig=sumsig+sign((x(sumnumber+vertices(k))-x(sumnumber+vertices(k)-1))*(y(sumnumber+1)-y(sumnumber+vertices(k)))-(x(sumnumber+1)-x(sumnumber+vertices(k)))*(y(sumnumber+vertices(k))-y(sumnumber+vertices(k)-1)));
        sumsig=sumsig+sign((x(sumnumber+1)-x(sumnumber+vertices(k)))*(y(sumnumber+2)-y(sumnumber+1))-(x(sumnumber+2)-x(sumnumber+1))*(y(sumnumber+1)-y(sumnumber+vertices(k))));  
        if(abs(sumsig)==vertices(k))
            concave(k)=0;
        else
            concave(k)=1;
        end
    end

end

%Calculate the area of a single valve
Stotal=0;
for q=1:length(X)-1
    Stotal=Stotal+0.5*(X(q)*Y(q+1)-X(q+1)*Y(q));
end
Stotal=Stotal+0.5*(X(length(X))*Y(1)-X(1)*Y(length(X)));
Stotal=abs(Stotal);  %should be bigger than sum(S) because polygons with ambiguous vertices aren't measured

%Normalization
S_normalized=S./Stotal;
S_normalized2=S./max(S);
xc_normalized=(xcenter-0.5*(max(X)+min(X)))/(max(X)-min(X));
yc_normalized=(ycenter-0.5*(max(Y)+min(Y)))/(max(Y)-min(Y));

Results1(:,1)=S_normalized; Results1(:,2)=S; Results1(:,3)=vertices; Results1(:,4)=concave; 
Results1(:,5)=S_normalized2; Results1(:,6)=Roundness; Results1(:,7)=xc_normalized;  Results1(:,8)=yc_normalized;
Results2(1)=mean(S_normalized); Results2(2)=std(S_normalized); Results2(3)=mean(Roundness); Results2(4)=std(Roundness);
Results2(5)=mean(vertices); Results2(6)=std(vertices); Results2(7)=sum(concave)/NN; Results2(8)=Stotal; Results2(9)=mean(S);  Results2(10)=std(S);
Results2(11)=En;

plot(xc_normalized,S_normalized,'o');
