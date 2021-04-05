%% general solution
clear all
% Loading of the data
load('proj_fit_14','-mat');

% The extraction of the identification data
x1=id.X{1};
x2=id.X{2};
y=id.Y;

% The extraction of the the validation data
x1val=val.X{1};
x2val=val.X{2};
yval=val.Y;

% Start of the variable degree function approximator
% count is used to iterrate through mse vector
count=1;

for index=1:30
    
Xflat(1:1681,1)=ones(1681,1);% making 1st coloumn of the regressor equal to 1

p=1;%used to generate the powers of x1 and x2 (up to index-1, this is not used for the powers of x1*x2)
k=1;%used to iterrate through x1
k1=1;%used to iterrate through x2

for j=0:id.dims-1
for i=id.dims*j+1:id.dims*j+id.dims
    %the two for loops are used to iterrate throught the vectors x1 and x2 in designd way, 41 by
    %41 for x1, and 1 by 1 for x2 for every iterration of x1
 for jcol=2:2*index+1 % this loop generates the powers of x1 and x2 (from 1 to index-1) for the regressor and places them in the regressor vector Xflat
  if(mod(jcol,2)==0) 
  Xflat(i,jcol)=x1(k).^p;
  p=p+1;
  else
  Xflat(i,jcol)=x2(k1).^(p-1);
  end
 end
  p=1;
  jcol=2*index+2; % we continue to generate the rest of the terms, as combinatos of x1*x2 at different powers (p1+p2<=index)
  for p1=1:index-1% p1,p2 are used for the powers of x1*x2 (p1+p2<=index)
      for p2=1:index-1
         if(p1+p2<=index)
         Xflat(i,jcol)=x1(k).^p1.*x2(k1).^p2;
         jcol=jcol+1;
         end
      end
  end
  k1=k1+1;
end
k=k+1;
k1=1;
end

% we repeat the above process for the validation data set

Xflatval(1:5041,1)=ones(5041,1);

p=1;
k=1;
k1=1;

for j=0:val.dims-1
for i=val.dims*j+1:val.dims*j+val.dims
  for jcol=2:2*index+1
  if(mod(jcol,2)==0) 
  Xflatval(i,jcol)=x1val(k).^p;
  p=p+1;
  else
  Xflatval(i,jcol)=x2val(k1).^(p-1);
  end
  end
  p=1;
  jcol=2*index+2;
  for p1=1:index-1
      for p2=1:index-1
        if(p1+p2<=index)
         Xflatval(i,jcol)=x1val(k).^p1.*x2val(k1).^p2;
         jcol=jcol+1;
        end    
      end
  end
  k1=k1+1;
end
k=k+1;
k1=1;
end

poz=1;
for i=1:id.dims% transformation of the output indentification data to a coloumn vector
    for j=1:id.dims
    Yflat(poz,1)=y(i,j);
    poz=poz+1;
    end
end

poz=1;
for i=1:length(yval)% transformation of the output validation data to a coloumn vector
    for j=1:length(yval)
        yvalT(poz,1)=yval(i,j);
        poz=poz+1;
    end
end

theta=Xflat\Yflat;% computation of the paramaters of the regressor from the identification data
Yhatval=Xflatval*theta;% computation of the new output from the validation data and the regressor parameters obtained
Yhatid=Xflat*theta;% computation of the new ouput from the identification data and the regressor parameters obtained

MSEid=0;% calculation of the MSE from the approximator and the validation data
for i=1:length(Yflat)
   MSEid=MSEid+1/length(Yflat)*((Yhatid(i,1)-Yflat(i,1)).^2);
end
mseid(count)=MSEid;% addition of every MSE to a mse vector for the purpose of determining the best order fit

MSE=0;% calculation of the MSE from the approximator and the validation data
for i=1:length(yvalT)
   MSE=MSE+1/length(yvalT)*((Yhatval(i,1)-yvalT(i,1)).^2);
end

mseval(count)=MSE;% addition of every MSE to a mse vector for the purpose of determining the best order fit
count=count+1

end

[ minim i]=min(mseval)% computation of the lowest MSE value and its coresponding order
m=4:1:20;

figure,plot(m,mseid(4:20))% plot for the comparison of the MSE on the identification data versus the degree
xlabel('The degree of the polynomial approximator');ylabel('MSE from identification data');title('MSEid vs ORDER');

figure,plot(m,mseval(4:20))% plot for the comparison of the MSE on the validation data versus the degree
xlabel('The degree of the polynomial approximator');ylabel('MSE from validation data');title('MSEval vs ORDER');
%% After we select the degree for the best approximator polynomial function we repeat the same process as before to compute the model
clear all
load('proj_fit_14','-mat');

x1=id.X{1};
x2=id.X{2};
y=id.Y;
yval=val.Y;
x1val=val.X{1};
x2val=val.X{2};
index=4;

Xflat(1:1681,1)=ones(1681,1);

p=1;
k=1;
k1=1;

for j=0:id.dims-1
for i=id.dims*j+1:id.dims*j+id.dims
  for jcol=2:2*index+1
  if(mod(jcol,2)==0) 
  Xflat(i,jcol)=x1(k).^p;
  p=p+1;
  else
  Xflat(i,jcol)=x2(k1).^(p-1);
  end
  end
  p=1;
  jcol=2*index+2;
  for p1=1:index-1
      for p2=1:index-1
         if(p1+p2<=index)
         Xflat(i,jcol)=x1(k).^p1.*x2(k1).^p2;
         jcol=jcol+1;
         end
      end
  end
  k1=k1+1;
end
k=k+1;
k1=1;
end

Xflatval(1:5041,1)=ones(5041,1);

p=1;
k=1;
k1=1;
for j=0:val.dims-1
for i=val.dims*j+1:val.dims*j+val.dims
  for jcol=2:2*index+1
  if(mod(jcol,2)==0) 
  Xflatval(i,jcol)=x1val(k).^p;
  p=p+1;
  else
  Xflatval(i,jcol)=x2val(k1).^(p-1);
  end
  end
  p=1;
  jcol=2*index+2;
  for p1=1:index-1
      for p2=1:index-1
        if(p1+p2<=index)
         Xflatval(i,jcol)=x1val(k).^p1.*x2val(k1).^p2;
         jcol=jcol+1;
        end    
      end
  end
  k1=k1+1;
end
k=k+1;
k1=1;
end

k=1;

for i=1:id.dims
    for j=1:id.dims
    Yflat(k,1)=y(i,j);
    k=k+1;
    end
end

theta=Xflat\Yflat;
Yhatval=Xflatval*theta;

Yhatval=reshape(Yhatval,71,71);

figure(3);
mesh(x1val,x2val,Yhatval), xlabel("x1val"); ylabel("x2val"); zlabel("Yhatval")% plot of the obtained model for the approximator
title('Obtained model for the polynomial approximator');

[x1plot,x2plot]=meshgrid(x1,x2);
figure(1);
mesh(x1plot,x2plot,y);hold, xlabel("x1"); ylabel("x2"); zlabel("y");% plot of the identification data set
title('Identification data set');

[x1plotval,x2plotval]=meshgrid(x1val,x2val);
figure(2);
mesh(x1plotval,x2plotval,yval);hold, xlabel("x1val"); ylabel("x2val"); zlabel("yval")% plot of the validation data set
title('Validation data set');

figure(4);
surf(x1val,x2val,Yhatval);% used to compare the obtained model with the validation data
hold on

surf(x1val,x2val,yval','FaceAlpha',0.5), xlabel("x1val"); ylabel("x2val"); zlabel("yval")% used to compare the obtained model with the validation data
title('Comparison between the obtained model and the validation data set');