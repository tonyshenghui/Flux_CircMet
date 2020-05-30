function myvar=myvariance(X,M,Y)
 Y2=M*X;
 myvar=sum((Y2-Y).^2);