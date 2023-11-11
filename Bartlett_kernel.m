function y=Bartlett_kernel(x)

   y=(1-abs(x)).*(abs(x)<=1);
   
end