function Anggen = bgrape_make_ang_mom(s)
d = 2*s+1;
Anggen.jx=zeros(d);
for m=1:d
   for n=1:d
      if(m+1==n)
          Anggen.jx(d-m+1,d-n+1)=(1/2)*sqrt((d-m)*m);
           Anggen.jx(d-n+1,d-m+1)=(1/2)*sqrt((d-m)*m);
      end;
      
   end;
   
end;

Anggen.jy=zeros(d);
for m=1:d
   for n=1:d
      if(m+1==n)
          Anggen.jy(d-m+1,d-n+1)=-i*(1/2)*sqrt((d-m)*m);
           Anggen.jy(d-n+1,d-m+1)=i*(1/2)*sqrt((d-m)*m);
      end;
      
   end;
   
end;

Anggen.jz=diag([-s:s]);
clear m n d 
