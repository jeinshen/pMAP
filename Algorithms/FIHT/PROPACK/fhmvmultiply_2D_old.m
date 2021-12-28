function y=fhmvmultiply_2D(h,x,q1,q2)

[n1,n2]=size(h);
lh=n1*n2;
h=reshape(h,lh,1);

p1=n1+1-q1;
p2=n2+1-q2;

lx=length(x);

xrev=x(lx:-1:1);

ind=zeros(lx,1);
for i=1:q2
    ind((i-1)*q1+1:i*q1)=(i-1)*n1+1:(i-1)*n1+q1;
end

xx=zeros(lh,1);
xx(ind)=xrev;

yy=ifft(fft(h).*fft(xx));

ind=zeros(p1*p2,1);
for i=1:p2
    ind((i-1)*p1+1:i*p1)=(q2+i-2)*n1+q1:(q2+i-1)*n1;
end

y=yy(ind);