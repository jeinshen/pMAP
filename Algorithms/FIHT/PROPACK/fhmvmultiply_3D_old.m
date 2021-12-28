function y=fhmvmultiply_3D(h,x,ind1,ind2)

[n1,n2,n3]=size(h);
lh=n1*n2*n3;
h=reshape(h,lh,1);

%p1=n1+1-q1;
%p2=n2+1-q2;
%p3=n3+1-q3;

lx=length(x);

xrev=x(lx:-1:1);

% IND=zeros(lx,1);
% for i=1:q3
%     ix=(i-1)*q1*q2;
%     ixx=(i-1)*n1*n2;
%     for j=1:q2
%         IND(ix+(j-1)*q1+1:ix+j*q1)=ixx+(j-1)*n1+1:ixx+(j-1)*n1+q1;
%     end
% end

xx=zeros(lh,1);
xx(ind1)=xrev;

L=2^nextpow2(lh);
fft_h=fft(h,L);
fft_xx=fft(xx,L);

yy=ifft(fft_h.*fft_xx);

%yy=ifft(fft(h).*fft(xx));

% IND=zeros(p1*p2*p3,1);
% for i=1:p3
%     iy=(i-1)*p1*p2;
%     iyy=(q3+i-2)*n1*n2;
%     for j=1:p2
%         IND(iy+(j-1)*p1+1:iy+j*p1)=iyy+(q2+j-2)*n1+q1:iyy+(q2+j-1)*n1;
%     end
% end

y=yy(ind2);