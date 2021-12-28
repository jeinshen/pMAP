function y=fhmvmultiply_3D(h,x,q1,q2,q3,ind1,ind2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fast three-level Hankel matrix vector mulitplication.
%
%Inputs
% h: 3D array that generates the three-level Hankel matrix H.
% x: column vector to be left multiplied by H.
% q1: column dimension of level 1.
% q2: column dimension of level 2.
% q3: column dimension of level 3.
% ind1 (optional): index for x after padding and reversing order. 
% ind2 (optional): index for y, to be extracted after ifft.
%
%Outputs
% y: H*x.  
%
%Example
% h=rand(128,64,32);x=rand(65*33*17,1);y=fhmvmultiply_3D(h,x,65,33,17), 
% where y is the multiplication of the three-level Hankel matrix of size
% (64*32*16)*(65*33*17) formed by h with the vector x.
%
%Reference: Lu L, Xu W, Qiao S. A fast SVD for multilevel block Hankel
%matrices with minimal memory storage[J]. Numerical Algorithms, 2015, 
%69(4): 875-891.
%
%Last modified: 15-March-2017
%Please email tianming-wang@uiowa.edu for bug report.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n1,n2,n3]=size(h);

lh=n1*n2*n3;
h=reshape(h,lh,1);

lx=length(x);
xrev=x(lx:-1:1);

p1=n1+1-q1;
p2=n2+1-q2;
p3=n3+1-q3;

nargs=nargin;
if nargs==5
    ind1=zeros(q1*q2*q3,1);
    for i=1:q3
        ix=(i-1)*q1*q2;
        ixx=(i-1)*n1*n2;
        for j=1:q2
            ind1(ix+(j-1)*q1+1:ix+j*q1)=ixx+(j-1)*n1+1:ixx+(j-1)*n1+q1;
        end
    end
    ind2=zeros(p1*p2*p3,1);
    for i=1:p3
        iy=(i-1)*p1*p2;
        iyy=(q3+i-2)*n1*n2;
        for j=1:p2
            ind2(iy+(j-1)*p1+1:iy+j*p1)=iyy+(q2+j-2)*n1+q1:iyy+(q2+j-1)*n1;
        end
    end
end

xx=zeros(lh,1);
xx(ind1)=xrev;

L=2^nextpow2(lh);

fft_h=fft(h,L);
fft_xx=fft(xx,L);

yy=ifft(fft_h.*fft_xx);
y=yy(ind2);