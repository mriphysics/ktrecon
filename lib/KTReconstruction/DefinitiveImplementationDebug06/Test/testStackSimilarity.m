I=imread('peppers.png');
x=I(:,:,1);
x=im2single(x);

N=size(x);N(end+1:3)=1;
x=resampling(x,2*N(1:2),2);
N=size(x);N(end+1:3)=1;
T=zeros(1,1,1,1,6);
T(4)=pi/2+0.0001;

[rGrid,kGrid,rkGrid,~,cGrid]=generateTransformGrids(N,[],[],round(3*N/8));
et=precomputeFactorsSincRigidTransform(kGrid,rkGrid,T,[],[],[],[],cGrid);
y=sincRigidTransform(x,et);
eth=precomputeFactorsSincRigidTransform(kGrid,rkGrid,T,0,[],[],[],cGrid);
z=sincRigidTransform(y,eth,0);

%figure
%imshow(abs(x),[0 1])
figure
imshow(abs(y),[0 1])
%figure
%imshow(abs(z),[0 1])