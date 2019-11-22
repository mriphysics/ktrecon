function transformTime

load('/home/lcg13/Work/DataDefinitiveImplementationDebug06/transform.mat','FT','FTH','xp','et');

NN=128;

gpuF=2;
xp=xp(1:NN,1:NN,1:NN);
et{1}=et{1}(1:NN,1:NN,1:NN,:,:,1);
for n=2:3
    for m=1:3                
        N=size(et{n}{m});
        et{n}{m}=et{n}{m}(1:min(NN,N(1)),1:min(NN,N(2)),1:min(NN,N(3)),:,:,1);        
    end
end
for m=1:3
    FT{m}=FT{m}(1:NN,1:NN);
    FTH{m}=FTH{m}(1:NN,1:NN);
end
xp=complex(xp);
dev=gpuDevice;
for n=1:5
wait(dev);tic
y=real(sincRigidTransform(xp,et,1,FT,FTH));
wait(dev);toc
end