function [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz)

[r,c,d] = size(Dxx);
for x = 1:r
    for y = 1:c
        for z = 1:d
            hessian(1,1) = Dxx(x,y,z);
            hessian(1,2) = Dxy(x,y,z);
            hessian(1,3) = Dxz(x,y,z);
            hessian(2,1) = Dxy(x,y,z);
            hessian(2,2) = Dyy(x,y,z);
            hessian(2,3) = Dyz(x,y,z);
            hessian(3,1) = Dxz(x,y,z);
            hessian(3,2) = Dyz(x,y,z);
            hessian(3,3) = Dzz(x,y,z);
            w = diag(hessian);
            [~,idx] = sort(abs(w));
            result = w(idx);
            Lambda1(x,y,z) = result(1);
            Lambda2(x,y,z) = result(2);
            Lambda3(x,y,z) = result(3);
        end
    end
end