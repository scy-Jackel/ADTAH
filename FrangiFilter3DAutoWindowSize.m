function[Iout,whatScale,Voutx,Vouty,Voutz]= FrangiFilter3DAutoWindowSize(I,options,window_size) 
 
defaultoptions = struct('FrangiScaleRange', [1 10], 'FrangiScaleRatio', 2, 'FrangiAlpha', 0.3,'FrangiBeta', 0.5, 'FrangiC', 1000, 'verbose',true,'BlackWhite',true); 
 
% Process inputs 
if(~exist('options','var')) 
    options=defaultoptions;  
else  
    tags = fieldnames(defaultoptions); 
    for i=1:length(tags) 
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end  
    end
    if(length(tags)~=length(fieldnames(options)))
        warning('FrangiFilter3D:unknownoption','unknown options found'); 
    end 
end 
 
% Use single or double for calculations 
if(~isa(I,'double')), I=single(I); end 
 
sigmas = options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2);  
sigmas = sort(sigmas, 'ascend');  
 
% Frangi filter for all sigmas  
for i = 1:length(sigmas) 
    % Show p rogress 
    if(options.verbose) 
        disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]); 
    end     
     
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3DDifferenceAutoSize(I,sigmas(i),window_size);   
 
    if(sigmas(i)>0) 
        % Correct for scaling 
        c=(sigmas(i)^2); 
        Dxx = c*Dxx; Dxy = c*Dxy; 
        Dxz = c*Dxz; Dyy = c*Dyy; 
        Dyz = c*Dyz; Dzz = c*Dzz; 
    end 
     
    % Calculate eigen values  
    if(nargout>2) 
        [Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz); 
    else 
        [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz); 
    end 
     
    % Free memory 
    clear Dxx Dyy Dzz Dxy Dxz Dyz; 
 
    % Calculate absolute values of eigen values   
    LambdaAbs1=abs(Lambda1); 
    LambdaAbs2=abs(Lambda2); 
    LambdaAbs3=abs(Lambda3); 
 
    % The Vesselness Features 
    Ra=LambdaAbs2./LambdaAbs3; 
    Rb=LambdaAbs1./sqrt(LambdaAbs2.*LambdaAbs3); 
 
    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D 
    S = sqrt(LambdaAbs1.^2+LambdaAbs2.^2+LambdaAbs3.^2); 
    A = 2*options.FrangiAlpha^2; B = 2*options.FrangiBeta^2;  C = 2*options.FrangiC^2; 
    C = 30;
  
    % Free memory clear LambdaAbs1 LambdaAbs2 LambdaAbs3 
 
    %Compute Vesselness function 
    expRa = (1-exp(-(Ra.^2./A))); 
    expRb = exp(-(Rb.^2./B)); 
    expS = (1-exp(-S.^2./(C))); 
    % keyboard 
    % Free memory 
    clear S A B C Ra Rb 
 
    %Compute Vesselness function 
    Voxel_data = expRa.* expRb.* expS; 
     
    % Free memory 
    clear expRa expRb expRc; 
     
    if(options.BlackWhite) 
        Voxel_data(Lambda2 < 0)=0; Voxel_data(Lambda3 < 0)=0; 
    else 
        Voxel_data(Lambda2 > 0)=0; Voxel_data(Lambda3 > 0)=0; 
    end 
         
    % Remove N aN values 
    Voxel_data(~isfinite(Voxel_data))=0;    
     
    if(i==1) 
        Iout=Voxel_data; 
        if(nargout>1) 
            whatScale = ones(size(I),class(Iout)); 
        end 
        if(nargout>2) 
            Voutx=Vx; Vouty=Vy; Voutz=Vz; 
        end 
    else 
        if(nargout>1) 
            whatScale(Voxel_data>Iout)=i; 
        end 
        if(nargout>2) 
            Voutx(Voxel_data>Iout)=Vx(Voxel_data>Iout); 
            Vouty(Voxel_data>Iout)=Vy(Voxel_data>Iout); 
            Voutz(Voxel_data>Iout)=Vz(Voxel_data>Iout); 
        end 
        % Keep maximum filter response 
        Iout=max(Iout,Voxel_data); 
    end 
end 
