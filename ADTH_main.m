clc
clear

%r = 5;
%ROW = 2048;
%COLUMN = 2048;
%DEPTH = 140;

fname= 'C:\Users\Jackel\Desktop\091201c1.tif.v3dpbd.tif';
 
info = imfinfo(fname);
num_images = numel(info);
for i=1:num_images
     a= imread(fname,i);
     I(:,:,i)=a;
end
Volume = I;
I1 = I;
[a,b,c] = size(Volume);
mean_volumn_bound = 15;

%  Ìî³ä°ûÌå 
mean_volumn_bound_cell =50; 

%  »ù´¡¾àÀë±ä»» 
I(Volume>mean_volumn_bound) = 0; 
I(Volume<mean_volumn_bound) = 1; 
D = bwdist(I);    %Çó¾àÀë±ä»» 
D = (D-min(D(:)))./(max(D(:))-min(D(:))); 
D = D*255+1;     
%  ¸ßãÐÖµ¾àÀë±ä»» 
I1(Volume>mean_volumn_bound_cell) = 0; 
I1(Volume<mean_volumn_bound_cell) = 1; 
D_high = bwdist(I1);    %Çó¾àÀë±ä»» 
D_high = (D_high-min(D_high(:)))./(max(D_high(:))-min(D_high(:))); 
D_high = D_high*255+1;     

%gai
D_Threshold = 5;
Volume(D_high>D_Threshold) = 255;
Volume(D_high<D_Threshold) = 0;

img_name = '.\fill\';
num_images = size(Volume);
for i=1:num_images(3)
     J = Volume(:,:,i);
     imwrite(uint8(J),[img_name,num2str(i),'.tif'],'WriteMode','append');
end

%gai

options.BlackWhite=false; %¼ì²â°×É« 

options.FrangiScaleRange=[1 1]; 
%-------------------------------------------------------------------- 
 
Ifiltered = FrangiFilter3DAutoWindowSize(Volume,options,window_size);

Ifiltered255 = Ifiltered*255;

%gai
img_name = '.\only\';
num_images = size(Ifiltered255);
for i=1:num_images(3)
     J = Ifiltered255(:,:,i);
     imwrite(uint8(J),[img_name,num2str(i),'.tif'],'WriteMode','append');
end
%gai

D_Threshold = 5;
Ifiltered255(D_high>D_Threshold) = 255;
img_name = '.\out5\';
num_images = size(Ifiltered255);
for i=1:num_images(3)
     J = Ifiltered255(:,:,i);
     imwrite(uint8(J),[img_name,num2str(i),'.tif'],'WriteMode','append');
end