
filename = 'lenac';
%for reading the RGB images 
image=rgb2gray(imread([filename '.png']));
%for reading bw images
% img=imread([filename '.png']);
image = double(image)./255;
[row1, col1] = size(image);
%%%%%%%%%%non-parallel partitioning%%%%%%%%%%%%%%
img1=image(1:(row1)/2,1:(col1)/2);
img2=image(1:(row1)/2,(col1)/2+1:col1);
img3=image((row1)/2+1:row1,1:(col1)/2);
img4=image((row1)/2+1:row1,(col1)/2+1:col1);
tic;
result1=edgedetect_partitioning(img1);
result2=edgedetect_partitioning(img2);
result3=edgedetect_partitioning(img3);
result4=edgedetect_partitioning(img4);
toc;
result(1:(row1)/2,1:(col1)/2)= result1;
result(1:(row1)/2,(col1)/2+1:col1)=result2;
result((row1)/2+1:row1,1:(col1)/2)=result3;
result((row1)/2+1:row1,(col1)/2+1:col1)=result4;
imwrite(result,gray(256),'partitioning_result.jpg','jpg');

