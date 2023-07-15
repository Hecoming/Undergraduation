clear;
cryscipt=fopen("1.txt");%open file
fileInfo=dir("1.txt");%get file size, to square decomposition
filesize=fileInfo.bytes;
[a,b]=splitIntoSquares(filesize);%square decomposition function,sometimes it error,but i have no time to get a better one
A2=fread(cryscipt,a^2,'int8');%devide the general data to 2 matrix
B2=fread(cryscipt,b^2,'int8');
fclose(cryscipt);%close file
A1=[a a];B1=[b b];%reshape the data
A=reshape(A2,A1);B=reshape(B2,B1);
[keyA,midA,ciphertextA,HVD_plainA,HVD_cipherA,ANPCR]=JCS(A);
[keyB,midB,ciphertextB,HVD_plainB,HVD_cipherB,BNPCR]=JCS(B);

function [key0,mid0,ciphertext0,HVD_plain,HVD_cipher,NPCRs]=JCS(image1)%this function provide the method and the indicators
    p1=imresize(image1,[64 64]);%to shorten the computing time
    p=sym(int16(p1));
    [key1,ciphertext]=jordan(p);
    [key,mid]=Demo_CS_OMP(key1);%1d compressed sensing,deal data 1by1 row
    key0=real(double(key));
    mid0=real(double(mid));
    ciphertext0=real(double(ciphertext));
    plaintext0=real(double(p1));
    histogram(ciphertext0);histogram(plaintext0);%generate histogram of plain and cipher
    plain_H=Correlation_of_adjacent_pixels(plaintext0,1,200);%Horizontal correlation,the third indicator is the number of parameters
    plain_V=Correlation_of_adjacent_pixels(plaintext0,2,200);%Vertical correlation
    plain_D=Correlation_of_adjacent_pixels(plaintext0,3,200);%Diagonal correlation
    cipher_H=Correlation_of_adjacent_pixels(ciphertext0,1,200);
    cipher_V=Correlation_of_adjacent_pixels(ciphertext0,2,200);
    cipher_D=Correlation_of_adjacent_pixels(ciphertext0,3,200);
    [i0,j0]=size(p1);
    i1=randi([1,i0],1,1);j1=randi([1,j0],1,1);
    p2=p1;p2(i1,j1)=p1(i1,j1)+randi([1,100],1,1);%use randi to change a random data in matrix
    pp=sym(double(p2));[key2,ciphertext2]=jordan(pp);
    [keykey,midmid]=Demo_CS_OMP(key2);
    key00=real(double(keykey));
    mid00=real(double(midmid));
    ciphertext00=real(double(ciphertext2));
    key_NPCR=NPCR(key0,key00);
    cipher_NPCR=NPCR(ciphertext0,ciphertext00);
    mid_NPCR=NPCR(mid0,mid00);%compute NPCR
    HVD_plain=[plain_H,plain_V,plain_D];HVD_cipher=[cipher_H,cipher_V,cipher_D];
    NPCRs=[key_NPCR,cipher_NPCR,mid_NPCR];
end

function [x, y] = splitIntoSquares(n)
    x = 0;
    y = int16(sqrt(n)); 
    while x <= y 
        if x^2 + y^2 == n
            return;
        elseif x^2 + y^2 < n 
            x = x + 1;
        else
            y = y - 1;
        end
    end
    x = -1;
    y = -1;
    error("file size can't make decomposition");
end

function [img_cs_1d,Phi]=Demo_CS_OMP(img)
img=double(int16(img));
[height,width]=size(img);
N=height-1;
Phi=randn(floor(0.7*height),width);  % 简化测量矩阵 
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(0.7*height),1]); % 规范每行
mat_dct_1d=zeros(height,height);  % 生成DCT基矩阵
for k=0:1:N 
    dct_1d=cos([0:1:N]'*k*pi/height);
    if k>0
        dct_1d=dct_1d-mean(dct_1d); 
    end
    mat_dct_1d(:,k+1)=dct_1d/norm(dct_1d);
end
img_cs_1d=Phi*img;  
sparse_rec_1d=zeros(height,width);  
Theta_1d=Phi*mat_dct_1d;%测量矩阵乘上基矩阵
 
for i=1:width
    column_rec=cs_omp(img_cs_1d(:,i),Theta_1d,height);
    sparse_rec_1d(:,i)=column_rec'; %  稀疏系数
end
end

function hat_x=cs_omp(y,T_Mat,m)

n=length(y);
s=floor(3*n/4); %  测量值维数
hat_x=zeros(1,m); %  待重构的谱域(变换域)向量                     

Aug_t=[];        %  增量矩阵(初始值为空矩阵)
r_n=y;  %  残差值 


for times=1:s %  迭代次数(稀疏度是测量的1/4)

    product=abs(T_Mat'*r_n);    
    [val,pos]=max(product);   %最大投影系数对应的位置
    Aug_t=[Aug_t,T_Mat(:,pos)];   %矩阵扩充
    T_Mat(:,pos)=zeros(n,1); %选中的列置零
    aug_x=(Aug_t'*Aug_t)^(-1)*Aug_t'*y;  % 最小二乘
    r_n=y-Aug_t*aug_x;   %残差
    pos_array(times)=pos;   %纪录最大投影系数的位置

end
end

function l=Correlation_of_adjacent_pixels(image,choose,n)
%抽取n对相邻像素
%choose 选择1水平，2垂直，3对角
%n 抽样对数
image=double(image);
[M,N]=size(image);%M行N列

x_coor(1,:)=randi([1  N],1,n);%x序列x坐标
x_coor(2,:)=randi([1  M],1,n);%x序列y坐标
y_coor=ones(2,n);%y序列坐标


if choose==1
%水平
for i=1:n
    if x_coor(1,i)==N
        y_coor(1,i)=1;
    end
    if x_coor(1,i)<N
        y_coor(1,i)=x_coor(1,i)+1;
    end
    y_coor(2,i)=x_coor(2,i);
end
end


if choose==2
%垂直
for i=1:n
    if x_coor(2,i)==M
        y_coor(2,i)=1;
    end
    if x_coor(2,i)<M
        y_coor(2,i)=x_coor(2,i)+1;
    end
    y_coor(1,i)=x_coor(1,i);
end
end


if choose==3
%对角
for i=1:n
    if x_coor(1,i)
        y_coor(1,i)=1;
    end
    if x_coor(1,i)<N
        y_coor(1,i)=x_coor(1,i)+1;
    end
    
    if x_coor(2,i)==M
        y_coor(2,i)=1;
    end
    if x_coor(2,i)<M
        y_coor(2,i)=x_coor(2,i)+1;
    end
end
end

x=ones(1,n);
y=ones(1,n);
%获取像素值
for i=1:n
    x(i)=image(x_coor(2,i),x_coor(1,i));
    y(i)=image(y_coor(2,i),y_coor(1,i));
end

l=coefficient_of_association(x,y);
end

function l=coefficient_of_association (x,y)
%求相关系数
mean_x=mean(x);
mean_y=mean(y);
n=length(x);
up=0;
sum_x=0;
sum_y=0;
for i=1:n
    up=up+(x(i)-mean_x)*(y(i)-mean_y);
    sum_x=sum_x+(x(i)-mean_x)^2;
    sum_y=sum_y+(y(i)-mean_y)^2;
end
down=sqrt(sum_x*sum_y);
l=up/down;
end

function NPC=NPCR(image1,image2)
    [iN,jN]=size(image1);
    m=0;
    for i=1:iN
        for j=1:jN
            if image1(i,j)==image2(i,j)
            m=m+1;
            end
        end
    end
    ij=iN*jN;
    NPC=(ij-m)/(ij);
end