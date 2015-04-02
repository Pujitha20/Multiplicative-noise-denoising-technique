clear all;clc;close all;
im=imread('C:\Users\AK PUJITHA\Desktop\lena.jpg');

[q,w]=size(size(im));
if w==3
    im=rgb2gray(im);
end
[m n]=size(im);
%Level of decoposition is 3
lvl=3;
%Type of wavelet
wavelet='db4';
%Set window size for computing variance,w1
w1=5;
%Set window size for detecting features,w2
w2=7;
%Set beta
beta=0.15;
%Set eta
eta=5;
[C,S] = wavedec2(im,lvl,wavelet);
figure()
sq1=plotwavelet2(uint16(C),S,lvl,wavelet,256,'square');
gamma2=[0.953 0.953 0.953];
gamma1=[0.1257 0.255 0.417];
for k = 1:lvl
    cA{k} = appcoef2(C,S,wavelet,k); % approx
    [cH{k} cV{k} cD{k}] = detcoef2('a',C,S,k); % details
 
    oA{k}=cA{k};
    oH{k}=cH{k};
    oV{k}=cV{k};
    oD{k}=cD{k};
    H=oH;
    V=oV;
    D1=oD;
    D2=oD;
    [m1,n1]=size(H{k});
    %Modifing H
    Hv=H{k}.^2;
    Hv=conv2(Hv,[zeros((w1-1)/2,w1);ones(1,w1);zeros((w1-1)/2,w1)]);

    for i=1:n1
        for j=1:m1
            if j==1
                R=Hv(j+1,i)/Hv(j,i);
            elseif j==m1
                R=Hv(j-1,i)/Hv(j,i);
            else
                R=max(Hv(j+1,i),Hv(j-1,i))/Hv(j,i);
            end
            if R<=gamma1(k)
                alpha=1;
            else
                if j==1
                    med1=NaN;
                    med2=median(Hv(j+1:j+(w2-1)/2,i));
                elseif j-(w2-1)/2<=0
                    med1=median(Hv(1:j-1,i));
                    med2=median(Hv(j+1:j+(w2-1)/2,i));
                elseif j==m1
                    med1=median(Hv(j-(w2-1)/2:j-1,i));
                    med2=NaN;
                elseif j+(w2-1)/2>m1
                    med1=median(Hv(j-(w2-1)/2:j-1,i));
                    med2=median(Hv(j+1:m1,i));
                else
                    med1=median(Hv(j-(w2-1)/2:j-1,i));
                    med2=median(Hv(j+1:j+(w2-1)/2,i));
                end
                R=min(med1,med2)/Hv(j,i);
                if R>gamma2(k)
                    alpha=beta;
                elseif R<=gamma1(k)
                    alpha=1;
                else
                    alpha=beta+(1-beta)*(1-exp(eta*(R-gamma2(k))/(gamma1(k)-gamma2(k))))/(1-exp(eta));
                end
            end
        H{k}(j,i)=alpha*H{k}(j,i);
        end
    end
    %Modifing V
    Vv=V{k}.^2;
    Vv=conv2(Vv,[zeros((w1-1)/2,w1);ones(1,w1);zeros((w1-1)/2,w1)]');
   
    for i=1:m1
        for j=1:n1
            if j==1
                R=Vv(i,j+1)/Vv(i,j);
            elseif j==m1
                R=Vv(i,j-1)/Vv(i,j);
            else
                R=max(Vv(i,j+1),Vv(i,j-1))/Vv(i,j);
            end
            if R<=gamma1(k)
                alpha=1;
            else
                if j==1
                    med1=NaN;
                    med2=median(Vv(i,j+1:j+(w2-1)/2));
                elseif j-(w2-1)/2<=0
                    med1=median(Vv(i,1:j-1));
                    med2=median(Vv(i,j+1:j+(w2-1)/2));
                elseif j==n1
                    med1=median(Vv(i,j-(w2-1)/2:j-1));
                    med2=NaN;
                elseif j+(w2-1)/2>n1
                    med1=median(Vv(i,j-(w2-1)/2:j-1));
                    med2=median(Vv(i,j+1:n1));
                else
                    med1=median(Vv(i,j-(w2-1)/2:j-1));
                    med2=median(Vv(i,j+1:j+(w2-1)/2));
                end
                R=min(med1,med2)/Vv(j,i);
                if R>gamma2(k)
                    alpha=beta;
                elseif R<=gamma1(k)
                    alpha=1;
                else
                    alpha=beta+(1-beta)*(1-exp(eta*(R-gamma2(k))/(gamma1(k)-gamma2(k))))/(1-exp(eta));
                end
            end
            V{k}(j,i)=alpha*V{k}(j,i);
        end
    end
    %Modifing D1
    Dv1=D1{k}.^2;
    Dv1=conv2(Dv1,diag(ones(w1)));
 
    for i=1+(w2-1)/2:m1-(w2-1)/2
        for j=1+(w2-1)/2:n1-(w2-1)/2
            R=max(Dv1(i+1,j+1),Dv1(i-1,j-1))/Dv1(i,j);
            if R<=gamma1(k)
                alpha=1;
            else
                med1=median(diag(Dv1(i+(w2-1)/2:-1:i+1,j-(w2-1)/2:j-1)));
                med2=median(diag(Dv1(i-1:-1:i-(w2-1)/2,j+1:j+(w2-1)/2)));
                R=min(med1,med2)/Dv1(j,i);
                if R>gamma2(k)
                    alpha=beta;
                elseif R<=gamma1(k)
                    alpha=1;
                else
                    alpha=beta+(1-beta)*(1-exp(eta*(R-gamma2(k))/(gamma1(k)-gamma2(k))))/(1-exp(eta));
                end
            end
            D1{k}(i,j)=alpha*D1{k}(i,j);
        end
    end
    %Modifing D2
    D2{k}=D2{k}(m1:-1:1,:);
    Dv2=D2{k}.^2;
    Dv2=conv2(Dv2,diag(ones(w1)));
 
    for i=1+(w2-1)/2:m1-(w2-1)/2
        for j=1+(w2-1)/2:n1-(w2-1)/2
            R=max(Dv2(i+1,j+1),Dv2(i-1,j-1))/Dv2(i,j);
            if R<=gamma1(k)
                alpha=1;
            else
                med1=median(diag(Dv2(i+(w2-1)/2:-1:i+1,j-(w2-1)/2:j-1)));
                med2=median(diag(Dv2(i-1:-1:i-(w2-1)/2,j+1:j+(w2-1)/2)));
                R=min(med1,med2)/Dv2(j,i);
                if R>gamma2(k)
                    alpha=beta;
                elseif R<=gamma1(k)
                    alpha=1;
                else
                    alpha=beta+(1-beta)*(1-exp(eta*(R-gamma2(k))/(gamma1(k)-gamma2(k))))/(1-exp(eta));
                end
            end
            D2{k}(i,j)=alpha*D2{k}(i,j);
        end
    end
    D2{k}=D2{k}(m1:-1:1,:);
    for i=1:m1
        for j=1:n1
            if(oD{k}(i,j)~=D1{k}(i,j) && oD{k}(i,j)~=D2{k}(i,j))
                oD{k}(i,j)=D1{k}(i,j);
                oH{k}(i,j)=H{k}(i,j);
                oV{k}(i,j)=V{k}(i,j);
            end
        end
    end
end
newC=reshape(oA{3},1,size(oA{k},1)^2);
for k=lvl:-1:1
    newC=horzcat(newC,reshape(oH{k},1,size(oH{k},1)^2),reshape(oV{k},1,size(oV{k},1)^2),reshape(oD{k},1,size(oD{k},1)^2));
end
out=waverec2(newC,S,'db4');
% % %Reconstructing image
figure()
imshow(uint8(out))
% figure()
% imshow(im)