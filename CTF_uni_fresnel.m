clc
clear all;
close all;
NA=0.42;
lambda=660*10^-9;
npad=0;
N=1200 + 2*npad;
no=N/2;
Mag=20;
dpix=(4.4*10^-6)/Mag;
NoiseLevel=10;
z=10*10^-6;

Kcut=1;
Kotf=NA/lambda;
k=2*pi/lambda;
x=(-N/2:N/2-1)*dpix;
y=(-N/2:N/2-1)*dpix;
[X,Y]=meshgrid(x,y);
du = 1/(N*dpix);
u =(-N/2:N/2-1)*du;
v=(-N/2:N/2-1)*du;
[U,V]=meshgrid(u,v);
K=sqrt(U.^2+V.^2);
K_sq=(U.^2+V.^2);

H=zeros(N,N);
H(K<=(0.5*Kotf))=1;
H=ifftshift(ifft2((abs(fft2((H))).^2)) );
H= H./(max(max(H)));

sigma=200000;
G=(1/sqrt(2*pi*sigma*sigma)).*(exp(-K_sq./(2*sigma*sigma)));
G=G./max(max(G));

figure
hold on
plot(u,abs(H(no+1,:)),'--','LineWidth',2,'MarkerSize',6)
plot(u,abs(G(no+1,:)),'--','LineWidth',2,'MarkerSize',6)
legend('H','G')
xline(Kotf)
xlabel("u");
grid on
box on

OTF=G;



Imagee=double(imread("gold.png"));
Imagee = padarray(Imagee,[npad npad],0,"both");

P= double(Imagee);
P=rescale(P,-0.2,0.3);
T=ones(N);
Object=sqrt(T).*exp(1j.*P);

figure
imagesc(x*10^3,y*10^3,(P));
colormap(gray)
axis off



% w = (2*pi/lambda).*( 1 - sqrt(1- (lambda^2.*(U.^2+V.^2))) )* z ;
% 
% prop_phs= 1i*(2*pi/lambda)*sqrt(1- (lambda^2.*(U.^2+V.^2)));
% prop_kernel  =exp(prop_phs * z);
% 
% prop_kernel3 = exp(1j* (2*pi/lambda)* z).*exp(-1j*w) ;
% 

% Ip = real(ifft2(ifftshift((fftshift(fft2(Object)).*prop_kernel3))));
% [Ip]=OTFatten(Ip,OTF);
% 
% Im = real(ifft2(ifftshift((fftshift(fft2(Object)).*conj(prop_kernel3)))));
% [Im]=OTFatten(Im,OTF);
% 
% figure
% imagesc(x*10^3,y*10^3,(Ip ));
% colormap(gray)
% axis off
% 
% 
% figure
% imagesc(x*10^3,y*10^3,(Im ));
% colormap(gray)
% axis off


prop_kernel  = exp(-1j*pi*lambda*z.*(U.^2+V.^2));

Ip = real(ifft2(ifftshift((fftshift(fft2(Object)).*prop_kernel))));
[Ip]=OTFatten(Ip,OTF);

Im = real(ifft2(ifftshift((fftshift(fft2(Object)).*conj(prop_kernel)))));
[Im]=OTFatten(Im,OTF);

figure
imagesc(x*10^3,y*10^3,(Ip ));
colormap(gray)
axis off

figure
imagesc(x*10^3,y*10^3,(Im ));
colormap(gray)
axis off


delI_FT = fftshift(fft2(Ip)) - fftshift(fft2(Im));
figure
imagesc(x*10^3,y*10^3,log(abs(delI_FT)));
colormap(gray)
axis off

D=(4.*sin(pi*lambda*z.*(U.^2+V.^2)));
D(isnan(D))=1;
D(isinf(D))=1;
figure
imagesc(D);
colormap(jet)
colorbar

%%%window 
W=zeros(N,N);
W=K>Kcut*Kotf;
W=double(W);

[Al,al]=Object_Para_Estimation(delI_FT,D,Kotf,U,V,W,K,OTF);
Objamp = Al.*(K.^-al);
Sigamp = Objamp.*D.*OTF;

figure
hold on
plot(u*10^-3,log(abs(delI_FT(no+4,:))),'k--')
plot(u*10^-3,log(abs(Sigamp(no+4,:))),'m-')
xlabel("u(cycles/mm)");
legend('Acutal signal power','Estimated. signal power')
title("Power spectrum estimation")
grid on
box on
%% 

[WF_P_uni,Noise_F_P_uni]=Weiner_filter_center(delI_FT,Al,al,D,K,W,OTF);
WP_uni=real(ifft2(ifftshift(WF_P_uni)));

figure
imagesc(x*10^3,y*10^3,(WP_uni));
xlabel("x(mm)");
ylabel("y(mm)");
title("Weiner filtered Phase using Gaussian process TIE")
colormap(gray)
colorbar


%%Uniform TIE
function [I]=FresProp(dpix,z,lambda,Hsize,Hcrop)
Xsize = Hsize*dpix; %Hsize is the number of pixel, dpix is the length of one pixel, Xsize is the total lenght of the image. 
du = 1/(Xsize);% the resolution of fourier frequency coordinates
umax = 1/(2*dpix); %define the k space 
u = -umax:du:umax-du;
[U,V]=meshgrid(u,u);
clear  u V  du;
%Evanescent cut-off 
uev = 1/lambda; %???????
%Nyquist cut-off for Fresnel Propagation Kernel
unp = uev*(Xsize/(2*abs(z)));
clear Xsize;
%Circular window
A = U.^2+(U').^2;
clear U;
if uev>=unp
    ucut = unp;
end
if unp>uev
    ucut = uev;
end
W= sqrt(A);
W = (W<=ucut); 
% disp(['Cutoff =',num2str(ucut),' Evansecent Cutoff =',num2str(uev),...
%' Nyquist Cutoff =', num2str(unp),'u max =',num2str(umax)])
clear ucut uev unp
%Fresnel kernel: paraxial approximation
H = exp((-1i*pi*lambda* z).*(A));
clear A;
%Truncate kernel
H = W.*H;
clear W;
%Hologram Spectrum
Htemp = fft2(Hcrop);
HH = fftshift(Htemp);
clear Htemp;
%Propagate field
RR = HH.*H;
clear H HH;
RR =ifftshift(RR);
Ri = ifft2(RR); 
I=abs(Ri).^2;
clear RR;
end
function [I_O]=Addnoise(I,NoiseLevel,N)
Noisepercentage=NoiseLevel/100;
SNR_image=1/Noisepercentage;
noise = random('norm', 0, Noisepercentage*std2(I), N,N);
I_O=I+noise;
end
function [I]=OTFatten(I,H)
I=real(ifft2(ifftshift(H.*fftshift(fft2(I)))));
end

function [Obj1,Obj2]=Object_Para_Estimation(delI_FT,D,Kotf,U,V,W,K,OTF)
CC=(K>0.3*Kotf).*(K<0.6*Kotf);
NSK=(delI_FT).*CC;
A =sum(sum(abs(NSK)))./sum(sum(CC));
alpha =-0.5;
OBJparam = [A alpha];
Objparaopt=@(OBJparam)Objoptfunc(OBJparam,delI_FT,Kotf,K,D,W,OTF)
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
[Objpara,fval]=fminsearch(Objparaopt,OBJparam,options);
Obj1=Objpara(1);
Obj2=-Objpara(2);
end
function C= Objoptfunc(OBJparam,delI_FT,Kotf,K,D,W,OTF)
sz=size(delI_FT,1)/2;
K(sz+1,sz+1) =1;
A= OBJparam(1,1);
alpha= OBJparam(1,2);
Objamp=A.*((K).^(alpha));
Signalamp=abs(Objamp).*D.*OTF;

Noisespectrum=delI_FT.*W;
NoisePower = sum(sum(abs(Noisespectrum).^2))./sum(sum(W));

Noisefreesignalamp=((abs(delI_FT).^2))-NoisePower;
Noisefreesignalamp=sqrt(Noisefreesignalamp);
Error=(Noisefreesignalamp)-Signalamp;
Zloop = (K<0.8*Kotf).*(K>0.2*Kotf);
invK=1./(K);
C=sum(sum(((Error.^2.*invK).*Zloop)));
end

function [WF_AF1,Noise_AF1]=Weiner_filter_center(delI_FT,A,alpha,D,K,W,OTF)
Noisespectrum=delI_FT.*W;
sz=size(delI_FT,1)/2;
K(sz+1,sz+1) =1;
NoisePower = sum(sum(abs(Noisespectrum).^2))./sum(sum(W));
Otfpower=abs(D).^2;
OBJamp= A.*abs(K).^(-alpha);
OBJpower= OBJamp.^2;
T1= (conj(D)./NoisePower);
T2= Otfpower./NoisePower;
T3= 1./OBJpower;
Filter= T1./(T2+T3);
WF_AF1= Filter.*delI_FT;
Noise_AF1=NoisePower;
end
