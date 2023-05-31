clc
close all
clear all;
NA=0.42;
lambda=660*10^-9;
npad=0;
N=1200 + 2*npad;
no=N/2;
Mag=20;
dpix=(4.4*10^-6)/Mag;
NoiseLevel=5;
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

% figure
% hold on
% plot(u,abs(H(no+1,:)),'--','LineWidth',2,'MarkerSize',6)
% plot(u,abs(G(no+1,:)),'--','LineWidth',2,'MarkerSize',6)
% legend('H','G')
% title("OTF approximation")
% xline(Kotf)
% xlabel("u");
% grid on
% box on

OTF=H;

Imagee=double(imread("gold.png"));
Imagee = padarray(Imagee,[npad npad],0,"both");

P= double(Imagee);
P=rescale(P,-0.2,0.3);
T=ones(N);
Object=sqrt(T).*exp(1j.*P);

figure
imagesc(x*10^3,y*10^3,(P));
title("Ground truth")
colormap(gray)
colorbar

prop_kernel  = exp(-1j*pi*lambda*z.*(U.^2+V.^2));


fr=600*10^3;
theta_ar= [0 pi/3 2*pi/3];
phi1=0;
phi2=2*pi/3;
phi3=4*pi/3;
P_c=0;


for i =1:numel(theta_ar)

theta=theta_ar(i);
phi= phi1;
I_ill=  (1/2)*(1 +  cos((2*pi*(fr)*cos(theta).*X) +  (2*pi*(fr)*sin(theta).*Y)  + phi));
U_ill=sqrt(I_ill);
MObject=Object.*U_ill;
Ip = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*prop_kernel))));
[Ip]=Addnoise(Ip,NoiseLevel,N);
Im = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*conj(prop_kernel)))));
[Im]=Addnoise(Im,NoiseLevel,N);
delI_FT1 = fftshift(fft2(Ip)) - fftshift(fft2(Im));

phi= phi2;
I_ill=  (1/2)*(1 +  cos((2*pi*(fr)*cos(theta).*X) +  (2*pi*(fr)*sin(theta).*Y)  + phi));
U_ill=sqrt(I_ill);
MObject=Object.*U_ill;
Ip = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*prop_kernel))));
[Ip]=Addnoise(Ip,NoiseLevel,N);
Im = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*conj(prop_kernel)))));
[Im]=Addnoise(Im,NoiseLevel,N);
delI_FT2 = fftshift(fft2(Ip)) - fftshift(fft2(Im));

phi= phi3;
I_ill=  (1/2)*(1 +  cos((2*pi*(fr)*cos(theta).*X) +  (2*pi*(fr)*sin(theta).*Y)  + phi));
U_ill=sqrt(I_ill);
MObject=Object.*U_ill;
Ip = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*prop_kernel))));
[Ip]=Addnoise(Ip,NoiseLevel,N);
Im = real(ifft2(ifftshift((fftshift(fft2(MObject)).*OTF.*conj(prop_kernel)))));
[Im]=Addnoise(Im,NoiseLevel,N);
delI_FT3 = fftshift(fft2(Ip)) - fftshift(fft2(Im));

M=[1/2, (1/4)*exp(1j.*phi1), (1/4)*exp(-1j.*phi1); 1/2, (1/4)*exp(1j.*phi2), (1/4)*exp(-1j.*phi2); 1/2, (1/4)*exp(1j.*phi3), (1/4)*exp(-1j.*phi3)];
IM=inv(M);

P_c= P_c + IM(1,1).*delI_FT1 +  IM(1,2).*delI_FT2 + IM(1,3).*delI_FT3;
%P_c= IM(1,1).*delI_FT1 +  IM(1,2).*delI_FT2 + IM(1,3).*delI_FT3;


P_c(isinf(P_c))=1;
P_c(isnan(P_c))=1;

P_m(:,:,i)=IM(2,1).*delI_FT1 +  IM(2,2).*delI_FT2 + IM(2,3).*delI_FT3;
P_p(:,:,i)=IM(3,1).*delI_FT1 +  IM(3,2).*delI_FT2 + IM(3,3).*delI_FT3;

   
end

P_c=P_c/numel(theta_ar);

D=(4.*sin(pi*lambda*z.*(U.^2+V.^2))).*OTF;
D(isnan(D))=1;
D(isinf(D))=1;

%%%window 
W=zeros(N,N);
W=K>Kcut*Kotf;
W=double(W);

[Al,al]=Object_Para_Estimation(P_c,D,Kotf,U,V,W,K,OTF);
Objamp = Al.*(K.^-al);
Sigamp = Objamp.*D;

% figure
% hold on
% plot(u*10^-3,log(abs(P_c(no+4,:))),'k--')
% plot(u*10^-3,log(abs(Sigamp(no+4,:))),'m-')
% xlabel("u(cycles/mm)");
% legend('Acutal signal power','Estimated. signal power')
% title("Power spectrum estimation")
% grid on
% box on
%% 

[P_c,Noise_P_c]=Weiner_filter_center(P_c,Al,al,D,K,W,OTF);
P_c_s=real(ifft2(ifftshift(P_c)));


OBJamp1=Al.*(K).^(-al);
Signamp1=OBJamp1.*OTF;
SNR_Centre=Signamp1.*conj(Signamp1)./Noise_P_c;
Fmerged_Num= SNR_Centre.*P_c;
Fmerged_Den= eps + SNR_Centre;
Apw= K>Kotf;
Index =0.2;


for i =1:numel(theta_ar)

theta=theta_ar(i);
Km=sqrt( (U-(fr*cos(theta))).^2 +  (V-(fr*sin(theta))).^2 );
Kp=sqrt( (U+(fr*cos(theta))).^2 +  (V+(fr*sin(theta))).^2 );    

[P_m(:,:,i),N_P_m]=Weiner_filter_side(P_m(:,:,i),U,V,Al,al,D,Kotf,Km,W);
[P_p(:,:,i),N_P_p]=Weiner_filter_side(P_p(:,:,i),U,V,Al,al,D,Kotf,Kp,W);

fvector=[fr*cos(theta) fr*sin(theta)];
[minValue,uclosestIndex] = min(abs(u-fvector(1)));
[minValue,vclosestIndex] = min(abs(v-fvector(2)));
up=uclosestIndex-(no+1);
vp=vclosestIndex-(no+1);
P_m(:,:,i)=circshift(P_m(:,:,i),-up,2);
P_m(:,:,i)=circshift(P_m(:,:,i),-vp,1);
P_p(:,:,i)=circshift(P_p(:,:,i),+up,2);
P_p(:,:,i)=circshift(P_p(:,:,i),+vp,1);

OBJamp_m=Al.*abs(Km).^(-al);
OBJamp_p=Al.*abs(Kp).^(-al);
Signamp_m=OBJamp_m.*OTF;
Signamp_p=OBJamp_p.*OTF;
  
Signamp_m=circshift(Signamp_m,-up,2);
Signamp_m=circshift(Signamp_m,-vp,1);
Signamp_p=circshift(Signamp_p,up,2);
Signamp_p=circshift(Signamp_p,vp,1);

SNR_m=Signamp_m.*conj(Signamp_m)./N_P_m;
SNR_p=Signamp_p.*conj(Signamp_p)./N_P_p;
   
Fmerged_Num= Fmerged_Num  + SNR_p.*P_p(:,:,i)  + SNR_m.* P_m(:,:,i);
Fmerged_Den= Fmerged_Den +  SNR_p + SNR_m;

Wm=double(Km>Kotf);
Wp=double(Kp>Kotf); 

Apw=Apw.*Wm.*Wp;

end

%% 
DistApoMask = bwdist(Apw);
maxApoMask = max(max(DistApoMask));
ApoFunc = double(DistApoMask./maxApoMask).^Index;


Fmerged=Fmerged_Num./Fmerged_Den;
Fmerged(isnan(Fmerged))=1;
Fmerged(isinf(Fmerged))=1;

FFF=Fmerged.*ApoFunc;
p = 10;
minL = min(min( abs(FFF).^(1/p) ));
maxL = max(max( abs(FFF).^(1/p) ));
figure;
imagesc(u*10^-3,v*10^-3,abs(FFF).^(1/p),[minL maxL])
xlabel("u(cycles/mm)");
ylabel("v(cycles/mm)");
title(" Merged frequency spectrum")
colormap(gray);
colorbar

SIMimage=real( ifft2(ifftshift(Fmerged)));

figure
imagesc(x*10^3,y*10^3,(SIMimage));
xlabel("x(mm)");
ylabel("y(mm)");
title("High resolution phase using Strucutred illumination")
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
disp(10 * log10(SNR_image));

noise = random('norm', 0, Noisepercentage*std2(I), N,N);
I_O=I+noise;
end
function [I]=OTFatten(I,H)
I=real(ifft2(ifftshift(H.*fftshift(fft2(I)))));
end

function [Obj1,Obj2]=Object_Para_Estimation(delI_FT,D,Kotf,U,V,W,K,OTF)
CC=(K>0.2*Kotf).*(K<0.8*Kotf);
NSK=(delI_FT).*CC;
A =sum(sum(abs(NSK)))./sum(sum(CC));
alpha =-0.5;
OBJparam = [A alpha];
Objparaopt=@(OBJparam)Objoptfunc(OBJparam,delI_FT,Kotf,K,D,W,OTF)
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',200,'MaxIter',200,'Display','notify');
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
Signalamp=abs(Objamp).*D;
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

function [WF_AF1,Noise_AF1]=Weiner_filter_side(AF1,U,V,A,alpha,Otf,Kotf,shifted_K,W)

element = 0; 
[~,ind] = min(reshape(abs(bsxfun(@minus,shifted_K,element)),numel(shifted_K),[]));
[x,y] = ind2sub(size(shifted_K),ind);


shifted_K(x,y)=1;


Noisespectrum=AF1.*W;
NoisePower = sum(sum(abs(Noisespectrum).^2))./sum(sum(W));
Otfpower=abs(Otf).^2;
OBJamp= A.*(shifted_K).^(-alpha);
OBJpower= OBJamp.^2;
T1= (conj(Otf)./NoisePower);
T2= Otfpower./NoisePower;
T3= 1./OBJpower;
Filter= T1./(T2+T3);
WF_AF1= Filter.*AF1;
Noise_AF1=NoisePower;
end
