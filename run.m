% ------------------------------MCML---------------------------------------
%
% Mo phong Monte Carlo cho da song, da lop
% Ngay: 02/07/2021, NIRS Team
% 
% Cac tep lien quan:
%   Tep du lieu         - Thong so quang hoc cua mo 400nm-1100nm [lamda(nm), mua(1/cm), mus(1/cm), g(1/cm) ])
%                         epi_interp.mat --> bieu bi
%                         der_interp.mat --> ha bi
%                         subf_interp.mat--> duoi da
%                         musc_interp.mat--> mo co
%                         blo_interp.mat --> mau
%   Tep thuc thi        - run.m (khoi tao ma tran 'tissuelayer' va so hat 'N', hien thi)
%   Cac tep chuc nang   - Monte3.m (nhan tissuelayer va N de mo phong Monte Carlo )
%                       - check.m (kiem tra photon truyen qua hay phan xa)
%                       - makec2f.m (colorbar)
%
% Su dung: 
%   run()
%
% Du lieu vao (thay doi trong tep run.m):
%   tissuelayer         - Khoi tao ma tran ca thuoc tinh cho mot hoac nhieu lop mo
%                         mua: he so hap thu [1/cm]
%                         mus: he so tan xa [1/cm]
%                         g: he so di huong [1/cm]
%                         n: chiet suat mo (ma tran chiet suat n=[n1 n2 ...nn])
%                         d: do day mo [cm] (ma tran do day d=[T1 T2 ... Tn])
%                         "tissuelayer=[t1;t2;t3;t4];
%                         tissuelayer=horzcat(tissuelayer,n',d')"
%                         >>tissuelayer=[mua ms g n d; lop 1
%                                        mua ms g n n; lop 2
%                                        mua ms g n n; lop 3
%                                        mua ms g n n] lop 4
%                       - Truong hop mo phong mot lop mo
%                         n=[n1];
%                         d=[d1];
%                         tissuelayer=[ti];
%                         tissuelayer=horzcat(tissuelayer,n,d);
%                         >>tissuelayer=[mua ms g n d]
%
%   lamdaminmax         - Khoang buoc song 
%                         min=400nm
%                         max=1100nm
%                         >>khoi tao ma tran lamda=[400 401 ...1100]
%   N                   - So luong photon mo phong cho tung buoc song
%                         N=10000, lamda=[400 1100]
%                         >>tong photon mo phong= 7000000 hat
%
% Du lieu ra (trong tep run.m):
%   C                   - Ma tran tong mat do nang luong hap thu cho moi buoc song A (J/cm3)
%                         >> ve do thi 2D Mat do nang luong hap thu (J/cm3) o moi buoc song
%   H1                  - Ma tran tong Mat do nang luong hap thu A cho dai song (J/cm3)
%                         >> tao hinh anh Mat do nang luong hap thu log10(A) (J/cm3)
%   H2                  - Ma tran tong Do luu loat F cho dai song (J/cm2)
%                         >> tao hinh anh Do luu thong log10(F) (J/cm2)
%                         >> ve do thi 2D Do lu thong F (J/cm2) theo chieu sau
%
% Mo phong 1000000 hat trong ~117 s
% -------------------------------------------------------------------------

%lay du lieu tu file.mat [lamda(nm), mua(1/cm), mus(1/cm), g(1/cm)]
global V
load('epi_interp.mat'); %lop bieu bi
load('der_interp.mat'); %lop ha bi
load('subf_interp.mat');%lop mo duoi da
load('musc_interp.mat');%lop co
load('blo_interp.mat'); %mau

n=[1.34 1.4 1.44 1.36 1.38]; %chiet suat cho tung lop mo n=[n1 n2 n3 n4]
d=[0.01 0.2 0.6 0.15 5];%do sau cho tung lop mo (cm) d=[T1 T2 T3 T4]
  
N=10000;%so hat cho tung lamda
dr=0.005; dz=0.005;
C=[];
H1=zeros(1200);
H2=zeros(1200);
lamdaminmax=[400 1100];%khoang buoc song min=400nm, max =1100nm
lamda=lamdaminmax(1):1:lamdaminmax(2);
disp('Please wait...');
%the tich cua phan thu luoi
V=1:1:1200;
V=(2*V+1)*pi*dr^2*dz; %(cm3)
tic
%//////////////////////////////////////////////////////////////////////////
for i=1:length(lamda)
    
    tissuelayer=[];
    mark=lamda(i)- musc_interp(1,1)+1;
    t1=epi_interp(mark,2:4); %lay mua, mus, g
    t2=der_interp(mark,2:4);
    t3=subf_interp(mark,2:4);
    t4=blo_interp(mark,2:4);
    t5=musc_interp(mark,2:4);
    
    tissuelayer=[t1;t2;t3;t4;t5]; % tissuelayer=[ti]; truong hop cho 1 lop
    tissuelayer=horzcat(tissuelayer,n',d');%tissuelayer=horzcat(tissuelayer,n,d); truong hop cho 1 lop
    
    [a,b]=Monte3(tissuelayer,N);%thuc hien mo phong
    C(i)=sum(sum(a));%tinh tong matran mat do nang luong hap thu (A)
    H1=H1+a; %ma tran Mat do nang luong hap thu A (J/cm3)
    H2=H2+b; %ma tran Do luu thong F (J/cm2)
end
disp('Done!.');
toc
%//////////////////////////////////////////////////////////////////////////   
x_plot= 0:0.005:(1200-0.005)*0.005;%ghi lai chieu sau (Zmax=1200*0.005 cm)
y_plot=sum(H2);
for i=0:length(y_plot)
    if y_plot(length(y_plot)-i)>=0.0001
        k=length(y_plot)-i;%lay gioi han k (k tai do do luu thong F>=0.0001(J/cm2))
        break;
    end
end
%//////////////////////////////////////////////////////////////////////////
%plot
figure(1)
   subplot(1,2,1)
   imagesc(log10(H1(1:600,1:600)));
   title('Mat do nang luong hap thu log10(A) (J/cm3)');
   xlabel('Do sau Z (cm)');
   ylabel('Ban kinh R (cm)');
   xt = get (gca, 'XTick' );
   yt = get (gca, 'YTick' );
   set (gca, 'XTick' , xt, 'XTickLabel' , xt*0.005);
   set (gca, 'YTick' , yt, 'YTickLabel' , yt*0.005);
   set(gca,'fontsize',12);
   colorbar
   colormap(makec2f);% <---------------makec2f.m
   set(colorbar,'fontsize',12)

   subplot(1,2,2)
   imagesc(log10(H2(1:600,1:600)));
   title('Do luu thong log10(F) (J/cm2)');
   xlabel('Do sau Z (cm)');
   ylabel('Ban kinh R (cm)');
   xt = get (gca, 'XTick' );
   yt = get (gca, 'YTick' );
   set (gca, 'XTick' , xt, 'XTickLabel' , xt*0.005);
   set (gca, 'YTick' , yt, 'YTickLabel' , yt*0.005);
   set(gca,'fontsize',12);
   colorbar
   colormap(makec2f);% <---------------makec2f.m
   set(colorbar,'fontsize',12)

figure(2)
   plot(x_plot(1:k),y_plot(1:k),'ko');
   title('Do luu thong F (J/cm2) theo chieu sau');
   xlabel('Do sau Z (cm)');
   ylabel('Do luu loat F (J/cm2)');
   set(gca,'fontsize',12);
   p= polyfit(x_plot,y_plot,9);
   zi= linspace(0,4);
   hold on
   pz=polyval(p,zi);
   plot(zi,pz,'b','LineWidth',3);
   legend('Duong cong luu thong F thuc te','Duong cong luu thong F min hoa');
 
figure(3)
   plot(lamda,C,'ko');
   title('Mat do nang luong hap thu (J/cm3) o moi buoc song');
   xlabel('Buoc song (nm)');
   ylabel('Mat do nang luong hap thu (J/cm3)');
   set(gca,'fontsize',12);
   v1=polyfit(lamda,C,9);
   vz=polyval(v1,lamda);
   hold on
   plot(lamda,vz,'b','LineWidth',3);
   legend('Duong cong mat do nang luong hap thu thuc te','Duong cong mat do nang luong hap thu min hoa');
