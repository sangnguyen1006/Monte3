function [H1,H2]=Monte3(tissuelayer,N)
global V
% Nguong dung photon
Wthr=0.001;
% ban kinh chum tia
spot=0.15;%cm
% nguong xo so
m=10; %{nm}
% Mo phong don gian cho n lop vat lieu ----------------------------
c1=1200;c2=1200;% so o phan chia theo r va z 
dr=0.005;dz=0.005; % kich thuoc tang dan moi o la 0.005 cm
%//////////////////////////////////////////////////////////////////////////
n0=1;%{chiet suat khong khi}
Rsp=((n0-tissuelayer(1,4))/(n0+tissuelayer(1,4)))^2;% he so phan xa cua song toi tinh W
%//////////////////////////////////////////////////////////////////////////
Rdif=Rsp; survive=0; zam=0;
siz= size(tissuelayer);
Q1(c1,c2)=0; Q2(c1,c2)=0;% Khoi dau ma tran bang 0  
% Bat dau m o phong N photon
%//////////////////////////////////////////////////////////////////////////
for nn=1:N
    loca=1; %ban dau cho vi tri mo la 1 (lop mo dau tien)
    W=1-Rsp; %trong luong photon (1-Rsp:he so truyen qua)
    %khoi tao vi tri va huong di ban dau
    x1=spot*sqrt(-log(rand)/2); %random vi tri ban dau theo spot
    y1=0;
    z1=0;
    mx=0;my=0;mz=1; %d(0,0,1) vec to don vi di thang goc
%//////////////////////////////////////////////////////////////////////////
    % Bat dau m o phong 1 photon
    while W ~= 0 
        %luu vi tri cu (x1_bef, y1_bef, z1_bef)
        x1_bef=x1;
        y1_bef=y1;
        z1_bef=z1;
        step=-log(rand)/(tissuelayer(loca,1)+tissuelayer(loca,2));% mo phong quang duong tu do S
        x1=x1+mx*step; %x1,y1,z1: vi tri moi cua photon (r=r+d*s),d la vec to huong di moi
        y1=y1+my*step; %mx,my,mz: huong di moi cua photon
        z1=z1+mz*step; %vi tri moi = vi tri cu + huong*quang duong tu do
%//////////////////////////////////////////////////////////////////////////
        %tinh do chenh lech delta (giua z1 va do sau) 
        if loca==1
            delta=z1;
            depth=0;
        else
            depth=0;
            for i1=1:loca-1
                depth=depth+tissuelayer(i1,5);
            end
            delta= z1-depth;
        end
        
        ni= tissuelayer(loca,4); %chiet suat mo lop hien tai 
%//////////////////////////////////////////////////////////////////////////
        % photon di chuyen sang lop mo tren
        if delta<=0
            if loca==1
               % photon thoat khoi mo thu thu nhat ra khong khi
               z1=-z1;
               zam=zam+1;
               ni=tissuelayer(loca,4);
               nt=n0;
               [sta,Reflectance,tetai]= check(ni, nt, mz);
               Rdif=Rdif+(1-Reflectance)*W;
               W=W*Reflectance;
            else
               % photon di chuyen sang lop mo tren khong phai la khong khi
               nt=tissuelayer(loca-1,4);  
               if (ni>nt)&&(abs(mz)~=1) % chiet suat n2<n1, goc toi tetai khac 0
                  [sta,Reflectance,tetai]= check(ni, nt, mz);
                  if sta==1 %phan xa
                     z1=z1+(tissuelayer(loca-1,5)-z1)*2;
                  else 
                        % truyen qua, thay doi huong di cua photon trong moi truong moi
                        ratio=abs(z1_bef-depth)/abs(z1-depth);
                        s1=(ratio*step)/(1+ratio);
                        x1=x1_bef+mx*s1; 
                        y1=y1_bef+my*s1; 
                        z1=z1_bef+mz*s1;
                        s2=(step-s1)*((tissuelayer(loca,1)+tissuelayer(loca,2))/(tissuelayer(loca-1,1)+tissuelayer(loca-1,2)));
                        if mz>0
                            sign=1;
                        else if mz<0
                                sign=-1;
                            else
                                sign=0;
                            end
                        end
                        tetar=asin(ni*sin(tetai)/nt);
                        mx=(mx*tissuelayer(loca,4))/tissuelayer(loca-1,4);
                        my=(my*tissuelayer(loca,4))/tissuelayer(loca-1,4);
                        mz=sign*cos(tetar);
                        x1=x1+mx*s2; 
                        y1=y1+my*s2; 
                        z1=z1+mz*s2; 
                
                        loca=loca-1;
                        ni= tissuelayer(loca,4);
                  end
               else
                  if (ni>nt)&&(abs(mz)==1) % chiet suat n2<n1, goc toi tetai bang 0
                      loca=loca-1;
                      ni= tissuelayer(loca,4);
                  end 
               end
               if (ni<nt)&&(abs(mz)~=1)% chiet suat n2>n1, goc toi tetai khac 0
                   ratio=abs(z1_bef-depth)/abs(z1-depth);
                   s1=(ratio*step)/(1+ratio);
                   x1=x1_bef+mx*s1; 
                   y1=y1_bef+my*s1; 
                   z1=z1_bef+mz*s1;
                   s2=(step-s1)*((tissuelayer(loca,1)+tissuelayer(loca,2))/(tissuelayer(loca-1,1)+tissuelayer(loca-1,2)));
                   if mz>0
                        sign=1;
                    else if mz<0
                            sign=-1;
                        else
                            sign=0;
                        end
                   end
                   tetar=asin(ni*sin(tetai)/nt);
                   mx=(mx*tissuelayer(loca,4))/tissuelayer(loca-1,4);
                   my=(my*tissuelayer(loca,4))/tissuelayer(loca-1,4);
                   mz=sign*cos(tetar);
                   x1=x1+mx*s2; 
                   y1=y1+my*s2; 
                   z1=z1+mz*s2; 
                
                   loca=loca-1;
                   ni= tissuelayer(loca,4);
               else
                   if (ni<nt)&&(abs(mz)==1)% chiet suat n2>n1, goc toi tetai bang 0
                       loca=loca-1;
                       ni= tissuelayer(loca,4);
                   end
               end
            end
        end
 %/////////////////////////////////////////////////////////////////////////
        % photon di chuyen sang lop mo duoi
        if delta>= tissuelayer(loca,5) 
           if loca>=siz(1)
               % photon thoat khoi mo cuoi cung ra khong khi
               nt=n0;
               zam=zam+1;
               [sta,Reflectance,tetai]= check(ni, nt, mz);
               z1=2*( depth+tissuelayer(loca,5))-z1;
               Rdif=Rdif+(1-Reflectance)*W;
               W=W*Reflectance;
           else
               % photon di chuyen sang lop mo duoi khong phai la khong khi
               nt= tissuelayer(loca+1,4);
               if (ni>nt)&&(abs(mz)~=1)% chiet suat n2<n1, goc toi tetai khac 0
                  [sta,Reflectance,tetai]= check(ni, nt, mz);
                  if sta==1 %phan xa
                     z1=2*(depth+tissuelayer(loca,5))-z1;
                  else
                     % truyen qua, thay doi huong di cua photon trong moi truong moi
                     ttd=depth+tissuelayer(loca,5);
                     ratio=abs(z1_bef-ttd)/abs(z1-ttd);
                     s1=(ratio*step)/(1+ratio);
                     x1=x1_bef+mx*s1; 
                     y1=y1_bef+my*s1; 
                     z1=z1_bef+mz*s1;
                     s2=(step-s1)*((tissuelayer(loca,1)+tissuelayer(loca,2))/(tissuelayer(loca+1,1)+tissuelayer(loca+1,2)));
                     if mz>0
                         sign=1;
                     else if mz<0
                             sign=-1;
                         else
                             sign=0;
                         end
                     end
                     tetar=asin(ni*sin(tetai)/nt);
                     mx=(mx*tissuelayer(loca,4))/tissuelayer(loca+1,4);
                     my=(my*tissuelayer(loca,4))/tissuelayer(loca+1,4);
                     mz=sign*cos(tetar);
                     x1=x1+mx*s2; 
                     y1=y1+my*s2; 
                     z1=z1+mz*s2;
                
                     loca=loca+1;
                     ni= tissuelayer(loca,4);
                  
                  end
               else
                   if (ni>nt)&&(abs(mz)==1)% chiet suat n2<n1, goc toi tetai bang 0
                      loca=loca+1;
                      ni= tissuelayer(loca,4);
                   end
               end  
               if (ni<nt)&&(abs(mz)~=1)% chiet suat n2>n1, goc toi tetai khac 0
                   tetai=acos(abs(mz));
                   ttd=depth+tissuelayer(loca,5);
                   ratio=abs(z1_bef-ttd)/abs(z1-ttd);
                   s1=(ratio*step)/(1+ratio);
                   x1=x1_bef+mx*s1; 
                   y1=y1_bef+my*s1; 
                   z1=z1_bef+mz*s1;
                   s2=(step-s1)*((tissuelayer(loca,1)+tissuelayer(loca,2))/(tissuelayer(loca+1,1)+tissuelayer(loca+1,2)));
                   if mz>0
                       sign=1;
                   else if mz<0
                           sign=-1;
                       else
                           sign=0;
                       end
                   end
                   tetar=asin(ni*sin(tetai)/nt);
                   mx=(mx*tissuelayer(loca,4))/tissuelayer(loca+1,4);
                   my=(my*tissuelayer(loca,4))/tissuelayer(loca+1,4);
                   mz=sign*cos(tetar);
                   x1=x1+mx*s2; 
                   y1=y1+my*s2; 
                   z1=z1+mz*s2;
                   
                   loca=loca+1;
                   ni= tissuelayer(loca,4);
               else
                   if (ni<nt)&&(abs(mz)==1) % chiet suat n2>n1, goc toi tetai bang 0
                       loca=loca+1;
                       ni= tissuelayer(loca,4);
                   end
               end
           end 
        end
%//////////////////////////////////////////////////////////////////////////
        %vi tri trong ma tran Q(C1,C2)
        ma=tissuelayer(loca,1);
        ms=tissuelayer(loca,2);
        g=tissuelayer(loca,3);
        
        r=sqrt(x1*x1+y1*y1);
        i=round(r/dr+0.5); %ban kinh
        j=round(z1/dz+0.5); %do sau
        if j<=0 
            j=1;
        end
        dQ1=W*ma/(ma+ms); %trong luong da duoc tich luy khi photon den 1 vi tri 
        dQ2=W/(ma+ms);
        W=W*ms/(ma+ms); %trong luong cua photon tiep tuc di chuyen
        if (i<=c1)&&(j<=c2) 
            Q1(i,j)= Q1(i,j)+dQ1; %tich luy trong luong tai 1 vi tri (x1,y1,z1)
            Q2(i,j)= Q2(i,j)+dQ2; 
        end
%//////////////////////////////////////////////////////////////////////////
        %tinh goc tan xa va goc phuong vi
        if g>0
           teta=acos((1+g^2-((1-g^2)/(1-g+2*g*rand))^2)/2/g); %goc tan xa
        else 
            if g==0
                teta= acos(2*rand -1);
            end
        end
        fi=2*pi*rand; %goc phuong vi
        %tinh toan vecto moi (huong di chuyen moi cua photon) sau khi tich
        %luy nang luong dQ va di tiep voi nang luong W 
        if abs(mz)>0.9999 %truong hop dat biet
            mx=sin(teta)*cos(fi);
            my=sin(teta)*sin(fi);
            mz=mz*cos(teta)/abs(mz);
        else
            mx1=mx;
            my1=my;
            mz1=mz;
            %truong hop tong quat
            mx=sin(teta)*(mx1*mz1*cos(fi)-my1*sin(fi))/sqrt(1-mz1*mz1)+mx1*cos(teta);
            my=sin(teta)*(my1*mz1*cos(fi)+mx1*sin(fi))/sqrt(1-mz1*mz1)+my1*cos(teta);
            mz=-sin(teta)*cos(fi)*sqrt(1-mz1*mz1)+mz1*cos(teta);     
        end
%//////////////////////////////////////////////////////////////////////////       
        %cham dut photon 
        if W<Wthr 
            if rand<=1/m %1 co hoi song sot trong m co hoi duoc trao
                W=m*W;
                survive=survive+1;
            else
                W=0;
            end
        end
    end
end
%CHUYEN NANG LUONG
%Sau khi thuc hien chuong trinh C:\MATLAB\bin\monte.m voi N photon, 
%trong so duoc ghi trong phan tu luoi Q(i,j).
%chuong trinh nay chuyen thanh mat do nang luong voi don vi J/cm3 
%voi gia thiet nang luong cua N photon la 1J
H1=zeros(1200,1200);
H2=zeros(1200,1200);
for i=1:1200
   H1(i,:)=Q1(i,:)/N/V(i); %mat do nang luong tich luy J/cm3
   H2(i,:)=Q2(i,:)/N/V(i); %su phan bo cua anh sang (do luu loat cua anh sang) J/cm2
end

end




