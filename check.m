function [sta,Reflectance,tetai]= check(ni, nt, mz)

          tetac=asin(nt/ni); % goc phan xa toan phan (goc khuc xa toi han)
          tetai=acos(abs(mz));
          if tetai<tetac
             tetar=asin(ni*sin(tetai)/nt);%{Snell's law}
             del=tetai-tetar;
             sum=tetai+tetar;
             Reflectance=0.5*(sin(del)/sin(sum))^2+0.5*(tan(del)/tan(sum))^2;
          else
             Reflectance=1;
          end
          if rand<= Reflectance
             sta= 1;%phan xa
          else
             sta= 0;%truyen qua
          end 
end

