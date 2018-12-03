function[] = Thermo
  m = 100;
  # table  temp(R) , pressure, Hf, Hg, Sg, Sf, Vf Vg
  T1 = [1440 .9256 280.56 1166.93 1.2358 .6203 .02224 515.1872];
  T3 = [2560 200.87 500.5869 1222.49 .2619 .7324 .0284 2.8936 ];
  T6 = [1680 5.45 325.76 1180.82 1.1983 .649 0.0233 79.805];
  Reheat = [2500 170.68 488.52 1219.47 1.8281 .7277 0.02801 3.3559];
  x1 = 0;
  % x3 =.85:.025:1;
  x3 = 1;
  
  h1 = (1-x1)*T1(1,3) +x1*T1(1,4);
  s1 = (1-x1)*T1(1,6) +x1*T1(1,5);
  v1 = (1-x1)*T1(1,7) +x1*T1(1,8);
  
  s2 = s1;
  h2 = h1 + v1*(T3(1,2)-T1(1,2));
  h3 = (1-x3)*T3(1,3) +x3*T3(1,4);
  s3 = (1-x3)*T3(1,6) +x3*T3(1,5);
  
  s4 = s3;
  x4 = (s4-Reheat(1,6))/(Reheat(1,5)-Reheat(1,6));
  h4 = (1-x4)*Reheat(1,3) + x4*Reheat(1,4);
  p4 = Reheat(1,2);
  
  %step 5 from superheated table
  h5 = 1219.47;
  s5 = 1.0201;
  
  %step 6
  s6 = s5;
  x6 = (s6-T6(1,6))/(T6(1,5)-T6(1,6));
  h6 = (1-x6)*T6(1,3) +x6*T6(1,4);
  
  Qr = m*(h5-h4);
  Qsg = m*(h3-h2);
  Qin = Qr + Qsg;
  Wp = m*(h1-h2);
  Wt = m*((h3-h4)+(h5-h6));
  n = (Wt+Wp)/(Qsg + Qr);
  
  SumQ = (Qsg/T3(1,1))+(Qr/Reheat(1,1));
  SumS = m*(s6-s1);
  Sigma = SumS - SumQ;
  %Sigma is negative so this rankine cycle is impossible
  
endfunction
