function Peq=Col66(fwav,ofname)
% v05 - Reduce DC
% v06 - scal Log + filter
% v07 - overlap
% v08 - gane for high freq - linear
% v09 - palitra
% v20 - norm with win
% v21 - work with Power (not with Amp)
% v22 - Exp for High freq gane
% v30 - matlab ver14b
% v60 - merge rgb+hvs, square picture option
% v60r - rand option
% v64 - pattern 3D
% v66 - zmeika

%TBD filter after frame sum -?

AudioPart = 1; %def 1
profile = 'MPEG-4';%'Uncompressed AVI'; %'MPEG-4'

FFT_len = 512; %512 , 256
Overlap = FFT_len - 256;    % 256 , 0
Samp_per_sec = 44100;
N_framesper_sec = 30;

FFT_Power_EXP = 2; %2 - for Energy (emph colors effect), 1- for Amp (wav effect)
FILTER_AFT_FFT_LEN = 15;%def 15 (smooth signals, reduce spikes)
FILTER_AFT_FFT_TYPE = 2; %1 -linear; 2 - median
FILTER_END_TIME = 5; %def 10
FILTER_END_FREQ = 1; %def 10
NORM_WIN_SIZE = 10; %in frames; def 10 ; for more dynimic changes in britness; if 1- auto level and constant lighting colors
HiFreqExpGane= 1;%Def 1(lenear gane)
PRE_WIN = 0;  %def 0 - no window, 1 - pre win
EFFECT = 0; % 0- Palitra, 1- MandelBrot, 2 - export
COLMAP = 2; % 0- rgb, 1 - hsv, 2- user
NORMALIZE_OUT = 1; %normailze contrast 0 - no norm, 1 -norm
RandomPattern = 0; %0 - not rand, 1 - rand
Pattern_mode = 0; %0 - lines; 1 -Squares
Pattern3D = 0;     %0 - no 3d , 1 - 3d 
Write_File_3D = 1; % 0 -no write, 1 -write
Rotate_3D = 1; %0- no rotate, 1 - rotate

%N_FFTs_per_frame=round(Samp_per_sec/N_framesper_sec/FFT_len); %TBD normilize

% Get wav and len
[Wch1ch2,fs_upd] = audioread(fwav);%Wch1ch2 = wavread(fwav);
Samp_per_sec = fs_upd;
Wch1ch2 = sum(Wch1ch2,2);

Wlench1ch2 = size(Wch1ch2);
Vlen = Wlench1ch2(1,1);
Vlen = fix(Vlen*AudioPart);

%Overlap H
blocksz = FFT_len - Overlap;
N_ffts_ov = ceil(Samp_per_sec/blocksz) ;
Samp_per_sec_ov = N_ffts_ov*blocksz+Overlap;

%Matrix - Samp_per_sec X N seconds
N_Sec=fix(Vlen/Samp_per_sec);

M_Sps_Nsec = buffer(Wch1ch2(1:Samp_per_sec*N_Sec),Samp_per_sec_ov,Samp_per_sec_ov-Samp_per_sec,'nodelay');
%M_Sps_Nsec = reshape(Wch1ch2(1:Samp_per_sec*N_Sec,1),Samp_per_sec,N_Sec);
Wch1ch2 = [];

%Ffts
N_ffts_pers=N_ffts_ov;%fix(Samp_per_sec/FFT_len);
M_Sps_Nsec_FFt = zeros(FFT_len,N_ffts_pers,N_Sec);
for i=1:N_Sec,
  M_Sps_Nsec_FFt(:,:,i) = buffer(M_Sps_Nsec(:,i),FFT_len,Overlap,'nodelay');
end
%M_Sps_Nsec_FFt = reshape(M_Sps_Nsec(1:FFT_len*N_ffts_pers,:),FFT_len,N_ffts_pers,N_Sec);
M_Sps_Nsec = [];

F_mat = single(M_Sps_Nsec_FFt);
M_Sps_Nsec_FFt = [];

if PRE_WIN == 1,
 W = single(hamming(FFT_len));
 Wn = zeros(FFT_len,N_ffts_pers,N_Sec);
 for t = 1:N_Sec, 
     for k=1:N_ffts_pers,
        Wn(:,k,t) = W(:);
     end
 end
 F_mat = F_mat.*Wn;
end

F_mat = fft(F_mat);
A_mat = abs(F_mat(1:FFT_len/2,:,:)).^FFT_Power_EXP;       % v21 - Power 
F_mat = [];

%Full_for_Filt = single(zeros(FFT_len/2,N_ffts_pers*N_Sec));
Full_for_Filt = A_mat(:,:);
if FILTER_AFT_FFT_TYPE == 1, 
    Full_for_Filt=filter(ones(1,FILTER_AFT_FFT_LEN),1,Full_for_Filt,[],2) ;                             
elseif FILTER_AFT_FFT_TYPE == 2,
    Full_for_Filt=single(medfilt1(double(Full_for_Filt),FILTER_AFT_FFT_LEN,[],2));
end
A_mat = reshape(Full_for_Filt,FFT_len/2,N_ffts_pers,N_Sec);
Full_for_Filt = [];
figure(1);mesh(double(A_mat(:,:,1)));

%Set grid for fft frames to video frames
Np=fix(N_ffts_pers/N_framesper_sec);
Rp = rem(N_ffts_pers,N_framesper_sec);
GrN(1:N_framesper_sec) = Np;
GrN(1:Rp) = GrN(1:Rp) + 1;
St(1)=1;En(1)=St(1)+GrN(1)-1;
for d=2:N_framesper_sec,
    St(d)= En(d-1)+1;
    En(d)=St(d)+GrN(d)-1;
end

Ps_m = zeros(FFT_len/2,N_framesper_sec,N_Sec);
for d=1:N_framesper_sec,
   Ps_m(:,d,:)= sum(A_mat(:,St(d):En(d),:),2);
end
figure(2);mesh(double(Ps_m(:,:,1)));

%M_f = reshape(A_mat(:,1:Np*N_framesper_sec,:),FFT_len/2,Np,N_framesper_sec,N_Sec);
%Ps_m=sum(M_f,2);
%Ps_m = squeeze(Ps_m);

%Ps_m = Ps_m/max(max(max(Ps_m))); %norm

%Log scale for channels
N_ch=FFT_len/2;
xx=[log10(1:N_ch)/log10(N_ch)*N_ch];

if COLMAP == 0,  %RGB   
 CM=colormap(jet(FFT_len/2)); 
 %CMZ=size(CM);CCN = CM; 
 CCN = CM(end:-1:1,:); %reverse for red in low, blue in high
 %YY=resample(CCN,FFT_len/2/CMZ(1),1);
 YY = CCN;
elseif COLMAP == 1, %HSV
 CM=hsv(FFT_len/2);    
 YY = CM;
else  %user      
 CM_rgb=colormap(jet(FFT_len/2));    
 CCN = CM_rgb(end:-1:1,:);
 YY = CCN;
 
 CM_hsv=hsv(FFT_len/2);
 
 len = length(CM_hsv(:,1));
 insert_vec  = [(len - len/4):(len-len/8)];
 YY(insert_vec+len/8,:) = CM_hsv(insert_vec,:);
 
 if RandomPattern == 0,
  RD = 1:len;
 else
  RD=round(rand([len,1])*len+1);
 end
 RD(find(RD >256))=256;
 YY = YY(RD,:);
end


YY(find(YY < 0))=0;

R = zeros(FFT_len/2,1);G = R;B = R;
R = YY(:,1); G = YY(:,2);B = YY(:,3);
%bl=fix(FFT_len/2/3);
%R(1:bl)=1;
%G(bl+1:bl*2) = 1;
%B(bl*2 + 1:end) =  1;

%all sum - true color rgb mxnx3 (0-1 each for double pres)

% Gane high freq v08,v22 
%Amp = (1:FFT_len/2).^HiFreqExpGane;
%H_ind = fix(FFT_len/2 - FFT_len/2/3 ); lnh=length(Amp(H_ind :end));
%Amp(H_ind :end)  = Amp(H_ind :end).*([1:lnh ].^HiFreqExpGane);
%Amp = Amp(:); 

k=1;
 for t=1:N_Sec,
      for i=1:N_framesper_sec
          Ps_m(:,i,t)=interp1(xx,Ps_m(:,i,t),1:N_ch);                     %change scale to log
          %Ps_m(1:3,i,t) = Ps_m(1:3,i,t)*0.1;                                            %reduce first filters     
          %Ps_m(:,i,t)  = Amp.*Ps_m(:,i,t);
          %Ps_m(:,i,t) = Ps_m(:,i,t)/max(max(max(Ps_m(:,i,t))));          %mormalize
          
          Peq(:,k)  = Ps_m(:,i,t);
          k=k+1;
      end
 end
  
 
 Pfrsz = k-1;
 %Equalize  beetween freq
  vm=mean(Peq,2);
  for i=1:Pfrsz,
   Peq(:,i) = Peq(:,i)./vm;
  end
 
   %last filter
  Peq=medfilt1(Peq,FILTER_END_TIME,[],2); 
  Peq=medfilt1(Peq,FILTER_END_FREQ,[],1); 
 
  
  %  %v20- norn by win - auto gain
  wl = NORM_WIN_SIZE;
 for k1= 1+wl:Pfrsz-wl,
  MX = max(max(Peq(:,k1-wl:k1+wl)));
     Peq(:,k1)  = Peq(:,k1)/ MX;
  end
 Peq(find(Peq > 1))=1; %limit values in start an in end
  
  
  if EFFECT == 0,  
    %set params and open file   
    vidObj = VideoWriter(ofname,profile);  %aviobj = avifile(ofname,'FPS',25,'COMPRESSION',Compession,'QUALITY',100);
    vidObj.FrameRate = N_framesper_sec;
    fl=strcmp(profile,'Uncompressed AVI'),
    if fl ~= 1,
     vidObj.Quality = 100;
    end
    open(vidObj);   
    
    %vector for matrix maltip and matrix building
    factor_gap = 1;
    width_factor=2; 
    pattern = Pattern_mode;
    Fr0 = ones(1,FFT_len/width_factor);
    
    %for RGB 
    %colormap
    frame = zeros(FFT_len/2*factor_gap,FFT_len/width_factor,3); %p(:,:,1) =a1;
    
    %write file
    gaps = zeros(length(Peq(:,1))*factor_gap,1);
    NS = length(gaps);
    N = sqrt(NS);
    
    if Pattern3D == 1,
     fract_fig_h=figure('Color',[0 0 0],'Position',[200,100,1280,720]);
     
      vid_frames_cnt =0;
      Vid_part_ind = 1;
      if Write_File_3D == 1,
       %aviobj = avifile([ofname int2str(Vid_part_ind)],'FPS',N_framesper_sec,'COMPRESSION','None','QUALITY',100);
       vid3dObj = VideoWriter([ofname int2str(Vid_part_ind)],profile);
       vid3dObj.FrameRate = N_framesper_sec;
       fll=strcmp(profile,'Uncompressed AVI'),
       if fll ~= 1,
        vid3dObj.Quality = 100;
       end
       open(vid3dObj);  
      end
    end
    
    for k1=1:Pfrsz, 
      Signal = Peq(RD,k1);
        
      frame(:,:,1) =Pattern_Build(Signal,R,factor_gap,pattern,Fr0,NS,N,gaps);
      frame(:,:,2) =Pattern_Build(Signal,G,factor_gap,pattern,Fr0,NS,N,gaps);
      frame(:,:,3) =Pattern_Build(Signal,B,factor_gap,pattern,Fr0,NS,N,gaps);
         
      if NORMALIZE_OUT == 1,
        maxval=max(max(max(frame)));  %normalize  
        frame = frame/maxval;
      end
            
      writeVideo(vidObj,double(frame));%aviobj = addframe(aviobj,double(frame));
         
      if Pattern3D == 1,
       %3d
       Npx = 1280;
       Npy = 720;
       Video_sz_frames = 700; %default 700
       
       color1 = sum(frame,3);
       maxval=max(max(max(frame)));
       color1 = color1/maxval;
         
       figure(fract_fig_h);h=surf(color1);

       cc = k1
       if Rotate_3D  == 1,
        view(-20+cc/2.5,40); %for rotation
       end
       zoom(1.2);
       axis('off');
       axis([1 256 1 256  0 2]);
       axis vis3d
       colormap(YY);
       shading interp;
       camlight left;camlight;lighting gouraud 
       
       vid_frames_cnt = vid_frames_cnt + 1;
       if Write_File_3D == 1,  
        F = getframe(gcf);   
        writeVideo(vid3dObj,F);%aviobj = addframe(aviobj,F);
   
        if vid_frames_cnt >= Video_sz_frames,
         close(vid3dObj);
         vid_frames_cnt = 0;
         Vid_part_ind = Vid_part_ind + 1;
         %aviobj = avifile([ofname int2str(Vid_part_ind)],'FPS',N_framesper_sec,'COMPRESSION','None','QUALITY',100);
         vid3dObj = VideoWriter([ofname int2str(Vid_part_ind)],profile);
         vid3dObj.FrameRate = N_framesper_sec;
         fll=strcmp(profile,'Uncompressed AVI'),
         if fll ~= 1,
          vid3dObj.Quality = 100;
         end
         open(vid3dObj);
        end
       end
      end 
    end %frame loop
    %close file
    close(vidObj);%aviobj = close(aviobj); 
    
    if Pattern3D == 1,
     if Write_File_3D == 1,
      close(vid3dObj);
     end
    end
      
  elseif EFFECT == 1,
     mandel06(ofname,double(Peq));
     
  else  %plot square
   
          
  end
  
Peq = double(Peq);
     
function frameT=Pattern_Build(Signal,ColorU,factor_gap,pattern,Fr0,NS,N,gaps)
   
       zmeika_Order =0;
       
       temp=(Signal.*ColorU);
         
         for g=1:factor_gap,
          gaps(g:factor_gap:end) = temp; %insert gaps
         end
         
         if pattern ==0,%lines
          frameT = gaps*Fr0;
         else      %squares
          vT=[];
          mT=zeros(NS,NS);
          for q=1:N,
              v0 = gaps((1:N)+(q-1)*N); 
              if mod(g,2)==0 && zmeika_Order == 1
               v0 = v0(N:-1:1);
              end
              
              m1 = v0*Fr0(1:N);
              m11 = m1';
              v11 = m11(:);
              v22 = [];
              for v=1:N,
                  v22 = [v22 ; v11];
              end
              vT = [vT ; v22];
          end
          mT(:)=vT;
          frameT = mT;          
         end  

