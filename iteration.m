clc;
clear;
p=0;close_C=[];phase_all=[];wrap_all=[]; h=4622; w=4282; F=[];phase_er=[];mask=[];close_B=[];close_F=[];FRR=[];build=[];close_b=[];xx=[];ee=[];LL=[];gai=[];sum_close=[];closeAB=[];
bs=load([pwd '\小基线对组合列表_37.gdh']);
nimg= size(bs,1);
bs=double(bs);
A=bs(:,1);
B=bs(:,2);
bs_name=num2str(bs);
workpath=pwd;
phasepath=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_corphase.gd',nimg,1)];
wrap_phasepath=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_diffF.gd',nimg,1)];
parfor i=1:nimg
    id1=find(A==B(i));
   for j=1:size(id1)
       id2=find(B==B(id1(j)));
    for k=1:size(id2)
       if(A(i)==A(id2(k)))
        p=p+1; 
        f=[i id1(j) id2(k)];
        F=[F;f];
        phase1 = multibandread(phasepath(i,:),[h,w,1],'float',0,'bsq','l');
        phase2 = multibandread(phasepath(id1(j),:),[h,w,1],'float',0,'bsq','l');
        phase3 = multibandread(phasepath(id2(k),:),[h,w,1],'float',0,'bsq','l');
        wrap_phase1 = multibandread(wrap_phasepath(i,:),[h,w,2],'float',0,'bip','l');
        wrap_phase2 = multibandread(wrap_phasepath(id1(j),:),[h,w,2],'float',0,'bip','l');
        wrap_phase3 = multibandread(wrap_phasepath(id2(k),:),[h,w,2],'float',0,'bip','l');
        wrap_phase1 =(complex(wrap_phase1(:,:,1),wrap_phase1(:,:,2)));
        wrap_phase2 =(complex(wrap_phase2(:,:,1),wrap_phase2(:,:,2)));
        wrap_phase3 =(complex(wrap_phase3(:,:,1),wrap_phase3(:,:,2)));    
        close_phase=phase1+phase2-phase3; 
        close_wrap=angle(wrap_phase1.*wrap_phase2 .*conj(wrap_phase3));
        close=round((close_phase-close_wrap)/(2*pi));
        
        close_A=(close(:))';
        close_B=[close_B;close_A];
         break;
      end
    end
   end
end
    closeA=ones(size(F,1),h*w)*9;
    closeA(close_B==0)=1;
    closeA(close_B~=0)=0;
for a=1:nimg
     [row,~]=find(F==a);
     Fr=F(row,:);
     ZFr=size(Fr,1);
     close_AA=zeros(1,h*w);
       for b=1:ZFr
         close_AA=close_AA+closeA(row(b,1),:);
       end
        close_AA=(close_AA)/ZFr;
        closeAB=[closeAB;close_AA];
end
num_blocks=round(size(close_B,2)/2000);
block_size=ceil(size(close_B,2)/num_blocks);
blocks=cell(num_blocks,1);
size_B=size(close_B,2);
closeAB_blocks=cell(num_blocks,1);
for i=1:num_blocks
    start_col=(i-1)*block_size+1;
    end_col=min(i*block_size,size_B);
    blocks{i}=close_B(:,start_col:end_col);
    closeAB_blocks{i}=closeAB(:,start_col:end_col);  
end
clear close_B close_A close_AB  close_AA closeA
for j=1:num_blocks
    close_BB=blocks{j};
    closeABB=closeAB_blocks{j};
    gai=[];
   parfor u=1:size(close_BB,2)
       ture=(ones(nimg,1)*88);
       IDz=find(closeABB(:,u)==1);
       ture(IDz)=0;
       IDt=find(ture==0);
       [xx,~]=find(sum(ismember(F,IDt),3));
       xx=unique(xx);
     for uu=1:size(xx,1)
        if closeABB(F(xx(uu,1),1),u)>0.5 && closeABB(F(xx(uu,1),2),u)>0.5 && closeABB(F(xx(uu,1),3),u)>0.5
            ture(F(xx(uu,1),1),1)=0;
            ture(F(xx(uu,1),2),1)=0;
            ture(F(xx(uu,1),3),1)=0;
        end
     end
       panduan=ones(nimg,1)*9999;
    while  ~isequal(ture,panduan)
          panduan=ture;
        for  q=1:p
            if  (ture(F(q,1),1)~=88) && (ture(F(q,2),1)~=88) && (ture(F(q,3),1)==88)
                ture(F(q,3),1)=-(close_BB(q,u)-ture(F(q,1),1)-ture(F(q,2),1));
            elseif  (ture(F(q,1),1)==88) && (ture(F(q,2),1)~=88) && (ture(F(q,3),1)~=88)
                ture(F(q,1),1)=(close_BB(q,u)-ture(F(q,2),1)+ture(F(q,3),1));
            elseif  (ture(F(q,1),1)~=88) && (ture(F(q,2),1)==88) && (ture(F(q,3),1)~=88)
                ture(F(q,2),1)=(close_BB(q,u)-ture(F(q,1),1)+ture(F(q,3),1));
            end
        end
    end
    gai=[gai ture];
   end
   sum_close=[sum_close gai];
end
headpath1=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_uintP.gdh',nimg,1)];
headpath2=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_ite.gdh',nimg,1)];
headpath3=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_uintP.hdr',nimg,1)];
headpath4=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_ite.hdr',nimg,1)];
name=[repmat(workpath,nimg,1) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('\\',nimg,1) bs_name(:,1:8) repmat('_',nimg,1) bs_name(:,11:18) repmat('_ite.gd',nimg,1)];
parfor n=1:nimg
    phase = multibandread(phasepath(n,:),[h,w,1],'float',0,'bsq','l');
    phase_r=(phase(:))';
    phase_all=[phase_all;phase_r];
end
   parfor gg=1:nimg
    ly=reshape(sum_close(gg,:),h,w);
    wx=reshape(phase_all(gg,:),h,w);
    ly(find(ly(:)==88))=0;
    qq=wx-ly*(2*pi);
    multibandwrite(qq,name(gg,:),'bsq','precision','float');
    oldname1=headpath1(gg,:);
    newname1=headpath2(gg,:);
    oldname2=headpath3(gg,:);
    newname2=headpath4(gg,:);
    copyfile(oldname1,newname1)
    copyfile(oldname2,newname2)     
   end