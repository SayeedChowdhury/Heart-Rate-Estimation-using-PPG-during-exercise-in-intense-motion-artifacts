function [ bpmf,ty,tf,tf2,signal] = bpmgeneration( signal,Fs,step,ha,f,dl,s,bpm,ty,tf,tf2,prev)
if s==1
    
    bpm(1)= firstbpm( signal(1,:),signal(2,:),signal(3,:),signal(4,:),signal(5,:) );
    % Sgolay filtering done with present window values only as input
    for j=1:50
        for t=1:5
            signal(t,:)=sgolayfilt(signal(t,:),3,9);
        end
    end
else
    % Concatenation of 250 new samples of present window with last 750
    % smoothed samples available from previous window for further processing
    % No future value is used anywhere
    signal=[prev(:,1+step:end) signal(:,end-step+1:end)];
    
    x=signal(1,:);
    t1=x;
    x=signal(2,:);
    t2=x;
    d=t1-1.05*t2; % difference signal used for reference noise
    for j=1:50
        for t=1:5
            signal(t,:)=sgolayfilt(signal(t,:),3,9);
        end
    end
    x=signal(1,:);
    
    z=sum(x)/length(x);
    x1=x-z;
    x1=normalize(x1);
    x11=x1(1,:);
    accz=signal(5,:);
    accy=signal(4,:);
    accx=signal(3,:);
    acczn=sum(accz(end:-1:end-11));
    accyn=sum(accy(end:-1:end-11));
    accxn=sum(accx(end:-1:end-11));
    rdirect=accxn+accyn+acczn;
    
    % concatenation of 12 samples from smoothed signal of previous
    % window and 988 smoothed signal of current window
    accxd=[prev(3,1+step-dl:step) signal(3,1:end-dl)];
    accyd=[prev(4,1+step-dl:step) signal(4,1:end-dl)];
    acczd=[prev(5,1+step-dl:step) signal(5,1:end-dl)];
    
    accxp=sum(accxd(1:dl));
    accyp=sum(accyd(1:dl));
    acczp=sum(acczd(1:dl));
    
    rdelay=accxp+accyp+acczp;
    
    if rdirect<rdelay
        
        accz=acczd;
        accy=accyd;
        accx=accxd;
    end
    
    racc=accx.^2+accy.^2+accz.^2;
    z=sum(accx)/length(accx);x1acc=accx-z;
    x1acc=normalize(x1acc);
    [yy,e1e] = filter(ha,x1acc,x11); % MA reduction using adaptive filter
    e1e=e1e-mean2(e1e);
    e1e=normalize(e1e);
    % BPM estimation using Periodogram
    [ya,f] = periodogram(e1e,[],f,Fs);
    ya(1:75)=0;
    ya=normalize(ya);
    
    [ypeak,loc]=findpeaks(ya,'SORTSTR','descend');z=loc(1);
    r(1)=f(z)*60;
    
    z=sum(accy)/length(accy);x2acc=accy-z;
    x2acc=normalize(x2acc);
    [yy,e1e] = filter(ha,x2acc,x11);
    e1e=e1e-mean2(e1e);
    e1e=normalize(e1e);
    [ya,f] = periodogram(e1e,[],f,Fs);ya(1:75)=0;
    ya=normalize(ya);
    [ypeak,loc]=findpeaks(ya,'SORTSTR','descend');z=loc(1);
    r(2)=f(z)*60;
    
    z=sum(accz)/length(accz);x3acc=accz-z;
    x3acc=normalize(x3acc);
    
    z=sum(d)/length(d);xd=d-z;
    xd=normalize(xd);
    
    [yy,e1e] = filter(ha,xd,x11);
    e1e=e1e-mean2(e1e);
    e1e=normalize(e1e);
    [ya,f] = periodogram(e1e,[],f,Fs);ya(1:75)=0;
    if s<4
        ya(220:300)=0;
    end
    ya=normalize(ya);
    
    [ypeak,loc]=findpeaks(ya,'SORTSTR','descend');z=loc(1);
    r(3)=f(z)*60;
    for j=1:3
        a(j)=abs(r(j)-bpm(s-1));
    end
    % finding the estimate closest to previous window's BPM
    gh=0;
    j=1:3;b=find(a(j)==min(a));
    if(length(b)>1)
        for l=1:length(b)
            if(b(l)==3)
                bpm(s)=r(3);gh=1;
            end
        end
        if gh==0
            bpm(s)=r(b(2));
        end
        
        
    else
        bpm(s)=r(b);
    end
    % from 11th window onwards some heuristic approaches are used
    % for tracking verification
    if(s>10&&(min(a)>8)&&(abs(bpm(s)-bpm(s-2))>10)&&(abs(bpm(s)-bpm(s-3))>10)&&(abs(bpm(s)-bpm(s-4))>12)&&(abs(bpm(s)-bpm(s-5))>12))
        [yae,f] = periodogram(racc,[],f,Fs);yae(1:75)=0;
        yaen=normalize(yae);
        u=3;
        [p,l]=findpeaks(ya,'SORTSTR','descend');
        p=p(1:u);l=l(1:u);
        for t=1:length(p)
            p(t)=p(t)/yaen(l(t));
            p(t)=p(t)/(abs(bpm(s-1)-60*f(l(t))));
            
        end
        m=max(p);
        for n=1:length(p)
            if p(n)==m
                z=n;
            end
        end
        r1=f(l(z))*60;
        bpm(s)=r1;
        if((abs(bpm(s)-bpm(s-2))>15)&&(abs(bpm(s)-bpm(s-3))>20)&&(abs(bpm(s)-bpm(s-4))>20)&&(abs(bpm(s)-bpm(s-5))>20)&&(abs(bpm(s)-bpm(s-6))>20)&&(abs(bpm(s)-bpm(s-7))>20)&&(abs(bpm(s)-bpm(s-8))>20)&&(abs(bpm(s)-bpm(s-9))>20)&&(abs(bpm(s)-bpm(s-1))>15)&&(abs(bpm(s)-bpm(s-10))>20)&&ty>=0)
            ty=ty-1;
            bpm(s)=bpm(s-1);
        end
    end
    
    if(s>5&&abs(bpm(s)-bpm(s-1))>9) % from 6th window onwards this check is done,
        % in earlier windows sudden jumps are not prohibited to allow
        % the estimate to settle
        x11=t2; % using PPG2
        [yy,e1e] = filter(ha,x3acc,x11);
        e1e=e1e-mean2(e1e);
        e1e=normalize(e1e);
        [ya,f] = periodogram(e1e,[],f,Fs);
        ya(1:75)=0;ya(220:300)=0; % neglecting values of very low and very high BPM
        ya=normalize(ya);
        
        [ypeakz,locz]=findpeaks(ya,'SORTSTR','descend');z=locz(1);
        jodi=f(z)*60;
        if(abs(jodi-bpm(s-1))<abs(bpm(s)-bpm(s-1)))
            bpm(s)=jodi;
            
        end
        
    end
    if((bpm(s-1)>=100)&&(bpm(s)-bpm(s-1))<=-17&&tf2>=0) % to prevent sudden fall
        % once BPM estimate has reached a certain limit
        bpm(s)=bpm(s-1);tf2=tf2-1;
    end
    % only BPM estimates from previous windows and estimate of current
    % window are used as inputs to Sgolay filter after some time has elapsed.
    % No future value is used
    if s>120
        mj=bpm(1:s);
        for gb=1:200
            mj=sgolayfilt(mj,3,9);
        end
        bpm(s)=mj(s);
    end
    
end


bpmf=bpm(s); % BPM estimate of current window

end


function [ x ] = normalize( x )
x=x/max(x);
end

function [ bpm ] = firstbpm( p1,p2,x,y,z ) % for evaluating first BPM of each dataset

d=p1-p2; % difference signal used for reference noise

r=sqrt(x.^2+y.^2+z.^2);

% Sgolay filtering done with present window values only as input
for j=1:40                              
    p1=sgolayfilt(p1,3,9);
    p2=sgolayfilt(p2,3,9);
    
    x=sgolayfilt(x,3,9);
    y=sgolayfilt(y,3,9);
    z=sgolayfilt(z,3,9);
end



[p1act,P1]=bpmfind(p1,60,200);
[p2act,P2]=bpmfind(p2,60,200 );
[ract,R]=bpmfind(r,60,200 );

if abs(p1act-p2act)<10 && (P1/R>2 || P2/R>2)
    bpmset=[p1act,p2act];
end

N=60;

bpmd=bpm_t(p1,d,N,60,200 );
bpmx=bpm_t(p1,x,N,60,200 );
bpmy=bpm_t(p1,y,N,60,200 );
bpmz=bpm_t(p1,z,N,60,200);

if abs(p1act-p2act)<10 && (P1/R>2 || P2/R>2)
    bpmset=[bpmset bpmd bpmx bpmy bpmz];
else
    bpmset=[bpmd bpmx bpmy bpmz];
end

bpmselect=bpmset;

moderange=10;         % range of histogram
bpmm=round(bpmselect/moderange)*moderange;
b=find(bpmm==mode(bpmm)) ;
bpf=zeros(1,length(bpmselect));
bpf(1,b)=1;
bpms=bpmselect.*bpf;
bpmnz=bpms(find(bpms~=0));

bpm=mode(bpmnz);

end

function [ bpm ] = bpm_t( p1,r,N,minbpm,maxbpm)

P0 = 10*eye(N);  % Initial sqrt correlation matrix inverse
lam = 0.99;            % RLS forgetting factor
ha = adaptfilt.rls(N,lam,P0);

x=p1-mean(p1);
x=x/max(x);


ref=r-mean(r);
ref=ref/max(ref);

[y1,x1] = filter(ha,ref,x);

[bpm X]=bpmfind(x1,minbpm,maxbpm);

end

function [ bpm, X ] = bpmfind( p1,minbpm,maxbpm )
p1=p1-mean(p1);
p1=p1/max(p1);

NFFT=16384;
Fs=125;
x1=p1;
[X1,F]=periodogram(x1,[],NFFT,Fs);
F=60*F;
X1=X1(F>minbpm&F<maxbpm);
F=F((F>minbpm&F<maxbpm));

X=max(X1);
b=find(X1==X);
fm=min(b);
bpm=F(fm);
end


