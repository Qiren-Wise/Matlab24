function [intx,intf] = DividePlane(A,c,b,baseVector)

% ���ø�ʽ��[intx,intf] = DividePlane(A,c,b,baseVector)
% ���У�A:Լ������
%        c:Ŀ�꺯��ϵ������
%        b:Լ���Ҷ�����
%        baseVector:��ʼ������
%        intx��Ŀ�꺯��ȡ��Сֵʱ���Ա���ֵ
%        intf��Ŀ�꺯������Сֵ

% ������������²���
%  A=[0.01 0.01 0.01 0.03 0.03 0.03 1 0 0 0;  %Լ������ ϵ��
%       0.02 0 0 0.05 0 0 0 1 0 0;
%       0 0.02 0 0 0.05 0 0 0 1 0;
%       0 0 0.03 0 0 0.08 0 0 0 1];
% %  c=[-20 -14 -16 36 -32 -30];   %���ֵȡ������Сֵ����ϵ��
%  b=[850 700 100 900];            %Լ������ �Ҷ�ֵ
%  [intx,intf]=DividePlane(A,c,b,[7 8 9 10]);  % X7 X8 X9 X10

sz = size(A);
nVia = sz(2);%��ȡ�ж��پ��߱���
n = sz(1);%��ȡ�ж���Լ������
xx = 1 : nVia;

if length(baseVector) ~= n
    disp('�������ĸ���Ҫ��Լ�������������ȣ�');
    mx = NaN;
    mf = NaN;
    return;
end

M = 0;
sigma = -[transpose(c) zeros(1,(nVia - length(c)))];
xb = b;

%���ȵ����α�������Ž�
while 1
    [maxs,ind] = max(sigma);
    %--------------�õ����η������Ž�-------------
    if maxs <= 0 %��У����С��0ʱ��������Ž�
        vr = find(c~=0,1,'last');
        for l = 1 : vr
            ele = find(baseVector == l,1);
            if (isempty(ele))
                mx(l) = 0;
            else
                mx(l) = xb(ele);
            end
        end
        if max(abs(round(mx) - mx)) < 1.0e-7  %�ж����Ž��Ƿ�Ϊ������ �������������
            intx = mx;
            intf = mx*c;
        else %������Žⲻ��������ʱ�������и��
            sz = size(A);
            sr = sz(1);
            sc = sz(2);
            [max_x,index_x] = max(abs(round(mx) - mx));
            [isB,num] = find(index_x == baseVector);
            fi = xb(num) - floor(xb(num));
            for i = 1 : (index_x - 1)
                Atmp(1,i) = A(num,i) - floor(A(num,i));
            end
            for i = (index_x + 1) : sc
                Atmp(1,i) = A(num,i) - floor(A(num,i));
            end
            Atemp(1,index_x) = 0; %������ż�����η��ĳ�ʼ���
            A = [A zeros(sr,1):-Atmp(1,:) 1];
            xb = [xb:-fi];
            baseVector = [baseVetor sc+1];
            sigma = [sigma 0];
            
            %----------��ż�����η��ĵ�������--------------
            while 1
                %-----------------------------------------
                if xb > 0   %�ж�����Ҷ���������0��������Ž�
                  if max(abs(round(xb)-xb)) < 1.0e-7  %����ö�ż�����η�����������⣬�򷵻�����������
                      vr = find(c~=0,1,'last');
                      for l = 1 : vr
                          ele = find(baseVector == l,1);
                          if (isempty(ele))
                              mx_1(1) = 0;
                          else
                              mx_1(1) = xb(ele);
                          end
                      end
                      intx = mx_1;
                      intf = mx_1*c;
                      return;
                  else %�����ż�����η���õ����Žⲻ�������⣬�и��
                      sz = size(A);
                      ar = sz(1);
                      ac = sz(2);
                      [max_x,index_x] = max(abs(round(mx_1) - mx_1));
                      [isB,num] = find(index_x == baseVector);
                      fi = xb(num) - floor(xb(num));
                      for i = 1 : (index_x - 1)
                          Atemp(1,i) = A(num,i) - floor(A(num,i));
                      end
                      for i = (index_x+1):sc
                          Atemp(1,i) = A(num,i) - floor(A(num,i));
                      end
                      Atemp(1,index_x) = 0;%
                      A = [A zeros(sr,1);-Atemp(1,:) 1];
                      xb = [xb;-fi];
                      baseVector = [baseVector sc+1];
                      sigma = [sigma 0];
                      continue;
                  end
                else %����Ҷ�������ȫ����0������ж�ż�����η��Ļ�����������
                    minb_1 = inf;
                    chagB_1 = inf;
                    sA = size(A);
                    [br,idb] = min(xb);
                    for j = 1 : sA(2)
                        if A(idb,j) < 0
                            bn = sigma(j)/A(idb,j);
                            if bn < minb_1
                                minb_1 = bn;
                                chagB_1 = j;
                                
                            end
                        end
                    end
                    sigma = sigma - A(idb,:)*minb_1;
                    xb(idb) = xb(idb)/A(idb,chagB_1);
                    A(idb,:) = A(idb,:)/A(idb,chagB_1);
                    for i = 1:sA(1)
                        if i~= idb
                            xb(i) = xb(i)-A(i,chagB_1)*xb(idb);
                            A(i,:) = A(i,:) - A(i,chagB_1)*A(idb,:);
                        end
                    end
                    baseVector(idb) = chagB_1;
                end
                
            end
            
        end      %��������
%     end
    
else %���У�����в�С��0�ģ�����е������㷨��������
        minb = inf;
        chagB = inf;
        for j = 1:n
            if A(j,ind) > 0
                bz = xb(j)/A(j,ind);
                if bz < minb
                    minb = bz;
                    chagB = j;
                end
            end
            
        end
        sigma = sigma - A(chagB,:)*maxs/A(chagB,ind);
        xb(chagB) = xb(chagB)/A(chagB,ind);
        A(chagB,:) = A(chagB,:)/A(chagB,ind);
        for i = 1 : n
            if i~=chagB
                xb(i) = xb(i) - A(i,ind)*xb(chagB);
                A(i,:) = A(i,:) - A(i,ind)*A(chagB,:);
                
            end
        end
        baseVector(chagB) = ind;
end
M = M + 1;
if (M == 1000000)
    disp('�Ҳ������Ž⣡');
    mx = NaN;
    minf = NaN;
    return;
end
        
end