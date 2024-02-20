H = [1 0 0 0 1 1 1 0 0 0 1 1 1 0 1;
    0 1 0 0 1 0 0 1 1 0 1 1 0 1 1;
    0 0 1 0 0 1 0 1 0 1 1 0 1 1 1;
    0 0 0 1 0 0 1 0 1 1 0 1 1 1 1];

x1 = [1 1 1 1 1 1 1 0 0 1 0 0 0 0 0];
x2 = [1 0 0 1 1 0 0 1 0 1 0 0 0 0 0];
x3 = [1 0 0 0 1 0 0 0 0 0 0 0 0 0 1];

r1 = checkCodewords(x1',H);
r2 = checkCodewords(x2',H);
r3 = checkCodewords(x3',H);

[E,S] = buildTable(H);

%% 6 b). 
x6 = [1 1 1 1 1 0 0 0 0 1 0 0 0 0 0];

for i = 1:15
    r6 = x6;
    if x6(i) ~= 1
        r6(i) = 1;
    else
        r6(i) = 0;
    end
    [e,x] = channelDecode(r6',H);
    disp(r6);
    disp(e');
    disp(x');
    disp('\n');
end

%% 6 c).
x6 = [1 1 1 1 1 0 0 0 0 1 0 0 0 0 0];

for i = 1:14
    r6 = x6;
    if x6(i) ~= 1
        r6(i) = 1;
    else
        r6(i) = 0;
    end

    if x6(i+1) ~= 1
        r6(i+1) = 1;
    else
        r6(i+1) = 0;
    end

    [e,x] = channelDecode(r6',H);
    disp(r6);
    disp(e');
    disp(x');
    disp('\n');
end
%%
% checkCodewords
function result = checkCodewords(x,H)
    temp = H * x;
    ones = 0;
    for i = 1:length(temp)
        if mod(temp(i),2) ~= 0
            ones = ones+1;
        end
    end

    if ones ~= 0
        result = 0;
        disp('The vector is not in the codebook.');
    else
        result = 1;
        disp('The vector is in the codebook');
    end
end

% buildTable
function [E,S] = buildTable(H)
    colm= size(H,2);
    E = zeros(colm,colm+1);

    for i = 2:colm+1
        E(i-1,i) = 1;
    end

    S = H * E;
end

% channelDecode
function [e,x] = channelDecode(r,H)
    [E,S] = buildTable(H); 
    temp = H*r; % normal matrix vector product
    s = zeros(size(temp)); 

   for i = 1: length(temp) % convert to binary
       if mod(temp(i),2) ~= 0
            s(i) = 1;
        end
   end
  loc = 1;

  while loc >= 1
    if ~isequal(S(:,loc),s)
        loc = loc + 1;
    else
        e = E(:,loc);
        loc = 0;
    end
  end
  
  x = max(r - e,0); % convert all -1 in x to 0
end