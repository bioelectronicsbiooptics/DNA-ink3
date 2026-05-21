% Textscan_240108_v1
%% Paired End Sequence Read for cultural heritage files

filename = 'Jikji6' 

tStart = tic;
tic 

for i = 1 : 2
    fid = fopen("6_S2_L001_R"+i+"_001"+".fastq");
    warning('off','all')
    data = textscan(fid,'%s','Delimiter','\n','CommentStyle',{'+','@'});
    warning('on','all')
    if i == 1
        seq_1 = data{1}(2:2:end);
    else
        seq_2 = data{1}(2:2:end);
    end
    fclose(fid);
end
seq = [seq_1;seq_2];
T(1)= toc; 
"NGS Read: " + T(1)
%% Reverse Complemantary Sequence
tic
ForwardSection = cell(length(seq),1);

chunkSize = 1000; 
numChunks = ceil(length(seq) / chunkSize);
ForwardSection = cell(size(seq)); 

for i = 1:numChunks    
    startIndex = (i - 1) * chunkSize + 1;
    endIndex = min(i * chunkSize, length(seq));
    for j = startIndex:endIndex
        ForwardSection{j} = seq{j}(1:21);
    end
end

primersData = struct('filename', {'Jikji1', 'Jikji2', 'Jikji3', 'Jikji4', 'Jikji5', 'Jikji6', 'Jikji7', 'tripitaka', 'liberty', 'seokga', 'pisa', 'venus', 'goldengate'}, ...
                     'F_primers', {'TTCGTTCGTCGTTGATTGGT', 'AAATCCTTTGTGCCTGCCAT', 'AATCATGGCCTTCAAACCGT', 'CTGTCCATAGCCTTGTTCGT', 'TGTATTTCCTTCGGTGCTCC', 'AGCCTTGTGTCCATCAATCC', 'GTCCAGGCAAAGATCCAGTT', 'TAGCCTCCAGAATGAAACGG', 'AAGGCAAGTTGTTACCAGCA', 'ATTCGCGTCGCCTAATTTGT', 'AATGGACGTTCCGCAATCAT', 'ATCCTGCAAACGCATTTCCT', 'TCCACCGTTCCTTGATTTCG'}, ...
                     'R_primers', {'AAACGGAGCCATGAGTTTGT', 'AAACTCAAGGCCGACCAATT', 'AACGCTCCGAAAGTCTTGTT', 'GCGGAAACGTAGTGAAGGTA', 'TTTCGACAACGGTCTGGTTT', 'TGCGCTATGGTTTGGCTAAT', 'ACCACCGTTAGGCTAAAGTG', 'TTCAAGCCAAACCGTGTGTA', 'TGCGACCGTAATCAAACCAA', 'AAACTGGAGGCGGCAAATTA', 'AGAGCCGTGGCAATGTAAAT', 'ATGCCTTTCCGAAGTTTCCA', 'AATCCGTTTGCCTGCCTTTA'}, ...
                     'OligoLen', {981, 1017, 1060, 1146, 895, 1018, 1339, 92626, 45381, 2586, 43138, 63952, 35464});

% Accessing the data
for i = 1:length(primersData)
    if strcmp(filename, primersData(i).filename)
        F_primers = primersData(i).F_primers;
        R_primers = primersData(i).R_primers;
        OligoLen = primersData(i).OligoLen;
        break;
    end
end

RC_R_primer = seqrcomplement(R_primers); 
S = cell2mat(ForwardSection);
F20 = S(:,1:20);
F21 = S(:,2:21);
ReversedSeqLog = zeros(length(seq),1);

parfor a = 1 : length(seq)
    if sum(F20(a,:) == RC_R_primer) >= 18 || sum(F21(a,:) == RC_R_primer) >= 18
        ReversedSeqLog(a) = 1;
    end
end

evertSeq = seq(logical(ReversedSeqLog));
tevertSeq = cell(1,length(evertSeq));

% option 1 (T-21)
parfor a = 1:length(evertSeq)
    seqLength = length(evertSeq{a});
    if strcmp(evertSeq{a}(1:20), RC_R_primer) && seqLength >= 150
        Reverse = evertSeq{a}(1:150);
        tevertSeq{a} = seqrcomplement(Reverse);
    end
    if strcmp(evertSeq{a}(2:21), RC_R_primer) && seqLength >= 151
        Reverse = evertSeq{a}(2:151);
        tevertSeq{a} = seqrcomplement(Reverse);
    end
end

seq(logical(ReversedSeqLog)) = tevertSeq;
T(2)= toc; 
"Reverse Complemantary Sequence: "+ T(2)
%% Primer Sorting
tic
PrimerSortingSeq = cell(length(seq),1);
for a = 1 : length(seq)
    W = seq{a};
    if length(W) >= 150
        if sum(W(1:20) == F_primers) >=18
            if sum(W(131:150) == R_primers) >= 18
                PrimerSortingSeq{a} = W(1:150);
            end
        end
    end
    if length(W) >= 151
        if sum(W(2:21) == F_primers) >=18
            if sum(W(132:151) == R_primers) >= 18

                PrimerSortingSeq{a} = W(2:151);
            end
        end

    end
end
PrimerSortingSeq = PrimerSortingSeq(~cellfun('isempty',PrimerSortingSeq));
T(3)=toc;
"Primer Sorting: "+T(3)
% tMul = sum(T)
%% Key finding
tic

DNAmatrix = cell2mat(PrimerSortingSeq);

payload_part = DNAmatrix(:,21:98);
KeyMIndex = zeros(length(payload_part),1);

parfor a = 1 : size(payload_part,1)
    if numel(strfind(payload_part(a,:),'AAAAA')) >= 1
        if numel(strfind(payload_part(a,:),'AAAA')) - ...
                numel(strfind(payload_part(a,:),'AAAAA')) + ...
                numel(strfind(payload_part(a,:),'GAAAT')) >= 3
            KeyMIndex(a) = 1;
        end
    end
end

pre_Key_sorting = PrimerSortingSeq(find(KeyMIndex),:);
mpre_Key_sorting = cell2mat(pre_Key_sorting);
Key_sorting_cell = cell(size(mpre_Key_sorting,1),1);
for i = 1 : size(mpre_Key_sorting,1)
    Key_sorting_cell{i} = mpre_Key_sorting(i,21:98);
end

moduleDNA = [repmat('TGCA',1,27) 'TG'];
KeyD = cell(length(Key_sorting_cell),1);
parfor a = 1 : length(Key_sorting_cell)
    temp = Key_sorting_cell{a};
    for b = 73 : -1 : 58
        keyi = length(temp(b+1:end));
        if strcmp(temp(b-4:b),repmat('A',1,5))
            if strcmp(temp(b-keyi-10:b-keyi-1),moduleDNA(b-keyi-10:b-keyi-1))
                KeyD{a} = temp;
            end
        end
    end
end
KeyD = KeyD(~cellfun('isempty',KeyD));

KeymatDNA = mode(cell2mat(KeyD),1);
k = 156;
n = k+8;
m = length(de2bi(n));
d = {'00','01','10','11'};
for b = 1 : size(KeymatDNA,2)-3
    if strcmp(KeymatDNA(b:b+3),'AAAA')
        KeyMaptemp = KeymatDNA(1:b-1);
        if (strfind(KeyMaptemp,'N'))>0
            KeyMaptemp = repmat('A',1,length(KeyMaptemp));
            continue
        end
        [~,x] = ismember(KeyMaptemp,'ATGC');
        chr_out = cell2mat(d(x));
        tempKeyMap = chr_out-'0';
        if ~isempty(tempKeyMap)
            KeyMap = bi2de(tempKeyMap);
        end
        break
    end
end

strkey = num2str(KeyMap);
if size(strkey,2) == 5
    ncolposz = str2double(strkey(1));
    nbitposz = str2double(strkey(2));
    ncolcosz = str2double(strkey(3:4));
    nbitcosz = str2double(strkey(5));
    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz;
elseif size(strkey,2) == 6
    ncolposz = str2double(strkey(1:2));
    nbitposz = str2double(strkey(3));
    ncolcosz = str2double(strkey(4:5));
    nbitcosz = str2double(strkey(6));
    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz;
end

[~,x] = ismember(KeymatDNA,'ATGC');
chr_out = cell2mat(d(x));
Keybin = chr_out-'0';

if KeyMap > 120000
    binkloc = Keybin(length(de2bi(KeyMap))+9:length(de2bi(KeyMap))+9+keysz);
else
    binkloc = Keybin(length(de2bi(KeyMap))+10:length(de2bi(KeyMap))+10+keysz);
end
ncolpo = bi2de(binkloc(1:ncolposz));
nbitpo = bi2de(binkloc(ncolposz+1:ncolposz+nbitposz));
ncolco = bi2de(binkloc(ncolposz+nbitposz+1:ncolposz+nbitposz+ncolcosz));
nbitco = bi2de(binkloc(ncolposz+nbitposz+ncolcosz+1:end));

nfile = 3*ncolpo*nbitpo+3*ncolco*nbitco;
tidx = 139;
tidxlen = ceil((nfile)/tidx*1.5);

while tidx+length(de2bi(ceil((nfile)/tidx*1.5))) < k
    tidx = tidx + 1;
end
idxlen = k-tidx;
keyidx = bi2de(Keybin(end-idxlen+1:end));
OligoLen =keyidx;
T(4)= toc; 
"Key finding: " + T(4)
% tMul = sum(T)
%% "Consensus Sequencing & binarization"
tic
IndexDNA = DNAmatrix(:,98-ceil(idxlen/2)+1:98);
Indexbin = zeros(length(IndexDNA),size(IndexDNA,2)*2);

for r = 1:size(IndexDNA, 1)
    if any(IndexDNA(r,:) == 'N')
         IndexDNA(r,:) = repmat('A', 1, size(IndexDNA, 2));
    end
end
parfor r = 1:size(IndexDNA, 1)
         [~, x] = ismember(IndexDNA(r,:), 'ATGC');
        chr_out = cell2mat(d(x));
        Indexbin(r,:) = chr_out - '0';
end

if rem(length(de2bi(keyidx)),2) ~= 0
    if rem(size(Indexbin,2),2) == 0
        Indexbin = Indexbin(:,2:end);
    end
end

Indexdec = bi2de(Indexbin);
indicesCell = cell(keyidx, 1);
parfor a = 1:keyidx 
    indicesCell{a} = find(Indexdec == a);
end

for a = 1:keyidx
    if ~isempty(indicesCell{a})
        MIndex{a} = DNAmatrix(indicesCell{a}, :);
    else
        MIndex{a} = []; 
    end
end

if strcmp(filename,'pisa') || strcmp(filename,'seokga') || strcmp(filename,'liberty')
    key_map = MIndex{end};
    KeyIndex = zeros(size(key_map,1),1);
    for a = 1 : size(key_map,1)
        for b = 21 : 98
            if sum(key_map(a,b:b+3) == 'AAAA') == 4
                KeyIndex(a) = 1;
                break
            end
        end
        for c = b+3 : 98
            if sum(key_map(a,c:c+3) == 'AAAA') == 4
                KeyIndex(a) = KeyIndex(a) +1;
                break
            end
        end
        for d = c+3 : 98
            if sum(key_map(a,d:d+3) == 'AAAA') == 4
                KeyIndex(a) = KeyIndex(a) +1;
                break
            end
        end
    end
    MapLocation = key_map(~(KeyIndex==3),:);
    KeyIndex = key_map((KeyIndex==3),:);
    MIndex(end) = {KeyIndex};
end  

dataMode = cell(length(MIndex),1);
for a = 1 : length(MIndex)
    if ~isempty(MIndex{a})
        dataMode{a} = mode(MIndex{a},1);
    end
end
dataMode(cellfun('isempty',dataMode)) = {repmat('A',1,150)};
DNAdata = cell2mat(dataMode);
T(5)=toc;
"Consensus Sequencing: " +T(5)
% tMul = sum(T)

MatData = DNAdata(:,21:130);
SeqtoBin = zeros(length(MatData),size(MatData,2)*2);
d = {'00','01','10','11'};
for r = 1 : size(MatData,1)
    Nnum = strfind(MatData(r,:),'N');
    if  numel(Nnum) >  0
        for i = 1 : numel(Nnum)
            MatData(r,Nnum(i)) = 'A';
        end
    end
    [~,x] = ismember(MatData(r,:),'ATGC');
    chr_out = cell2mat(d(x));
    SeqtoBin(r,:) = chr_out-'0';
end
binKey = SeqtoBin(end,:);
keyline = binKey(1:k);

T(6)=toc;
"Binarization: " +T(6)
% tMul = sum(T)

Depth = zeros(length(MIndex),1);
for i = 1 : length(MIndex)
    Depth(i) = size(MIndex{i},1);
end
"Depth mean: "+mean(Depth)
%%
tic

bin_payload = SeqtoBin(:, 1:k);
bin_redun = SeqtoBin(:, k+1:end);
dec_redun = zeros(size(SeqtoBin, 1), (size(SeqtoBin, 2) - k) / m);


for z = m : m : size(SeqtoBin, 2) - k
    if z == m 
        dec_redun(:, z/m) = bi2de(bin_redun(:, 1:m));
    else
        dec_redun(:, z/m) = bi2de(bin_redun(:, z-m+1:z));
    end
end


rscode = gf([bin_payload dec_redun], m);


chunkSize = 1000;
numRows = size(rscode, 1);
rsECC = zeros(numRows, size(rscode, 2) - (n - k));

errCount = zeros(numRows, 1);


for startIdx = 1:chunkSize:numRows
    endIdx = min(startIdx + chunkSize - 1, numRows);
    [rxcode, cnumerr] = rsdec(rscode(startIdx:endIdx, :), n, k);
    rxcode = rxcode.x;
    rsECC(startIdx:endIdx, :) = double(rxcode);
    errCount(startIdx:endIdx) = cnumerr;
end

GRS_err_indices = find(errCount<0);
rsECC(find(errCount), :) = zeros(length(find(errCount)), size(rsECC, 2));

idxECC = bi2de(rsECC(:, k - idxlen + 1:end));
payload = rsECC(:, 1:k - idxlen);
Temp_mat = payload(1:end-1,:);
T(7)=toc;
"RS Decoding: " +T(7)
% tMul = sum(T)
% tic
% % tMul = sum(T)
%% Null finding
tic
if rem(ncolpo*3*nbitpo+ncolco*3*nbitco,2) == 1
    Nullkey_1 = ncolpo*3*nbitpo+ncolco*3*nbitco+1;
else
    Nullkey_1 = ncolpo*3*nbitpo+ncolco*3*nbitco;
end
s = size(Temp_mat,2);

Nullkey_2 = Nullkey_1/2;
Null_count = 0;
while rem(Nullkey_1,s) ~= 0
    Nullkey_1 = Nullkey_1 + 1;
    Null_count = Null_count + 1;
end

Null_xor_count = 0;
while rem(Nullkey_2,s) ~= 0
    Nullkey_2 = Nullkey_2 + 1;
    Null_xor_count = Null_xor_count + 1;
end

if rem(s,2) ~= 0 && rem(Null_count,2) ~= 0 && rem(Null_xor_count,2) ~= 0
    data_Nrng = (s/2-Null_count/2+1:floor(s/2))+20;
    xor_Nrng = (s/2-Null_xor_count/2+1:floor(s/2))+20;
elseif (rem(s,2) ~= 0 && rem(Null_count,2) == 0 && rem(Null_xor_count,2) == 0)...
    || (rem(s,2) == 0 && rem(Null_count,2) == 0 && rem(Null_xor_count,2) == 0)
    data_Nrng = (floor(s/2)-Null_count/2+1:floor(s/2))+20;
    xor_Nrng = (floor(s/2)-Null_xor_count/2+1:floor(s/2))+20;
else
    data_Nrng = (floor(s/2)-Null_count/2+1:floor(s/2))+20;
    xor_Nrng = (floor(s/2)-ceil(Null_xor_count/2)+1:floor(s/2))+20;
end
data_polyA = repmat('A',1,length(data_Nrng));
xor_polyA = repmat('A',1,length(xor_Nrng));

data_Noligo = zeros(length(seq),1);
xor_Noligo = zeros(length(seq),1);
for a = 1 : length(seq)
    if length(seq{a}) >= s/2+20
        if sum(seq{a}(data_Nrng) == data_polyA) == length(data_polyA)
            if sum(seq{a}(1:20) == F_primers) >= 15
                data_Noligo(a) = 1;
            end
        end

        if sum(seq{a}(xor_Nrng) == xor_polyA) == length(xor_polyA)
            if sum(seq{a}(1:20) == F_primers) >= 15
                xor_Noligo(a) = 1;
            end
        end
    end
end

data_Nsort = seq(logical(data_Noligo))';
xor_Nsort = seq(logical(xor_Noligo))';
if gt(numel(xor_Nrng),numel(data_Nrng))
    for a = 1 : length(data_Nsort)
        if sum(data_Nsort{a}(xor_Nrng) == xor_polyA) == length(xor_polyA)
            data_Nsort{a} = {};
        end
    end
    cellidx = ~cellfun('isempty',data_Nsort);
    data_Nsort = data_Nsort(cellidx);
else
    for a = 1 : length(xor_Nsort)
        if sum(xor_Nsort{a}(data_Nrng) == data_polyA) == length(data_polyA)
            xor_Nsort{a} = {};
        end
    end
    cellidx = ~cellfun('isempty',xor_Nsort);
    xor_Nsort = xor_Nsort(cellidx);
end

dataNcons = cell(length(sum(data_Noligo)),1);
for a = 1 : length(data_Nsort)
    dataNcons{a} = data_Nsort{a}(1:min(data_Nrng)-1);
end

xorNcons = cell(length(sum(xor_Noligo)),1);
for a = 1 : length(xor_Nsort)
    xorNcons{a} = xor_Nsort{a}(1:min(xor_Nrng)-1);
end

data_Ndecode = mode(cell2mat(dataNcons'),1);
xor_Ndecode = mode(cell2mat(xorNcons'),1);
T(8) = toc;
"Null finding " + T(8)


fid = fopen(filename+".txt");
[scan,count] = fscanf(fid,'%s');
orgseq = reshape(scan,[length(scan)/count count])';

if sum(data_Ndecode(21:end) == orgseq(floor((numel(Temp_mat)*2/3)/s),21:length(data_Ndecode))) == length(data_Ndecode)-20
    "Data part Null Decoding Complete"
end

if sum(xor_Ndecode(21:end) == orgseq(end-1,21:length(xor_Ndecode))) == length(xor_Ndecode)-20
    "XOR part Null Decoding Complete"
end
%% XOR
tic
[~,x] = ismember(data_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(x));
bindataN = chr_out-'0';

[~,y] = ismember(xor_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(y));
binxorN = chr_out-'0';

Temp_mat(floor((numel(Temp_mat)*2/3)/s),1:length(bindataN)) = bindataN; 
Temp_mat(end,1:length(binxorN)) = binxorN;

tri_part = Temp_mat(1:floor(numel(Temp_mat)*2/3/s),:);
xor_part = Temp_mat(floor(numel(Temp_mat)*2/3/s)+1:end,:);

row_tri_part = reshape(tri_part',[numel(tri_part) 1])';
row_tri_part(end-Null_count+1:end) = [];

if rem(ncolpo*3*nbitpo+ncolco*3*nbitco,2) ~= 0
    row_tri_part(end) = [];
end

XOR_ECC_C = reshape(xor_part',[numel(xor_part) 1])';
XOR_ECC_C(end-Null_xor_count+1:end) = [];

XOR_ECC_A = row_tri_part(1:ceil(length(row_tri_part)/2));
XOR_ECC_B = row_tri_part(length(XOR_ECC_A)+1:end);
if rem(length(row_tri_part),2) ~= 0
    XOR_ECC_B(end+1) = 0;
end
XOR_ECC_matrix = [xor(XOR_ECC_B,XOR_ECC_C) xor(XOR_ECC_A,XOR_ECC_C)];

ECC_example = [XOR_ECC_matrix zeros(1,Null_count)];
ECC_ex_matrix = reshape(ECC_example,[s numel(ECC_example)/s])';
missing_data = find(sum(tri_part')==0);
tri_part(missing_data,:) = ECC_ex_matrix(missing_data,:);
ECC_tri_part = tri_part;
% T(9)= toc;
% "XOR " + T(9)
%% G-RS tri
tic
ECC_payload = reshape(ECC_tri_part',[numel(ECC_tri_part) 1])';
ECC_payload(end-Null_count+1:end) = [];

if rem(ncolpo*3*nbitpo+ncolco*3*nbitco,2) ~=0
    ECC_payload(end) = [];
end

po2row = ECC_payload(1:ncolpo*nbitpo*3);
co2row = ECC_payload(ncolpo*nbitpo*3+1:end);

po_bin = reshape(po2row',[nbitpo ncolpo*3'])';
po_dec = bi2de(double(po_bin))/(10^16);

tri_points = reshape(po_dec,[3 ncolpo])';
if strcmp(filename,'seokga')||strcmp(filename,'liberty')||strcmp(filename,'pisa')...
         ||strcmp(filename,'venus')||strcmp(filename,'goldengate')||strcmp(filename,'tripitaka')    
        tri_points = reshape(po_dec,[ncolpo 3]);
end

co_bin = reshape(co2row',[nbitco ncolco*3])';
co_dec = bi2de(double(co_bin));
tri_connects = reshape(co_dec,[3 ncolco])';

if numel(find(~tri_points')) > 0
    XOR_rep = find(~tri_points');
    XOR_prng = [XOR_rep*nbitpo-nbitpo+1 XOR_rep*nbitpo];
    eccxorpo = zeros(size(XOR_prng,1),nbitpo);
    for a = 1 : size(XOR_prng,1)
        eccxorpo(a,:) = XOR_ECC_matrix(XOR_prng(a,1) : XOR_prng(a,2));
    end
    respo = bi2de(eccxorpo)/10^16;
    transposePO = tri_points';
    transposePO(XOR_rep) = respo;
    tri_points = transposePO';
end
% T(10)= toc;
% "Global RS XOR " + T(10)

if numel(find(~tri_connects')) > 0 || (numel(find((tri_connects'>size(tri_points,1)))))
    XOR_req = [find(~tri_connects'); find(tri_connects'>size(tri_points,1))];
    XOR_rng = [XOR_req*nbitco-nbitco+1 XOR_req*nbitco];
    coXOR1 = XOR_ECC_B(ncolpo*nbitpo*3-length(XOR_ECC_A)+1:end);
    coXOR2 = XOR_ECC_C(ncolpo*nbitpo*3-length(XOR_ECC_A)+1:end);
    coXOR3 = XOR_ECC_A(ncolpo*nbitpo*3-length(XOR_ECC_A)+1:end);
    if rem(length(row_tri_part),2) ~= 0
        coXOR1(end) = [];
        coXOR2(end) = [];
        coXOR3(end) = [];
    end
    XORrest = xor(coXOR1,coXOR2);
    XORrest2 = xor(coXOR2,coXOR3);
    XR = [];
    for xr = 1 : size(XOR_rng,1)
        XC = XORrest2(XOR_rng(xr,1):XOR_rng(xr,2));
        XR = [XR; XC];
    end
    transposeCO = tri_connects';
    transposeCO(XOR_req) = bi2de(XR);
    tri_connects = transposeCO';

    nullsz = [find(tri_connects > ncolpo); find(tri_connects == 0)];
    for a = 1 : length(nullsz)
        while nullsz(a) > ncolco
            nullsz(a) = nullsz(a) -ncolco;
        end
    end
    size(tri_connects(unique(nullsz),:),1)
    tri_connects(unique(nullsz),:) = [];
end

trid = triangulation(tri_connects,tri_points);
trimesh(trid)
T(9)= toc;
"tri " + T(9)
tMul = sum(T)
tMul/60