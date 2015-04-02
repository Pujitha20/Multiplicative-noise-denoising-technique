function sq=plotwavelet2(C,S,level,wavelet,rv,mode)

A = cell(1,level); H = A; V = A; D = A;

for k = 1:level
    A{k} = appcoef2(C,S,wavelet,k); % approx
    [H{k} V{k} D{k}] = detcoef2('a',C,S,k); % details  
    
    A{k} = wcodemat(A{k},rv);
    H{k} = wcodemat(H{k},rv);
    V{k} = wcodemat(V{k},rv);
    D{k} = wcodemat(D{k},rv);
end

if strcmp(mode,'tree')
    
    aff = 0;
    
    for k = 1:level
        subplot(level,4,aff+1); image(A{k});
        title(['Approximation A',num2str(k)]);
        subplot(level,4,aff+2); image(H{k});
        title(['Horizontal Detail ',num2str(k)]);
        subplot(level,4,aff+3); image(V{k});
        title(['Vertical Detail ',num2str(k)]);
        subplot(level,4,aff+4); image(D{k});
        title(['Diagonal Detail ',num2str(k)]);
        aff = aff + 4;
    end
    
elseif strcmp(mode,'square')
    
    dec = cell(1,level);
    dec{level} = [A{level} H{level} ; V{level} D{level}];
    
    for k = level-1:-1:1
        dec{k} = [imresize(dec{k+1},size(H{k})) H{k} ; V{k} D{k}];
    end
    
    imshow(uint8(dec{1}));
    sq=uint8(dec{1});
    
end

end