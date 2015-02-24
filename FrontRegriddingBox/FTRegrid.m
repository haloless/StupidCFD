
function [] = FTRegrid()

FTRegridGlobals;

ipassmax = 6;

for ipass = 1:ipassmax
    
    neladd = FTRegridAddElem();
    
    neldel = FTRegridDelElem();
    
    fprintf('IPASS=%d,NADD=%d,NDEL=%d\n', ipass, neladd,neldel);
    
    if (neladd==0 && neldel==0); break; end
end

return
end

