
function [ newid num head tail free next prev ] = LinkListAdd(...
num,head,tail,free,next,prev, maxnum)

if num+1 > maxnum
    fprintf('Link List: MAXNUM=%d\n',maxnum);
    error('Link List overflow');
end

% sieze new ID
newid = free;
free = next(free);

next(newid) = next(tail);
next(tail) = newid;
prev(newid) = tail;
prev(head) = newid;

tail = newid;

num = num + 1;

return
end


