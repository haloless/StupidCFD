
function [ num head tail free next prev ] = LinkListDel(oldid, ...
num,head,tail,free,next,prev, maxnum)

% adjust the used list
if oldid == head
    head = next(head);
end
if oldid == tail
    tail = prev(oldid);
end
if prev(oldid)==0 || next(oldid)==0
    fprintf('ID=%d,PREV=%d,NEXT=%d\n',oldid,prev(oldid),next(oldid));
end
next(prev(oldid)) = next(oldid);
prev(next(oldid)) = prev(oldid);

% adjust the free list
next(oldid) = free;
prev(oldid) = prev(free);
prev(free) = oldid;
free = oldid;

num = num - 1;

return
end

