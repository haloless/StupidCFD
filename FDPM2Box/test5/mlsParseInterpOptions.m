
kvpairs = varargin;
for i = 1:2:length(kvpairs)
    key = kvpairs{i};
    val = kvpairs{i+1};
    
    switch key
        case 'mls_order'
            mls_order = val;
        case 'mls_volume'
            mls_volume = val;
        otherwise
            error('Unknown %s',kvpairs{i});
    end
end


