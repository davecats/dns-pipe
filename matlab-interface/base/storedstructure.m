function mmap=storedstructure(fname,field,writable)

FSIZE=0;
for i=1:numel(field)/3
    if numel(field)/3 == 1; c=field{1}; else c=field{i,1}; end
    switch c
        case 'double'
            bs=8;
        case 'complex'
            bs=16;
        case 'uint8'
            bs=1;
    end
    if numel(field)/3 == 1; s=field{2}; else s=field{i,2}; end
    FSIZE=FSIZE+bs*prod(s);
end
command = sprintf('%s%u %s', 'truncate --size=', FSIZE, fname);
unix(command);

mmap=memmapfile(fname, ...
                 'Format', field, ...
                 'Writable',writable);
