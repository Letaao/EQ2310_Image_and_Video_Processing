function dctblock=dct16(block)
dctblock=blkproc(block,[8 8],@dct2);
end