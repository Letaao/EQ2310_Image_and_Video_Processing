function quant = midtread(x,stepsize)
% uniform mid-tread quantizer function
quant = round(x/stepsize)*stepsize;
end