function getAmpPhs(dataValue::Array)

rdata = real(dataValue)
idata = imag(dataValue)

amp = sqrt.(rdata.^2 + idata.^2)
pha = atan.(rdata, idata)

DSP.unwrap!(pha)
pha = pha*180/pi

return amp, pha
end
