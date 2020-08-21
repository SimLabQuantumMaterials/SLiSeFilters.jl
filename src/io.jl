"Import filter from file"
function readFilter(filename, skipstart = 2)
	saved = readdlm(filename,',', skipstart = skipstart)
	z1  = saved[1,:] + im*saved[2,:]
	a1  = saved[3,:] + im*saved[4,:]

	return z1,a1
end

"Import weight function from file"
function readWS(filename)
	saved = readdlm(filename,',')
	intv  = saved[1,:]
	wght  = saved[2,:]

	return intv,wght
end
