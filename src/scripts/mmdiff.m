% read 2 mm files, take the difference, and write as diff.mm
function  [out] = mmdiff(file1,file2)
	[a] = mmread(file1)
	[b] = mmread(file2)
	c = a - b
	mmwrite("diff.mm",c)
end

