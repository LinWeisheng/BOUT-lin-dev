from boututils.datafile import DataFile
from boutdata.collect import collect
from boututils.showdata import showdata

T = collect("T")
showdata(T[:,0,:,0])
