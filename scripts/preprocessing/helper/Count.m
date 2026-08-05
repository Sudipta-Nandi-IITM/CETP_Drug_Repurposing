liste = dir("*.txt");
files = {liste.name};
VALUE = []
TOTAL = []
for k = 1:(length(files))
    table = readtable(string(files(k)))
    a = table{:,3}
    num = groupcounts(a>7)
    value = num(1) - 1
    VALUE = [VALUE,value]
    total = length(a) -1 
    TOTAL = [TOTAL,total]
end

