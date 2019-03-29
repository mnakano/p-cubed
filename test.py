str = '------HHHHHHH---HHHHHHHHHHH--------HHHHHH-------HHHH--'

list = [[str[0]]]
num = 0
print(str)
for i in range(1,len(str)):
	if(str[i-1] == str[i]):
		list[num].append(str[i])
	else:
		list.append([str[i]])
		num += 1

for l in list:
	print(l)