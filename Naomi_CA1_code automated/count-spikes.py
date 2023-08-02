threshold = 0
timediff = 40

numspikes = 0
lastspike = -100000
with open("input.txt", "r") as df:
	lines = df.readlines()
	for line in lines:
		time, value = map(float, line.rstrip().split('\t'))
		if value > threshold and time - lastspike > timediff:
			lastspike = time
			numspikes += 1
print(f'There were {numspikes} spikes, with threshold {threshold}')