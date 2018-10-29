import sys

for line in sys.stdin:
	if line.startswith(">"):
		print(line.rstrip("\n")[:50])
	else:
		print(line, end="")
