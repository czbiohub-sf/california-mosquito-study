import pysam
import os.path

culex_names = []
with open ("../mito_cms1_Culex.txt") as f:
  for line in f:
    culex_names.append(line)

culex_names = [x.strip() for x in culex_names] 

NM_scores = []

for x in culex_names:
  if (os.path.isfile(x+".sam")):
    samfile = pysam.AlignmentFile(x+".sam", "rb")
    current_NM = []
    unmapped = 0
    count = 0
    for read in samfile.fetch():
      if (count >= 1000):
        break
      if (read.has_tag("NM")):
        current_NM.append(read.get_tag("NM"))
        count = count + 1
      else:
        unmapped = unmapped + 1
    prop = current_NM.count(0)/len(current_NM)
    print (x+' '+str(prop)+'\n')
    NM_scores.append(current_NM)
  else:
    NM_scores.append(0)


