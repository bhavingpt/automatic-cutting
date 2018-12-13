import os, sys

os.chdir(sys.argv[1])
os.system("rm -rf *.asc")
for i in range(1, 6):
    for j in range(1, int(sys.argv[2]) + 1):
        os.system("touch cut" + str(i) + "_" + str(j) + ".asc")
        os.system("touch wall" + str(i) + "_" + str(j) + ".asc")
