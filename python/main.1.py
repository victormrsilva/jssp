from instance import Instance
from flow import Flow
import sys
import xlsxwriter
import os, fnmatch
import subprocess

if __name__ == "__main__":
  instance = Instance(sys.argv[1], int(sys.argv[2]))
  
  flow = Flow(instance)
  flow.constructProblem()

  listOfFiles = os.listdir('.')  
  pattern = "{}*.sol".format(sys.argv[1])
  # for entry in listOfFiles:  
  #   if fnmatch.fnmatch(entry, pattern):
  #     print(entry)
  #     cat = subprocess.Popen(["cat",entry],stdout=subprocess.PIPE)
  #     grep = subprocess.Popen(["grep", "x"],stdin=cat.stdout,stdout=subprocess.PIPE)
  #     grep2 = subprocess.Popen(["grep", " 1"],stdin=cat.stdout,stdout=subprocess.PIPE)
  #     endOfPipe = grep.stdout
  #     for line in endOfPipe:
  #       print(line) 
