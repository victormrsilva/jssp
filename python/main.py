from instance import Instance
import sys
import xlsxwriter
import os, fnmatch
import subprocess

if __name__ == "__main__":
  instance = Instance(sys.argv[1], int(sys.argv[2]))
  instance.print()

  # colors to use
  colors = ['red','blue','green','brown','gray','pink','yellow','black','cyan','gray','lime','magenta','navy','orange','purple','silver','white']

  # Create a workbook and add a worksheet.
  filename = "{}.xlsx".format(sys.argv[1])
  workbook = xlsxwriter.Workbook(filename)
  worksheet = workbook.add_worksheet()
  worksheet.set_column(1,int(sys.argv[2])+1,3)


  listOfFiles = os.listdir('.')  
  pattern = "{}*.sol".format(sys.argv[1])
  row = 0
  col = 0
  for entry in listOfFiles:  
    if fnmatch.fnmatch(entry, pattern):
      print(entry)
      worksheet.write(row,col,entry)
      row = row+1
      for i in range(int(sys.argv[2])+1):
        col = col+1
        worksheet.write(row,col,i)
      col = 0
      for i in range(instance.getQtdMachines()):
        row = row+1
        worksheet.write(row,col,"m{}".format(i))
      row = row - instance.getQtdMachines()+1
      with open(entry,"r") as solution_file:
        for line in solution_file:
          tokens = line.split()
          if (tokens[1].find("x") >= 0 and float(tokens[2]) > 0):
            values = tokens[1]
            values = values.replace("x","")
            values = values.replace("(","")
            values = values.replace(")","")
            values = values.split(",")
            print(values)
            task,machine,time = int(values[0])-1,int(values[1])-1,int (values[2])
            worksheet.write(row+machine,col,"m{}".format(machine))
            for i in range(instance.getTime(task,machine)):
              cell_format = workbook.add_format()
              cell_format.set_pattern(1)
              cell_format.set_bg_color(colors[task])
              worksheet.write(row+machine,time+1+i,"",cell_format)

            
      row = row+instance.getQtdMachines()+2

  workbook.close()
