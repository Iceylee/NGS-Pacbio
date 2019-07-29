#!/usr/bin/python
"""
txt2xls.py [text-file-folder]
text必须都是UTF-8编码（注意COG和summary）
否则会报错：UnicodeDecodeError: 'utf-8' codec can't decode byte 0x91 in position 803: invalid start byte
"""

import sys
mypath = sys.argv[1]

from os import listdir
from os.path import isfile, join
textfiles = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) and '.txt' in  f]

#判断是否为数字类型的字符（"1.11" "1" "3e-5"）
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

import xlwt
import xlrd

font = xlwt.Font()
font.name = 'Times New Roman'
#font.bold = True

#整数
styleInt = xlwt.XFStyle()
styleInt.num_format_str = '#,###'
styleInt.font = font

#文字
styleALL = xlwt.XFStyle()
styleALL.font = font

#小数
styleFloat= xlwt.XFStyle()
styleFloat.num_format_str = '#,###0.00'
styleFloat.font = font

#科学记数
styleEE= xlwt.XFStyle()
styleEE.num_format_str = '0.00E+0'
styleEE.font = font

#header
styleHeader  =xlwt.easyxf('font: color-index green, name Times New Roman, bold on'); 

for textfile in textfiles:
#     f = open(textfile, 'r+')
    row_list = []
#     for row in f:
#         row_list.append(row.split('\t'))
    
    with open(textfile, 'r') as f:   
        for row in f.readlines():
            row_list.append(row.split('\t'))
    
        column_list = zip(*row_list)
        # for column_list in f:
        #     column_list.append(column.split('|'))
        workbook = xlwt.Workbook()
        worksheet = workbook.add_sheet('Sheet1')
        i = 0
        for column in column_list:
            for item in range(len(column)):
                value = column[item].strip()
                if item == 0:
                    worksheet.write(item, i, value, style=styleHeader)
                else:
                    if is_number(value): #数字类型的字符
                        if "-" in value : #科学计数（因为也可能包含小数点。放在第一个筛选）
                            worksheet.write(item, i, float(value), style=styleEE)
                        elif "." in value: #小数
                            worksheet.write(item, i, float(value), style=styleFloat)
                        elif value == "0": #零值 显示为”0“ 
                            worksheet.write(item, i, int(value), style=styleALL)
                        else:  #整数
                            worksheet.write(item, i, int(value), style=styleInt)
                    else:
                        worksheet.write(item, i, value, style=styleALL)
            i+=1
        workbook.save(textfile.replace('.txt', '.xls'))


