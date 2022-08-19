###########################################################################                     Contains all the functions required for executing the pipeline
import time
import requests
import pandas as pd
import os
from Bio import SeqIO
from io import StringIO
import numpy as np
import iedb
import re

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.firefox.options import Options

###########################################################################                     Parse data from input file
def get_data_from_file(folder_name):
	file_input = folder_name + '/Input.txt'
	fp = open(file_input,'r')
	data = []
	for ln in fp:
		ln = ln.replace('\n','')
		data.append(ln)
	fp.close()
	return data

###########################################################################                     Parse pdb id and chain id from input file
def get_pdb_chain_from_file(folder_name):
	file_input = folder_name + '/Input.txt'
	fp = open(file_input,'r')
	pdb = []
	chain = []
	for ln in fp:
		ln = ln.replace('\n','')
		tabs = ln.split(' ')
		pdb.append(tabs[0])
		chain.append(tabs[1])
	return pdb, chain

###########################################################################                     Get sequence data based on the Swiss-Prot ID
def get_sequence(swissprot):
	baseurl = 'https://www.uniprot.org/uniprot/'
	url = baseurl + swissprot + '.fasta'
	response = requests.post(url)
	data =''.join(response.text)
	postString = data.split("\n",1)[1]
	postString = postString.replace('\n','')
	return postString

###########################################################################                     Process fasta format as a single-line sequence and get the Swiss-Prot ids from the fasta format itself
def process_fasta(folder_name):
	data = []
	swissprot = []
	li = []
	filename = folder_name + '/Input.txt'
	fp = open(filename,'r')
	flag = 0
	for ln in fp:
		ln = ln.replace('\n','')
		if ('>' in ln) & (flag == 0):
			tabs = ln.split('|')
			swissprot.append(tabs[1])
		elif ('>' not in ln) & (flag == 1):
			li.append(ln)
		elif ('>' not in ln) & (flag == 0):
			li.append(ln)
			flag = 1
		elif ('>' in ln) & (flag == 1):
			tabs = ln.split('|')
			swissprot.append(tabs[1])
			string = ''.join(li)
			data.append(string)
			li.clear()
			flag = 0
	string = ''.join(li)
	data.append(string)
	# print(swissprot)
	# print(data)
	return swissprot, data

###########################################################################                     Run Bebipred with input as Swiss-Prot ids
def run_bebipred_swissprot():
	data = get_data_from_file('Bebipred')
	detailed_output_file = 'Bebipred/Detailed_'
	output_file = 'Bebipred/Output_'
	fronturl = 'curl --data "method=Bepipred&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Bebipred for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		position = []
		residue = []
		l = 0
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & (tabs[3] == 'E'):
				flag = 1
				l = l + 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == 'E'):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == '.'):
				flag = 0
				string = ''.join(residue)
				fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
				position.clear()
				residue.clear()
		fp_write.close()
		print('\n\n-------------------Bebipred executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Bebipred with input as protein sequences
def run_bebipred_sequence():
	swissprot, data = process_fasta('Bebipred')
	detailed_output_file = 'Bebipred/Detailed_'
	output_file = 'Bebipred/Output_'
	fronturl = 'curl --data "method=Bepipred&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Bebipred for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		position = []
		residue = []
		l = 0
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & (tabs[3] == 'E'):
				flag = 1
				l = l + 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == 'E'):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == '.'):
				flag = 0
				string = ''.join(residue)
				fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
				position.clear()
				residue.clear()
		fp_write.close()
		print('\n\n-------------------Bebipred executed successfully for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		count = count + 1

###########################################################################                     Run Chou & Fasman with input as Swiss-Prot ids
def run_choufasman_swissprot():
	data = get_data_from_file('ChouFasman')
	detailed_output_file = 'ChouFasman/Detailed_'
	output_file = 'ChouFasman/Output_'
	fronturl = 'curl --data "method=Chou-Fasman&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Chou-Fasman for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Chou-Fasman executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Chou & Fasman with input as protein sequences
def run_choufasman_sequence():
	swissprot, data = process_fasta('ChouFasman')
	detailed_output_file = 'ChouFasman/Detailed_'
	output_file = 'ChouFasman/Output_'
	fronturl = 'curl --data "method=Chou-Fasman&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Chou-Fasman for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Chou-Fasman executed successfully for ' + str(swissprot[i]) + '-----------------\n\n')

###########################################################################                     Run Emini with input as Swiss-Prot ids
def run_emini_swissprot():
	data = get_data_from_file('Emini')
	detailed_output_file = 'Emini/Detailed_'
	output_file = 'Emini/Output_'
	fronturl = 'curl --data "method=Emini&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Emini for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >= 1): # Threshold is always 1.0
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= 1):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < 1):
				if (len(position) >= 6): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Emini executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Emini with input as protein sequences
def run_emini_sequence():
	swissprot, data = process_fasta('Emini')	
	detailed_output_file = 'Emini/Detailed_'
	output_file = 'Emini/Output_'
	fronturl = 'curl --data "method=Emini&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in range(0,len(data)):
		print('\n\n-------------------Running Emini for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >= 1):
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= 1):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < 1):
				if (len(position) >= 6): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Emini executed successfully for Protein ' + str(swissprot[i]) + '-----------------\n\n')

###########################################################################                     Run Karplus & Schulz with input as Swiss-Prot ids
def run_karplusschulz_swissprot():
	data = get_data_from_file('KarplusSchulz')
	detailed_output_file = 'KarplusSchulz/Detailed_'
	output_file = 'KarplusSchulz/Output_'
	fronturl = 'curl --data "method=Karplus-Schulz&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Karplus-Schulz for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Karplus-Schulz executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Karplus-Schulz with input as protein sequences
def run_karplusschulz_sequence():
	swissprot, data = process_fasta('KarplusSchulz')
	detailed_output_file = 'KarplusSchulz/Detailed_'
	output_file = 'KarplusSchulz/Output_'
	fronturl = 'curl --data "method=Karplus-Schulz&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Karplus-Schulz for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Karplus-Schulz executed successfully for ' + str(swissprot[i]) + '-----------------\n\n')

###########################################################################                     Run Kolaskar-Tongaonkar with input as Swiss-Prot ids
def run_kolaskartongaonkar_swissprot():
	data = get_data_from_file('KolaskarTongaonkar')
	detailed_output_file = 'KolaskarTongaonkar/Detailed_'
	output_file = 'KolaskarTongaonkar/Output_'
	fronturl = 'curl --data "method=Kolaskar-Tongaonkar&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Kolaskar-Tongaonkar for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Kolaskar-Tongaonkar executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Kolaskar-Tongaonkar with input as protein sequences
def run_kolaskartongaonkar_sequence():
	swissprot, data = process_fasta('KolaskarTongaonkar')
	detailed_output_file = 'KolaskarTongaonkar/Detailed_'
	output_file = 'KolaskarTongaonkar/Output_'
	fronturl = 'curl --data "method=Kolaskar-Tongaonkar&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Kolaskar-Tongaonkar for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Kolaskar-Tongaonkar executed successfully for ' + str(swissprot[i]) + '-----------------\n\n')

###########################################################################                     Run Parker with input as Swiss-Prot ids
def run_parker_swissprot():
	data = get_data_from_file('Parker')
	detailed_output_file = 'Parker/Detailed_'
	output_file = 'Parker/Output_'
	fronturl = 'curl --data "method=Parker&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in data:
		print('\n\n-------------------Running Parker for ' + i + '-----------------\n\n')
		sequence = get_sequence(i)
		url = fronturl + sequence + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + i + '.tsv'
		output_file_2 = output_file + i + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Parker executed successfully for ' + i + '-----------------\n\n')

###########################################################################                     Run Parker with input as protein sequences
def run_parker_sequence():
	swissprot, data = process_fasta('Parker')
	detailed_output_file = 'Parker/Detailed_'
	output_file = 'Parker/Output_'
	fronturl = 'curl --data "method=Parker&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Parker for Protein ' + str(swissprot[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot[i] + '.tsv'
		output_file_2 = output_file + swissprot[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Parker executed successfully for ' + str(swissprot[i]) + '-----------------\n\n')

###########################################################################                     Run Discotope-1.1 with input as protein structures (pdb ids)
def run_discotope1():
	pdb, chain = get_pdb_chain_from_file('Discotope-1.1')
	options = Options()
	options.headless = True
	discotope_url = 'http://tools.iedb.org/discotope/'
	detailed_output_file = 'Discotope-1.1/Detailed_'
	output_file = 'Discotope-1.1/Output_'
	for j in range(0,len(pdb)):
		discotope_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues']
		discotope_epitopes = [discotope_columns]
		print("\n\n------------------- Running Discotope-1.1 for Protein " + str(pdb[j]) + " -----------------\n\n")
		driver = webdriver.Firefox(options=options, executable_path = './geckodriver')
		driver.get(discotope_url)
		driver.find_element(By.ID, "id_pdb").send_keys(pdb[j])
		driver.find_element(By.ID, "id_chain").send_keys(chain[j])
		driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[3]/td[2]/select/option[1]").click() # Option-1 is for Discotope-1.1
		driver.find_element(By.NAME, "submit").click()
		wait = WebDriverWait(driver, 600)
		wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form")))
		driver.find_element(By.XPATH, "/html/body/div[3]/form/a[1]/button").click()
		wait.until(ec.visibility_of_element_located((By.ID, "result_table")))
		discotope_table = driver.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
		discotope_table_columns = driver.find_element(By.XPATH, '/html/body/div[3]/table/thead').text.splitlines()
		driver.close()
		# print(discotope_table)
		output_file_1 = detailed_output_file + pdb[j] + '_' + chain[j] + '.tsv'
		output_file_2 = output_file + pdb[j] + '_' + chain[j] + '.tsv'
		fp_write = open(output_file_2,'w')
		fp_write_1 = open(output_file_1,'w')
		fp_write.write('')
		fp_write_1.write('Chain id' + '\t' + 'Residue id' + '\t' + 'Residue name' + '\t' + 'Contact number' + '\t' + 'Propensity score' + '\t' + 'Discotope score' + '\n')
		for k in range(0,len(discotope_table)):
			tabs = discotope_table[k].split(' ')
			fp_write_1.write(str(tabs[0]) + '\t' + str(tabs[1]) + '\t' + str(tabs[2]) + '\t' + str(tabs[3]) + '\t' + str(tabs[4]) + '\t' + str(tabs[5]) + '\n')
			# print(tabs)
		fp_write_1.close()
		discotope_epitopes_split_rows = []
		for row in discotope_table:
			split_row = row.split(" ")
			discotope_epitopes_split_rows.append(split_row)
		df = pd.DataFrame(discotope_epitopes_split_rows, columns=discotope_table_columns)
		df["Discotope Score"] = pd.to_numeric(df["Discotope Score"])
		df = df.loc[df['Discotope Score'] >= -7.7] # Threshold for Discotope-1.1
		if df.empty == True:
			continue
		else:
			residue_id = df.iloc[:, 1].to_list()
			residue_name = df.iloc[:, 2].to_list()
		discontinous_epitope = []
		for i in range(len(residue_name)):
			if residue_name[i] == 'ARG':
				ARG = 'R'
				pos = ARG + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'ASN':
				ASN = 'N'
				pos = ASN + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'ASP':
				ASP = 'D'
				pos = ASP + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'GLN':
				GLN = 'Q'
				pos = GLN + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'GLU':
				GLU = 'E'
				pos = GLU + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'LYS':
				LYS = 'K'
				pos = LYS + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'PHE':
				PHE = 'F'
				pos = PHE + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'TRP':
				TRP = 'W'
				pos = TRP + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'TYR':
				TYR = 'Y'
				pos = TYR + str(residue_id[i])
				discontinous_epitope.append(pos)
			else:
				pos = residue_name[i][0] + str(residue_id[i])
				discontinous_epitope.append(pos)
		epitope = ','.join(discontinous_epitope)
		start_pos = residue_id[0]
		end_pos = residue_id[-1]
		row_to_append = [pdb[j]] + [chain[j]] + [epitope] + [start_pos] + [end_pos] + [len(residue_name)]
		discotope_epitopes.append(row_to_append)
		# print(discotope_epitopes)
		# print('\n')
		tabs = discotope_epitopes[1]
		values = tabs[2].split(',')
		# print(values)
		for k in range(0,len(values)):
			fp_write.write(str(values[k]) + '\n')
		fp_write.close()
		print("\n\n------------------- Discotope-1.1 executed successfully for Protein " + str(pdb[j]) + " -----------------\n\n")
		discotope_epitopes.clear()

###########################################################################                     Run Discotope-2.0 with input as protein structures (pdb ids)
def run_discotope2():
	pdb, chain = get_pdb_chain_from_file('Discotope-2.0')
	options = Options()
	options.headless = True
	discotope_url = 'http://tools.iedb.org/discotope/'
	detailed_output_file = 'Discotope-2.0/Detailed_'
	output_file = 'Discotope-2.0/Output_'
	for j in range(0,len(pdb)):
		discotope_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues']
		discotope_epitopes = [discotope_columns]
		print("\n\n------------------- Running Discotope-2.0 for Protein " + str(pdb[j]) + " -----------------\n\n")
		driver = webdriver.Firefox(options=options, executable_path = './geckodriver')
		driver.get(discotope_url)
		driver.find_element(By.ID, "id_pdb").send_keys(pdb[j])
		driver.find_element(By.ID, "id_chain").send_keys(chain[j])
		driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[3]/td[2]/select/option[2]").click() # Option-2 is for Discotope-2.0
		driver.find_element(By.NAME, "submit").click()
		wait = WebDriverWait(driver, 600)
		wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div[3]/form")))
		driver.find_element(By.XPATH, "/html/body/div[3]/form/a[1]/button").click()
		wait.until(ec.visibility_of_element_located((By.ID, "result_table")))
		discotope_table = driver.find_element(By.XPATH, '/html/body/div[3]/table/tbody').text.splitlines()
		discotope_table_columns = driver.find_element(By.XPATH, '/html/body/div[3]/table/thead').text.splitlines()
		driver.close()
		# print(discotope_table)
		output_file_1 = detailed_output_file + pdb[j] + '_' + chain[j] + '.tsv'
		output_file_2 = output_file + pdb[j] + '_' + chain[j] + '.tsv'
		fp_write = open(output_file_2,'w')
		fp_write_1 = open(output_file_1,'w')
		fp_write.write('')
		fp_write_1.write('Chain id' + '\t' + 'Residue id' + '\t' + 'Residue name' + '\t' + 'Contact number' + '\t' + 'Propensity score' + '\t' + 'Discotope score' + '\n')
		for k in range(0,len(discotope_table)):
			tabs = discotope_table[k].split(' ')
			fp_write_1.write(str(tabs[0]) + '\t' + str(tabs[1]) + '\t' + str(tabs[2]) + '\t' + str(tabs[3]) + '\t' + str(tabs[4]) + '\t' + str(tabs[5]) + '\n')
			# print(tabs)
		fp_write_1.close()
		discotope_epitopes_split_rows = []
		for row in discotope_table:
			split_row = row.split(" ")
			discotope_epitopes_split_rows.append(split_row)
		df = pd.DataFrame(discotope_epitopes_split_rows, columns=discotope_table_columns)
		df["Discotope Score"] = pd.to_numeric(df["Discotope Score"])
		df = df.loc[df['Discotope Score'] >= -3.7] # Threshold for Discotope-2.0
		if df.empty == True:
			continue
		else:
			residue_id = df.iloc[:, 1].to_list()
			residue_name = df.iloc[:, 2].to_list()
		discontinous_epitope = []
		for i in range(len(residue_name)):
			if residue_name[i] == 'ARG':
				ARG = 'R'
				pos = ARG + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'ASN':
				ASN = 'N'
				pos = ASN + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'ASP':
				ASP = 'D'
				pos = ASP + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'GLN':
				GLN = 'Q'
				pos = GLN + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'GLU':
				GLU = 'E'
				pos = GLU + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'LYS':
				LYS = 'K'
				pos = LYS + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'PHE':
				PHE = 'F'
				pos = PHE + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'TRP':
				TRP = 'W'
				pos = TRP + str(residue_id[i])
				discontinous_epitope.append(pos)
			elif residue_name[i] == 'TYR':
				TYR = 'Y'
				pos = TYR + str(residue_id[i])
				discontinous_epitope.append(pos)
			else:
				pos = residue_name[i][0] + str(residue_id[i])
				discontinous_epitope.append(pos)
		epitope = ','.join(discontinous_epitope)
		start_pos = residue_id[0]
		end_pos = residue_id[-1]
		row_to_append = [pdb[j]] + [chain[j]] + [epitope] + [start_pos] + [end_pos] + [len(residue_name)]
		discotope_epitopes.append(row_to_append)
		# print(discotope_epitopes)
		# print('\n')
		tabs = discotope_epitopes[1]
		values = tabs[2].split(',')
		# print(values)
		for k in range(0,len(values)):
			fp_write.write(str(values[k]) + '\n')
		fp_write.close()
		print("\n\n------------------- Discotope-2.0 executed successfully for Protein " + str(pdb[j]) + " -----------------\n\n")
		discotope_epitopes.clear()

###########################################################################                     Run Ellipro with input as protein structures (pdb ids)
def run_ellipro():
	pdb = get_data_from_file('Ellipro')
	options = Options()
	options.headless = True
	ellipro_url = 'http://tools.iedb.org/ellipro/'
	detailed_output_file_le = 'Ellipro/Linear_Epitope_Detailed_'
	output_file_le = 'Ellipro/Linear_Epitope_Output_'
	detailed_output_file_de = 'Ellipro/Discontinuous_Epitope_Detailed_'
	output_file_de = 'Ellipro/Discontinuous_Epitope_Output_'
	for m in range(0,len(pdb)):
		linear_columns = ['pdb_id','chain','start','end','peptide','nr_of_residues','score']
		discontinous_columns = ['pdb_id','chain','peptide','start','end','nr_of_residues','score']
		linear_epitopes = [linear_columns]
		discontinous_epitopes = [discontinous_columns]
		print("\n\nPredicting linear and discontinous epitopes of protein " + pdb[m] + " with Ellipro\n\n")
		driver = webdriver.Firefox(options=options, executable_path = './geckodriver')
		driver.get(ellipro_url)
		driver.find_element(By.NAME, "pdb_id").send_keys(pdb[m])
		driver.find_element(By.NAME, "submit").click()
		wait = WebDriverWait(driver, 180)
		wait.until(ec.visibility_of_element_located((By.ID, "result_table")))
		driver.find_element(By.NAME, "chain").click()
		chain = driver.find_element(By.XPATH, "/html/body/div[3]/form/table/tbody/tr[1]/td[3]").text
		driver.find_element(By.NAME, "submit").click()
		wait.until(ec.visibility_of_element_located((By.CLASS_NAME, "output_title")))
		table_of_linear_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[2]/tbody').text.splitlines()
		columns_of_linear_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[2]/thead').text.splitlines()
		table_of_discontinous_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/tbody').text.splitlines()
		columns_of_discontinous_epitopes = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/thead').text.splitlines()
		linear_epitopes_split_rows = []
		# print(table_of_linear_epitopes)
		# print(table_of_discontinous_epitopes)
		for row in table_of_linear_epitopes:
			split_row = row.split(" ")
			linear_epitopes_split_rows.append(split_row)
		df = pd.DataFrame(linear_epitopes_split_rows, columns=columns_of_linear_epitopes[:-1])
		df["Score"] = pd.to_numeric(df["Score"])
		df = df.loc[df['Score'] >= 0.5]
		if df.empty == True:
			continue
		else:
			rows = [[i for i in row[1:]] for row in df.itertuples()]
			for i in rows:
				i = [pdb[m]] + i[1:]
				linear_epitopes.append(i)
		rows = []
		for i in range(len(table_of_discontinous_epitopes)):
			row = []
			for e in range(len(columns_of_discontinous_epitopes[:-1])):
				cell = driver.find_element(By.XPATH, '/html/body/div[3]/table[3]/tbody/tr[' + str(i+1) + ']/td[' + str(e+1) + ']').text
				row.append(cell)
			rows.append(row)
		driver.close()
		for arow in rows:
			if float(arow[3]) >= 0.5:
				unproc_aa_list = arow[1].split(', ')
				aa_list = []
				for position in unproc_aa_list:
					aa = position[2:]
					aa_list.append(aa)
				epitope = ','.join(aa_list)
				start_pos = unproc_aa_list[0][3:]
				end_pos = unproc_aa_list[-1][3:]
				nr_of_residues = len(aa_list)
				score = float(arow[3])
				row_to_append = [pdb[m]] + [chain] + [epitope] + [start_pos] + [end_pos] + [nr_of_residues] + [score]
				discontinous_epitopes.append(row_to_append)
		# print(linear_epitopes)
		# print('\n\n')
		# print(discontinous_epitopes)
		fp_write = open(detailed_output_file_le,'w')
		for k in range(0, len(linear_epitopes)):
			tabs = linear_epitopes[k]
			for l in range(0, len(tabs)):
				fp_write.write(str(tabs[l]) + '\t')
			fp_write.write('\n')
		fp_write.close()
		fp_write = open(detailed_output_file_de,'w')
		for k in range(0, len(discontinous_epitopes)):
			tabs = discontinous_epitopes[k]
			for l in range(0, len(tabs)):
				fp_write.write(str(tabs[l]) + '\t')
			fp_write.write('\n')
		fp_write.close()
		print("\n\n------------------- Ellipro executed successfully for Protein " + str(pdb[m]) + " -----------------\n\n")
		linear_epitopes.clear()
		discontinous_epitopes.clear()

###########################################################################                     Run MSA with input as UniProt accession numbers
def run_msa(matrix='bl62', gapopen=1.53, gapext=0.123, order='aligned', nbtree=2, treeout='true', maxiterate=2, ffts='none'): 
	data = get_data_from_file('MSA')
	print(data)
	seq_string = ''
	for i in data:
		seq_string = seq_string + 'sp:' + i + ','
	seq_string = seq_string[:-1]
	print(seq_string)
	os.system(
		'python3.9 msa_algos/mafft.py --email barnali.das@tum.de --stype protein --sequence ' + seq_string + ' --outfile MSA/MAFFT/mafft --format fasta --matrix ' + matrix + ' --gapopen ' + str(gapopen) + ' --gapext ' + str(gapext) 
		+ ' --order ' + order + ' --nbtree ' + str(nbtree) + ' --treeout ' + treeout + ' --maxiterate ' + str(maxiterate) + ' --ffts ' + ffts)
	os.system(
		'python3.9 msa_algos/muscle.py --email barnali.das@tum.de --sequence ' + seq_string + ' --outfile MSA/MUSCLE/muscle --format fasta')

###########################################################################                     Run GBlocks to extract the conserved sequences from the MSA
def run_gblocks(path_to_file, folder, min_seq_conserved_pos='default', min_seq_flank_pos='default', max_contigous_nonconserved_pos = 8, min_length_block= 10, allowed_gap_pos='None'):
	options = Options()
	options.headless = True
	gblocks_url = 'https://ngphylogeny.fr/tools/tool/276/form'
	driver = webdriver.Firefox(options=options, executable_path = './geckodriver')
	driver.maximize_window()
	driver.get(gblocks_url)
	upload_file = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[1]/div/input")
	upload_file.send_keys(path_to_file)
	data_type = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[2]/div[1]/div/select")
	data_type.send_keys('Protein')
	b1 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[3]/div/input")
	b1.clear()
	b1.send_keys(min_seq_conserved_pos)
	b2 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[4]/div/input")
	b2.clear()
	b2.send_keys(min_seq_flank_pos)
	b3 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[5]/div/input")
	b3.clear()
	b3.send_keys(max_contigous_nonconserved_pos)
	b4 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[6]/div/input")
	b4.clear()
	b4.send_keys(min_length_block)
	b5 = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/fieldset/div/div/div[7]/div/select")
	b5.send_keys(allowed_gap_pos)
	submit = driver.find_element(By.XPATH, "/html/body/div/div[1]/div[1]/div/div/form/input")
	submit.click()
	wait = WebDriverWait(driver, 180)
	wait.until(ec.visibility_of_element_located((By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]")))
	sequence_info_html = driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[2]/td[5]/a[5]")
	sequence_info_html.click()
	wait.until(ec.visibility_of_element_located(
		(By.XPATH, "/html/body/div/div[1]/div[2]/pre/pre[3]/b[1]")))
	positions = driver.find_element(By.XPATH, '/html/body/div/div[1]/div[2]/pre/pre[3]').get_attribute("innerHTML").splitlines()[1]
	positions = positions[8:]
	list_positions = [int(s) for s in re.findall(r'\d+', positions)]
	driver.execute_script("window.history.go(-1)")
	wait.until(ec.visibility_of_element_located(
		(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]")))
	cleaned_sequences_fasta = driver.find_element(By.XPATH, "/html/body/div/div[1]/div/div/div[6]/div/table/tbody/tr[3]/td[4]/a[5]")
	cleaned_sequences_fasta.click()
	wait.until(ec.visibility_of_element_located(
		(By.XPATH, "/html/body/div/div[1]/div[2]/pre")))
	cleaned_seq = driver.find_element(By.XPATH, '/html/body/div/div[1]/div[2]/pre').text
	cleaned_seq_list = cleaned_seq.split('>')
	cleaned_seq_list.pop(0)
	conserved_sequences_dictionary = {}
	for i in range(len(cleaned_seq_list)):
		lines = cleaned_seq_list[i].splitlines()
		protein_id = lines[0].split(' ')[1]
		cleaned_fasta_all = ''.join(lines[1:])
		list_of_fasta = []
		for i in range(1, len(list_positions), 2):
			cut_point = list_positions[i] - list_positions[i - 1] + 1
			cons_seq = cleaned_fasta_all[:cut_point]
			list_of_fasta.append(cons_seq)
			cleaned_fasta_all = cleaned_fasta_all[cut_point:]
		conserved_sequences_dictionary[protein_id] = list_of_fasta
	driver.close()
	file_write = folder + 'gblocks_conserved_sequences.txt'
	fp = open(file_write,'w')
	for key, value in conserved_sequences_dictionary.items(): 
		fp.write('%s:%s\n' % (key, value))
	print(conserved_sequences_dictionary)
	fp.close()

###########################################################################                     Run Bebipred with input as conserved subsequences of MSA
def run_bebipred_msa(swissprot, data, folder):
	detailed_output_file = folder + 'Bebipred/Detailed_'
	output_file = folder + 'Bebipred/Output_'
	fronturl = 'curl --data "method=Bepipred&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Bebipred for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		position = []
		residue = []
		l = 0
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & (tabs[3] == 'E'):
				flag = 1
				l = l + 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == 'E'):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & (tabs[3] == '.'):
				flag = 0
				string = ''.join(residue)
				fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
				position.clear()
				residue.clear()
		time.sleep(3) # Bebipred fails sometimes for batch job submission
		fp_write.close()
		print('\n\n-------------------Bebipred executed successfully for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		count = count + 1

###########################################################################                     Run Chou & Fasman with input as conserved subsequences of MSA
def run_choufasman_msa(swissprot, data, folder):
	detailed_output_file = folder + 'ChouFasman/Detailed_'
	output_file = folder + 'ChouFasman/Output_'
	fronturl = 'curl --data "method=Chou-Fasman&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Chou-Fasman for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Chou-Fasman executed successfully for ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')

###########################################################################                     Run Emini with input as conserved subsequences of MSA
def run_emini_msa(swissprot, data, folder):
	detailed_output_file = folder + 'Emini/Detailed_'
	output_file = folder + 'Emini/Output_'
	fronturl = 'curl --data "method=Emini&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	for i in range(0,len(data)):
		print('\n\n-------------------Running Emini for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		# print(result)
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >= 1):
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= 1):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < 1):
				if (len(position) >= 6): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Emini executed successfully for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')

###########################################################################                     Run Parker with input as conserved subsequences of MSA
def run_parker_msa(swissprot, data, folder):
	detailed_output_file = folder + 'Parker/Detailed_'
	output_file = folder + 'Parker/Output_'
	fronturl = 'curl --data "method=Parker&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Parker for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Parker executed successfully for ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')

###########################################################################                     Run Karplus-Schulz with input as conserved subsequences of MSA
def run_karplusschulz_msa(swissprot, data, folder):
	detailed_output_file = folder + 'KarplusSchulz/Detailed_'
	output_file = folder + 'KarplusSchulz/Output_'
	fronturl = 'curl --data "method=Karplus-Schulz&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Karplus-Schulz for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Karplus-Schulz executed successfully for ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')

###########################################################################                     Run Kolaskar-Tongaonkar with input as conserved subsequences of MSA
def run_kolaskartongaonkar_msa(swissprot, data, folder):
	detailed_output_file = folder + 'KolaskarTongaonkar/Detailed_'
	output_file = folder + 'KolaskarTongaonkar/Output_'
	fronturl = 'curl --data "method=Kolaskar-Tongaonkar&sequence_text='
	backurl = '" http://tools-cluster-interface.iedb.org/tools_api/bcell/'
	count = 1
	for i in range(0,len(data)):
		print('\n\n-------------------Running Kolaskar-Tongaonkar for Protein ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')
		url = fronturl + data[i] + backurl
		# print(url)
		result = os.popen(url).read()
		output_file_1 = detailed_output_file + swissprot + '_' + data[i] + '.tsv'
		output_file_2 = output_file + swissprot + '_' + data[i] + '.tsv'
		fp_write_1 = open(output_file_1,'w')
		for line in result.splitlines():
			# print(line)
			fp_write_1.write(line)
			fp_write_1.write('\n')
		fp_write_1.close()
		fp_write = open(output_file_2,'w')
		fp_write.write('Sl. No' + '\t' + 'Start position' + '\t' + 'End position' + '\t' + 'Epitope sequence' + '\t' + 'Length' + '\t' + '\n')
		with open(output_file_1) as fp:
			lines = fp.readlines()
		del lines[0]
		total = 0
		count = 0
		for j in range(0,len(lines)):
			count = count + 1
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			total = total + float(tabs[5])
		threshold = total/count
		# print(threshold)
		position = []
		residue = []
		l = 1
		flag = 0
		for j in range(0, len(lines)):
			ln = lines[j].replace('\n','')
			tabs = ln.split('\t')
			if (flag == 0) & ((float(tabs[5])) >=  threshold): 
				flag = 1
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) >= threshold):
				position.append(tabs[0])
				residue.append(tabs[1])
			elif (flag == 1) & ((float(tabs[5])) < threshold):
				if (len(position) >= 7): # default window size
					flag = 0
					string = ''.join(residue)
					fp_write.write(str(l) + '\t' + str(position[0]) + '\t' + str(position[-1]) + '\t' + string + '\t' + str(len(position)) + '\n')
					l = l + 1
					position.clear()
					residue.clear()
				else:
					flag = 0
					position.clear()
					residue.clear()
		fp_write.close()
		print('\n\n-------------------Kolaskar-Tongaonkar executed successfully for ' + str(swissprot) + ' for conserved subsequence-' + str(data[i]) + '-----------------\n\n')

###########################################################################                     Reading the conserved sequences from file and predicting the linear epitopes from them
def run_conserved_sequences(file, folder):
	fp_read = open(file,'r')
	for ln in fp_read:
		ln = ln.replace(':',',')
		ln = ln.replace('[','')
		ln = ln.replace(']','')
		ln = ln.replace('\n','')
		ln = ln.replace('\'','')
		tabs = ln.split(',')
		swissprot = tabs[0]
		tabs.pop(0)
		tabs = [x.strip(' ') for x in tabs]
		print(swissprot)
		print(tabs)
		run_bebipred_msa(swissprot,tabs,folder)
		run_choufasman_msa(swissprot,tabs,folder)
		run_emini_msa(swissprot,tabs,folder)
		run_parker_msa(swissprot,tabs,folder)
		run_karplusschulz_msa(swissprot,tabs,folder)
		run_kolaskartongaonkar_msa(swissprot,tabs,folder)
		tabs.clear()
	fp_read.close()