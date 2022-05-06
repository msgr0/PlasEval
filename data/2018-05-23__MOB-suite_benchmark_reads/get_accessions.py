from sys import argv
import os

samples_file = "samples.csv"

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

accessions_file = "accessions.txt"
accessions = open(accessions_file,"w")

samples = read_file(samples_file)

for line in samples[1:]:
	acc_nos = line.split(";")[-1].split(",")
	for acc_no in acc_nos:
		accessions.write(acc_no+"\n")	