import os
import argparse

#USAGE: python format_input_v2.py --i input_dir/putative_plasmid_contigs.fasta --o output_dir --f formatted_file.tsv

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

if __name__ == "__main__":
	#Parsing arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--i", help="Path to input (HyAsP output) file")
	parser.add_argument("--o", help="Path to formatted file output directory")
	parser.add_argument("--f", help="Name of formatted file")

	args = parser.parse_args()	

	prediction_file = args.i
	output_dir = args.o
	formatted_file = args.f

	pred_details = []
	c_id, p_id, c_len = '', '', 0

	string_list = read_file(prediction_file)
	for line in string_list:
		if line[0] == '>':
			c_id = line[1:].split("|")[0]
			p_id = line[1:].split("|")[1].split("_")[2]
		else:
			c_len = len(line)
			pred_details.append(('plasmid_'+p_id, c_id, str(c_len)))
			#print(c_id, p_id, str(c_len))

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)	
	out_file = open(os.path.join(output_dir, formatted_file), "w")				

	out_file.write("#Plasmid\tContig_ID\tLength\n")
	for x in pred_details:
		out_file.write(x[0]+"\t"+x[1]+"\t"+x[2]+"\n")

