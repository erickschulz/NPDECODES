#!/usr/bin/env python
# coding=utf-8
# Sets up an empty file structure for a NumPDE homework problem
# Invoke with: new_problem <ProblemName> <problemname>
# author: Liaowang Huang

import os
import re
import sys

# to be replaced
Old_name = 'NewProblem'
old_name = 'newproblem'


def replace_string(string, Replace_name, replace_name):
	"""
	replace every occurrence of Old_name and old_name in string by Replace_name and replace_name respectively

	:param string: string needs replacement
	:param Replace_name:
	:param replace_name:
	"""
	string = re.sub(Old_name, Replace_name, string)
	string = re.sub(old_name, replace_name, string)
	return string


def copy_replace(input_name, outdir, Replace_name, replace_name):
	"""
	Recursively copy the file or directory input to outdir and do the replacement

	:param input_name : a file or a directory
	:param outdir: the output directory
	:param Replace_name: replace all the OldName by Replace_name
	:param replace_name: replace all the oldname by replace_name
	"""

	if not outdir.endswith("/"):
		outdir += "/"

	if os.path.isfile(input_name):
		# new file name after replacement
		filename = input_name.split('/')[-1]
		outputname = replace_string(filename, Replace_name, replace_name)

		# replace in the file content
		file_handle = open(input_name, 'r')
		file_string = file_handle.read()
		file_handle.close()

		file_string = replace_string(file_string, Replace_name, replace_name)

		# write to the new file
		new_file = open(outdir + outputname, "w")
		new_file.write(file_string)
	else:
		# recursively copy the directory
		if not input_name.endswith("/"):
			input_name += "/"

		dir_name = input_name.split("/")[-2]
		dir_output = replace_string(dir_name, Replace_name, replace_name)

		if not os.path.exists(outdir + dir_output):
			os.mkdir(outdir + dir_output)

		for fileOrFolders in os.listdir(input_name):
			if fileOrFolders.startswith('.'):  # ignore hidden files
				continue
			copy_replace(input_name + fileOrFolders, outdir + dir_output, Replace_name, replace_name)


if __name__ == '__main__':
	assert len(sys.argv) == 3, "There should be two arguments."
	Replace_name = sys.argv[1]  # e.g. BurgersEquation
	replace_name = sys.argv[2]  # e.g. burgersequation

	this_dir = os.path.dirname(os.path.realpath(__file__))  # ".../scripts"
	os.chdir(this_dir)  # change the current working directory to this_dir

	develop_path = "../developers/"

	copy_replace("./NewProblem", develop_path, Replace_name, replace_name)
