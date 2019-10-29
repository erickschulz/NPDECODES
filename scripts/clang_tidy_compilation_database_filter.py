# file: clang_tidy_compilation_database_filter.py
# brief: removes files from the compilation database that shall not be checked by clang tidy (i.e. homework templates)
# author: Till Ehrengruber <tille@student.ethz.ch>
# date: 18.03.2019

import os
import sys
import json
import re

prefix=os.path.dirname(os.path.abspath(__file__))+"/../build"
# a set of patterns for files that shall be removed from the database
skip_patterns=["templates/", "template/", "mysolution/"]

with open(prefix + "/compile_commands.json", "r") as read_file:
    data = json.load(read_file)
    new_data = []
    for entry in data:
      if not any(map(lambda p: re.search(p, entry["file"]), skip_patterns)):
      	new_data.append(entry)

    with open(prefix + "/compile_commands.json", "w") as write_file:
      json.dump(new_data, write_file, indent=4, separators=(',', ': '))
