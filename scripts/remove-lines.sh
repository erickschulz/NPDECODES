#! /bin/bash

usage() {
	echo;
  echo "  Usage:";
  echo "  $0 filename";
  echo "   -- will remove text between '/* SOLUTION_BEGIN */' & '/* SOLUTION_END */'"
  echo "      from file and replace it with:"
  echo "      '/* TODO Your implementation goes here! */'";
  echo;
  echo "  $0 STR1 STR2 filename";
  echo "   -- will remove text between STR1 & STR2 from file";
  echo;
  echo "  $0 STR1 STR2 filename REPL";
  echo "   -- will remove text between STR1 & STR2 from file and replace";
  echo "      it with REPL";
  echo;
  echo "  Suggested workflow for creating templates and 'mysolution' from";
  echo "  mastersolution:";
  echo "  1) Mark the solution parts in your code with";
  echo "     '/* SOLUTION_BEGIN */' mastersolution code '/* SOLUTION_END */'";
  echo "  2) Create the directory structure:"
  echo "     homeworks";
  echo "     └── YourProblem";
  echo "         ├── mastersolution";
  echo "         │   └── main.cc <- your solution goes here";
  echo "         ├── mysolution";
  echo "         └── templates";
  echo " 3) In the directory 'YourProblem' do:";
  echo "    $ ../../scripts/remove-lines.sh mastersolution/main.cc > mysolution/main.cc";
  echo "    $ cp mysolution/* templates/";
  echo "    The remove-lines script automatically replaces the content in";
  echo "    between the flags '/* SOLUTION_BEGIN(END) */' with";
  echo "    '/* TODO Your implementation goes here! */'."
  echo "    For other options see Usage above.";
  echo;
	exit;
}

# check parameters
if [ $# == 1 ]
then
  FILE=$1;
  STR1="/* SOLUTION_BEGIN */";
  STR2="/* SOLUTION_END */";
  REPL="/* TODO Your implementation goes here! */";
  # fix regex characters (put escape \ infront)
  STR1=$(echo "$STR1" | sed 's/\*/\\\\*/g');
  STR2=$(echo "$STR2" | sed 's/\*/\\\\*/g');
elif [ $# == 3 ]
then 
  STR1=$(echo "$1" | sed 's/\*/\\\\*/g');
  STR2=$(echo "$2" | sed 's/\*/\\\\*/g');
  FILE=$3
  REPL=""
elif [ $# == 4 ]
then 
  STR1=$(echo "$1" | sed 's/\*/\\\\*/g');
  STR2=$(echo "$2" | sed 's/\*/\\\\*/g');
  FILE=$3
  REPL=$4
else
  usage
fi

# check if file exists
[ ! -f "$FILE" ] && echo "File '$FILE' not found..." && exit;

# main - using awk
# Purpose: Convert lines of format 
# ```   ... something ...
#       STR1 
#       ... something ...
#       STR2                
#       ... something ...   ```
# Into
# ```   ... something ...
#       STR1 
#       (REPL, if it's been set)
#       STR2              
#       ... something ...   ```
# In addition leading whitespaces before STR2 are preserved and 
# also added before REPL, to ensure nice formatting. 
awk -v STR1="$STR1" -v STR2="$STR2" -v REPL="$REPL" '
BEGIN{
     	found1=0; # false
}
{
     	line=$0;
	skipLine=0; # false

	if (! found1 && match(line, STR1".*"STR2)){
		sub(STR1".*"STR2, STR1" "STR2, line);
		if (line=="") skipLine=1; # true
	} else if (! found1 && match(line, STR1)) {
		#sub(STR1".*", STR1"", line);
		found1=1; # true
		if (line=="") skipLine=1;
	} else if (found1 && match(line, STR2)) {
    	# if text before STR2 is all whitespaces 
    	# keep them for nice formatting
    	split(line, a, STR2);
    	if (a[1] ~ "^ *$") {
      		leading_ws = a[1];
      		# There is nothing before STR2 to be 
      		# replaced, so no call to `sub`
		# sub(".*"STR2, ""STR2, line);
    	}
    	else { 
      		leading_ws = "";
      		# remove whatever is before STR2
		sub(".*"STR2, ""STR2, line); 
    	}
	found1=0; # false
    	# print replacement text with leading 
    	# whitespaces, if they exist
    	if (REPL != "") print leading_ws REPL
		if (line=="") skipLine=1; # true
	} else if (found1==1) {
		skipLine=1; # true
		line="";
	}

	if (!skipLine) print line;

}' "$FILE" | cat -s
