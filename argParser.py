#
# ArgParser.py
#
# by Victor Hanson-Smith
# victorhansonsmith@gmail.com
#
# This is a super useful class to include in your Python project that
# needs many command-line options passed to the program.
#
# Usage in your code:
#
# 1. Instantiate an instance of this class:
# ap = ArgParser(sys.argv)
#
# 2. Get the values of required arguments:
# x = ap.getArg("--inputfile")
# --> where "--inputfile" can be really anything.
#
# 3. Get the values of optional arguments:
# x = ap.getOptionalArg("--verbose")
# --> getOptionalArg will return False if this argument was not specified by the user.
#
#


import sys
import os
import re

from debugTools import *

class ArgParser:
	def __init__(self, cliArgs):
		self.args = cliArgs
		self.params = {}
	
	# use this method to grab REQUIRED command-line parameters:
	def getArg(self, flag):
		if self.args.__contains__(flag):
			i = self.args.index(flag)
			return self.args[i+1]
		else:
			printError( argDoesNotExist )
			printTip("argument: " + flag)
			exit(1)
			
	def getOptionalArg(self, flag):
		if self.args.__contains__(flag):
			i = self.args.index(flag)
			return self.args[i+1]
		else:
			return False
	
	# This method will return a list of tokens following 'flag', but not including
	# tokens which start with "--"
	def getList(self, flag):
		if self.args.__contains__(flag):
			i = self.args.index(flag)
			returnList = []
			flagPattern = re.compile("^\-\-.*")
			for j in range( i+1, self.args.__len__() ):
				if re.match(flagPattern, self.args[j] ):
					return returnList
				else:
					returnList.append( self.args[j] )
			return returnList
		else:
			printError( argDoesNotExist )
			printTip("argument: " + flag)
			exit(1)		
		
	
	# use this method to grab OPTIONAL command-line toggles (boolean on/off switches)
	def getOptionalToggle(self, flag):
		if self.args.__contains__(flag):
			return True
		else:
			return False
	
	def doesContainArg(self, flag):
		return self.args.__contains__(flag)

	def getOptionalList(self, flag, type=str):
		#print self.args, flag
		if self.args.__contains__(flag):
			i = self.args.index(flag)
			returnList = []
			flagPattern = re.compile("^\-\-.*")
			for j in range( i+1, self.args.__len__() ):
				if re.match(flagPattern, self.args[j] ):
					return returnList
				else:
					returnList.append( type(self.args[j]) )
			return returnList
		else:
			return None
	