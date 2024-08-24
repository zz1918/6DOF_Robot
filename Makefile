# Makefile for Delta robot
#
#  HOW TO USE THIS Makefile:
#
# (A)	To compile run this demo, you can use this Makefile and just have to type:
#
# 	       > make 		-- compiles the main program into "main"
# 	       > make help	-- usage help
# 	       > make test	-- runs "main" non-interactively
# 	       			   (this is Core Library's automatic self-test)
#
# (B)	But you probably want to try the various interactive demos of this program:
#
# 	       > make r	-- runs "main" interactively
#
# 	The default search-mode is "width".  We encourage you to just
# 	hit the "Run me" button to see different runs of the program.
#
#	This is equivalent to typing:
#
#		./build/Debug/main \
#			$(interactive) $(Qtype) $(epsilon) $(envrange) $(edgeL) \
#			$(startOx) $(startOy) $(startOz) \
#			$(startAphi) $(startAtheta) $(startBtheta) \
#			$(goalOx) $(goalOy) $(goalOz) \
#			$(goalAphi) $(goalAtheta) $(startBtheta) \
#
#	where the arguments to main correspond to the values $(var) of
#	various Makefile variables.
#   The variables $(Aphi) $(Atheta) $(Btheta) defines triangle AOB as  
#	point A = O + edgeL * (cos(Aphi) , cos(Atheta) * sin(Aphi), sin(Atheta) * sin(Aphi)) and
#	point B = O + edgeL * ( -sin(Aphi) * cos(Btheta), 
#						  	cos(Atheta) * cos(Aphi) * cos(Btheta) - sin(Atheta) * sin(Btheta), 
#							sin(Atheta) * cos(Aphi) * cos(Btheta) + cos(Atheta) * sin(Btheta))
#
# (C)	You can override any of these values at the command line.
#	For instance, to change the length of the right edge of the robot
#	to 55, and to change the start position of point O of the robot to (20,12,0), you can type:
#
#		> make eg edgeL=55 startOx=20 startOy=12 startOz=0
#
# (D)	Instead of the target "eg" you can use any of these targets:
#
#		ega, egb, egc
#		eg1,				-- no path example
#		eg2, eg2a			-- like input2.txt, but with bounding box
#		eg3, eg3a, eg3b			-- 100 random triangles
#		eg4				-- no path example
#
# 	These targets uses different input files, all taken from
# 	the subdirectory "inputs".
#
#
# August 19th, 2024
# --Zhaoqi Zhang, Chee Yap, Yijen Chiang
#


#=================================================
# User variables (you can change them in the command line)
#=================================================
interactive = 0		# -1=interactive, >=0 is non-interactive

startOx = -2.5		# start configuration
startOy = -2.5
startOz = -2.5
startAphi = 0	# set A to be O + (1,0,0)
startAtheta = 0
startBtheta = 0		# set B to be O + (0,1,0)
goalOx = 2.5		# goal configuration
goalOy = 2.5
goalOz = 2.5
goalAphi = 0	# set A to be O + (1,0,0)
goalAtheta = 0
goalBtheta = 0		# set B to be O + (0,1,0)

epsilon = 0.05		# resolution parameter
envrange = 4		# environment boundary
edgeL = 1.0		# robot radius

Qtype = "width"
ExpandLimit = 40000

file1 = "wall1"
file2 = "wall2"

#=================================================
# Define target folder
#=================================================

# Default is to initialize, compile and run the program with default arguments.
d default: r
# Initial build up the program from cmake.
b build:
	mkdir build
	(cd build; cmake ../main; cmake --build .)
	./build/Debug/main \
		$(interactive) $(Qtype) $(epsilon) $(envrange) $(edgeL) \
		$(startOx) $(startOy) $(startOz) \
		$(startAphi) $(startAtheta) $(startBtheta) \
		$(goalOx) $(goalOy) $(goalOz) \
		$(goalAphi) $(goalAtheta) $(startBtheta) \
		$(ExpandLimit) $(file1) $(file2)
# Initialize the program from cmake
i initialize:
	rm -rf build
	mkdir build
	(cd build; cmake ../main)
# Compile the program from cmake
c compile: 
	(cd build; cmake --build .)
# Run program by default settings
r run:
	./build/Debug/main \
		$(interactive) $(Qtype) $(epsilon) $(envrange) $(edgeL) \
		$(startOx) $(startOy) $(startOz) \
		$(startAphi) $(startAtheta) $(startBtheta) \
		$(goalOx) $(goalOy) $(goalOz) \
		$(goalAphi) $(goalAtheta) $(goalBtheta) \
		$(ExpandLimit) $(file1) $(file2)
# Eliminate the program
e eliminate:
	rm -rf build
# Test 1
t1 test1:
	make r startAphi=0.785398 startAtheta=0 startBtheta=1.570796
# note: this target is the standard target that Core Library uses
#       to test its subdirectories.   So the program must run in a
#       non-interactive mode (i.e., the first argument to "main" is "1").
#

help:
	@echo "USAGE:  The main program is called main.  Demos can be invoked thus:"
	@echo "         > make eg"
	@echo "         > make egX"
	@echo "     where X is replaced by 0,1,2,3,4,5,6,100,200,300, etc."
	@echo "     Some demos have variants, such as egXa, egXb or egXc:"
	@echo "     LIST OF TARGETS:"
	@echo "	ega, egb, egc"
	@echo "	eg0,				-- no path example"
	@echo "	eg1, eg1a, (bug, buga)		-- bugtrap"
	@echo "	eg2, eg2a			-- like input2.txt, but with bounding box"
	@echo "	eg3, eg3a, eg3b			-- 100 random triangles"
	@echo "	eg4				-- maze"
	@echo "	eg5, eg5a (bug2, bug2a)		-- double bugtrap"
	@echo "	eg6, eg6a 			-- example from Kavraki's OOPSMP"
	@echo "	eg100, eg100a, eg100b		-- 100 random triangles"
	@echo "	eg200, eg200a			-- 200 random triangles"
	@echo "	eg300				-- 300 random triangles"



#=================================================
# Processing environment
#=================================================
f=input11

process:
	cp inputs/$(f).txt output-tmp.txt
	python input_interpreter.py

#=================================================
# Temp
#=================================================
m:
	gvim Makefile
vi:
	gvim main.cpp

#=================================================
# Rules
#=================================================
%: %.o
	${CXX} $(OBJ_FILES) $(LDFLAGS) $< $(CORE_LIB) -o $@

.cpp.o:
	${CXX} -c -O3 $(CXXFLAGS) $(CORE_INC) $< -o $@

#=================================================
# Clean object files
#=================================================
clean:
	-@test -z "*.o" || rm -f *.o

#=================================================
# Remove executable files
#=================================================
EXEPROGS=$(TARGETS:=$(EXETYPE))

vclean veryclean: clean
	-@test -z "$(EXEPROGS)" || rm -f $(EXEPROGS)

#=================================================
# END Makefile
#=================================================