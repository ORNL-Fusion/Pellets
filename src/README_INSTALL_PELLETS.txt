
1. Install Xcode, open the program and let it install “components”

2. Open terminal, type: xcode-select --install    You might get a warning that these are already installed

3.Install brew: https://brew.sh

4. From the terminal install these programs:
	brew install NetCDF
	brew install pig-config
	brew install cmake

5. Install gfortran from: https://github.com/fxcoudert/gfortran-for-macOS/releases. You can probably also install this with brew: brew install gcc. You can always check what the current command is in brew by typing: brew search gfortran

6. Check to see if you have gfortran under /usr/local/

7. Open terminal, use cd to get to the pellets/src directory. 

8. Make sure that you have the CMakeLists.txt already in there. You'll need to make sure to include the directory which has 'netcdf.mod' If you used brew, it should be under /usr/local/include. This should be in your CMakeLists.txt file as: include_directories(/usr/local/include)

9. Type cmake ./

10. Type make -f MakeFile

Done!
-Jeff

